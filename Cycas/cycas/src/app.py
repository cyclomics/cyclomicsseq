from cgitb import reset
from curses import meta
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Dict

import click
import pysam
from loguru import logger

from src.alignment import AlignmentGroup
from src.barcode_extractor import UNKNOWN_BARCODE
from src.consensus_caller import ConsensusCallerMetadata
from src.filters import (
    FilterMinimunMappingQuality,
    FilterMinimunRawReadLength,
    FilterSecondaryAlignments,
)
from src.plots.piechart import create_piechart
from tqdm import tqdm

from src.classification_rules import RuleFactory
from src.metadata_classifier import MetadataClassifier


def _get_read_names(bam: pysam.AlignmentFile) -> str:
    processed = set()
    i = 0
    for read in bam.fetch():
        read_name = read.query_name
        if read_name not in processed:
            i += 1
            # yeilding does not get all read names for some reason
            # yield read_name
            processed.add(read_name)
    return processed


def _create_fastq_entry(read_name, sequence, quality) -> str:
    return f"@{read_name}\n{sequence}\n+\n{quality}\n"


def _backbone_first_sorter(x):
    # alignment_chromosomes_present returns a set
    try:
        y = x.alignment_chromosomes_present.pop()
    # pop raises keyerror for pop on empty list
    except KeyError:
        return False
    return y[0].startswith("BB")


def _plot_specific_type(classifier, bam, read_type: str, amount: int = 5, filters=[]):
    bam_indexed = pysam.IndexedReads(bam)
    bam_indexed.build()

    available_reads = []
    for r_name, r_type in classifier.result.items():
        if r_type == read_type:
            available_reads.append(r_name)

    for i in range(min(len(available_reads), amount)):
        selected_read = available_reads[i]
        logger.info(f"plotting {i}: {selected_read}")

        try:
            group = AlignmentGroup(
                selected_read, bam_indexed, [FilterSecondaryAlignments()]
            )
        except TypeError:
            logger.critical(
                f"failed to create unfiltered group object for read {selected_read}"
            )
            continue

        try:
            group.create_spagettogram(title=f"{read_type}:{selected_read}")
        except TypeError:
            logger.critical(f"failed to plot unfiltered for read {selected_read}")
            continue

        logger.info(f"class: {classifier.result[selected_read]}")

        try:
            for k, v in classifier.result_details[selected_read].items():
                logger.info(f"{k}: {v}")
        except KeyError:
            logger.info("No details available")
        except AttributeError:
            logger.info("No details available on object.")

        group = AlignmentGroup(selected_read, bam_indexed, filters)
        group.create_spagettogram(title=f"{read_type}:{selected_read} Filtered")


def create_Y_tags(metadata, sequence):
    """
    Convert the metadata to follow the Y Bam tag specification
    """
    barcode = "_".join(
        [
            x["barcode"]
            for x in metadata["consensus_reads"].values()
            if x["barcode"] != UNKNOWN_BARCODE
        ]
    )

    consensus = sequence
    classification = metadata["classification"]
    original_readlen = metadata["raw_length"]
    partner_information = [""]
    ancestor_mapping = "|".join(
        [x["alignment_position"] for x in metadata["consensus_reads"].values()]
    )
    alignment_count = metadata["alignment_count"]

    Y_bam_tags = {
        "YB": barcode,
        "YL": len(consensus),
        "YC": classification,
        "YT": original_readlen,
        "YP": partner_information,
        "YA": ancestor_mapping,
        "YM": alignment_count,
    }
    return Y_bam_tags


def _create_fastq_reads(metadata) -> str:
    """
    BB_0_BBCR
    I_0_TP53
    """
    reads = []
    consensus_data = metadata["consensus_reads"]
    for k, v in consensus_data.items():
        if v["filtered"]:
            continue
        # extract the int from the string
        tag_id = re.findall(r"\d+", k)[0]
        if k.startswith("backbone"):
            assembly_type = "BB"
        else:
            assembly_type = "I"
        unique_readname = (
            f"{metadata['readname']}_{assembly_type}_{tag_id}_{v['alignment_position']}"
        )
        v["readname_unique"] = unique_readname
        reads.append(_create_fastq_entry(unique_readname, v["seq"], v["qual"]))

    return reads


def _create_metadata(group, rules, empty_group=False):
    result = {}
    result["readname"] = group.read_name
    result["alignment_count"] = len(group)
    result["alignment_count_pre_filter"] = group.alignment_count_prefilter

    if empty_group:
        return result

    result["raw_length"] = group.alignments[0].raw_read_length
    result["aligned_bases_before_consensus"] = group.count_aligned_bases()
    result["unaligned_segments"] = group.find_unaligned_regions()
    result["longest_unaligned_segment"] = group.find_longest_unaligned_region()
    result["mean_length_unaligned_segment"] = group.find_mean_unaligned_region()
    result["median_length_unaligned_segment"] = group.find_median_unaligned_region()
    result["alternative representation"] = group.create_intermediate_read_structure()
    # put the classification rules in a seperate entity
    result["rules"] = {}
    for rule in rules:
        result["rules"][rule.name] = rule.score(group)
    return result


def _process_read(
    read,
    bam_indexed,
    filters,
    rules,
    full_metadata,
    classifier,
    consensus_caller,
    fastq,
):
    """
    Code for each seperate read
    """
    logger.debug(f"readname: {read}")
    group = AlignmentGroup(read, bam_indexed, filters)

    logger.debug(f"found alignments #: {len(group)}")
    if len(group) == 0:
        metadata = _create_metadata(
            AlignmentGroup(read, bam_indexed, []), rules, empty_group=True
        )
        full_metadata[read] = metadata
        metadata["classification"] = "Filtered"
    # Create the metadata, contains the outcome of all the rules for this group
    metadata = _create_metadata(group, rules)
    # group._print_group_properties(4)

    # classify based on the rule outcome
    metadata["classification"] = classifier.classify(metadata)
    logger.debug(f"{read}: {metadata['classification']}")

    # feed it back to itself to separate the alignments correctly
    blocks = classifier.separate_alignments(metadata["classification"], group)

    # Feed the separated data to the consensus caller
    metadata["consensus_reads"] = consensus_caller.call(blocks)

    # write to the outputs
    fastq_reads = _create_fastq_reads(metadata)

    # process results
    fastq.write("\n".join(fastq_reads))
    full_metadata[read] = metadata
    logger.debug(f"Finished {read}")


def create_classifications_consensus(
    bam_file,
    output,
    plot_readtypes,
    plot_piechart,
    limit_calls,
    metadata_json=None,
    filters=[FilterSecondaryAlignments()],
    classifier=MetadataClassifier(),
    consensus_caller=ConsensusCallerMetadata(),
    rule_factory=RuleFactory(),
):
    """
    Main function for consensus generation
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    bam_indexed = pysam.IndexedReads(bam)
    bam_indexed.build()

    fastq_file = Path(output)

    # obtain all rules used for feature extraction used by classification
    rules = rule_factory.get_rules()

    full_metadata = {}
    with open(output, "w+") as fastq:
        for i, read in tqdm(enumerate(_get_read_names(bam))):
            _process_read(
                read,
                bam_indexed,
                filters,
                rules,
                full_metadata,
                classifier,
                consensus_caller,
                fastq,
            )

            # stop if we have done enough calls, usefull for debugging and comparisons
            if limit_calls != 0 and i > limit_calls:
                err_msg = f"Stopping to limit runtime, caused by --limit-calls {limit_calls}, set to zero to disable"
                logger.critical(err_msg)
                break

    # show the count of each class
    logger.info(f"{classifier.show_counts()}")

    if metadata_json:
        with open(metadata_json, "w") as out_json:
            out_json.write(json.dumps(full_metadata))

    if plot_readtypes:
        for read_type in set(classifier.result.values()):
            print(f"plotting {read_type}")
            _plot_specific_type(classifier, bam, read_type, amount=4, filters=filters)


def classify_single_read(
    bam_file,
    result_fastq,
    read,
    metadata_json=None,
    filters=[FilterSecondaryAlignments()],
    classifier=MetadataClassifier(),
    consensus_caller=ConsensusCallerMetadata(),
    rule_factory=RuleFactory(),
):
    bam = pysam.AlignmentFile(bam_file, "rb")
    bam_indexed = pysam.IndexedReads(bam)
    bam_indexed.build()

    fastq_file = Path(result_fastq)

    # obtain all rules used for feature extraction used by classification
    rules = rule_factory.get_rules()

    full_metadata = {}
    with open(fastq_file, "w+") as fastq:
        _process_read(
            read,
            bam_indexed,
            filters,
            rules,
            full_metadata,
            classifier,
            consensus_caller,
            fastq,
        )

    if metadata_json:
        with open(metadata_json, "w") as out_json:
            out_json.write(json.dumps(full_metadata))

    for read_type in set(classifier.result.values()):
        print(f"plotting {read_type}")
        _plot_specific_type(classifier, bam, read_type, amount=4, filters=filters)
