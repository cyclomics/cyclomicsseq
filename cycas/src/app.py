import gc
import json
import re
import sys
from pathlib import Path
from typing import Iterable, Iterator, List, Optional

import numpy as np
import pysam
from loguru import logger
from tqdm import tqdm

from .alignment import AlignmentGroup
from .barcode_extractor import UNKNOWN_BARCODE
from .calibrated_consensus.algorithms.fraction import FractionCalibratedConsensus
from .classification_rules import RuleFactory as MetadataRuleFactory
from .cluster_metrics import ClusterMetrics
from .config import ReadClassificationConfig
from .consensus_caller import ConsensusCallerMetadata
from .filters import Filter, FilterSecondaryAlignments
from .hdf5_parse import Hdf5Writer
from .metadata_classifier import MetadataClassifier
from .splitting_strategy import RelativeIdSplit
from .structure_classifier import StructureClassifier

# Add parent directory to Python path for local imports
# sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


# Silence all logs except ERROR and CRITICAL
logger.remove()
logger.add(sys.stderr, level="ERROR")


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


def convert_to_ascii(quality: int, offset: int = 33) -> str:
    return chr(quality + offset)


def _create_fastq_entry(read_name: str, sequence: str, quality: List[int]) -> str:
    quality_str = "".join(map(convert_to_ascii, quality))
    return f"@{read_name}\n{sequence}\n+\n{r'{}'.format(quality_str)}\n"


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


def _write_split_bam(
    output_bam_path: str,
    input_bam: pysam.AlignmentFile,
    processed_alignments: Iterable[List[pysam.AlignedSegment]],
):
    """
    Writes processed alignments to a new BAM file using the header from the input BAM.

    Args:
        output_bam_path: Path to the output BAM file.
        input_bam: Original BAM file (for header).
        processed_alignments: Iterable of lists of processed alignments.
    """
    with pysam.AlignmentFile(
        output_bam_path, "wb", header=input_bam.header.to_dict()
    ) as out_bam:
        for group in processed_alignments:
            for aln in group:
                out_bam.write(aln.to_pysam())


def _create_metadata(group, rules, empty_group=False):
    result = {}
    result["readname"] = group.read_name
    result["alignment_count"] = len(group)
    result["alignment_count_pre_filter"] = group.alignment_count_prefilter

    if empty_group:
        return result

    result["contig"] = group.alignments[0].alignment_chromosome
    result["raw_length"] = group.alignments[0].raw_read_length
    result["aligned_bases_before_consensus"] = group.count_aligned_bases()
    result["unaligned_segments"] = group.find_unaligned_regions()
    result["longest_unaligned_segment"] = group.find_longest_unaligned_region()
    result["mean_length_unaligned_segment"] = group.find_mean_unaligned_region()
    result["median_length_unaligned_segment"] = group.find_median_unaligned_region()
    result["alternative_representation"] = group.create_intermediate_read_structure()
    # put the classification rules in a seperate entity
    result["rules"] = {}
    for rule in rules:
        result["rules"][rule.name] = rule.score(group)
    return result


def _quick_process_read(
    read,
    bam_indexed,
    filters,
    rules,
    full_metadata,
    classifier,
    consensus_caller,
    fastq,
    hdf5_file,
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
    if hdf5_file:
        read_blocks_corrected = []
        for block in blocks:
            consensus_alignments = consensus_caller.create_consensus_alignment_objects(
                block.alignments
            )
            consensus_caller.update_insert_locations(consensus_alignments)
            read_blocks_corrected.append(consensus_alignments)

        # Save to HDF5
        hdf5_file.save_consensus_blocks(metadata, read_blocks_corrected)

    logger.debug(f"Finished {read}")


def create_quick_consensus(
    bam_file,
    output,
    plot_readtypes,
    plot_piechart,
    limit_calls,
    metadata_json=None,
    filters=[FilterSecondaryAlignments()],
    classifier=MetadataClassifier(),
    consensus_caller=ConsensusCallerMetadata(),
    rule_factory=MetadataRuleFactory(),
    hdf5_file=None,
):
    """
    Main function for consensus generation
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    bam_indexed = pysam.IndexedReads(bam)
    bam_indexed.build()

    # obtain all rules used for feature extraction used by classification
    rules = rule_factory.get_rules()
    hdf5_file = Hdf5Writer(hdf5_file) if hdf5_file else None
    full_metadata = {}
    read_names = _get_read_names(bam)

    total_calls = len(read_names)
    if limit_calls != 0:
        if total_calls < limit_calls:
            logger.warning(
                f"There are less reads ({total_calls}) than the requested --limit-calls ({limit_calls}). Will process all available reads."
            )
        else:
            total_calls = limit_calls
            read_names = read_names[:total_calls]

    with open(output, "w+") as fastq:
        for i, read in tqdm(enumerate(read_names), total=total_calls):
            _quick_process_read(
                read,
                bam_indexed,
                filters,
                rules,
                full_metadata,
                classifier,
                consensus_caller,
                fastq,
                hdf5_file,
            )

    # show the count of each class
    logger.info(f"{classifier.show_counts()}")

    if hdf5_file:
        hdf5_file.close()

    if metadata_json:
        with open(metadata_json, "w") as out_json:
            out_json.write(json.dumps(full_metadata))

    if plot_readtypes:
        for read_type in set(classifier.result.values()):
            print(f"plotting {read_type}")
            _plot_specific_type(classifier, bam, read_type, amount=4, filters=filters)


def create_quick_consensus_sr(
    bam_file,
    result_fastq,
    read,
    metadata_json=None,
    filters=[FilterSecondaryAlignments()],
    classifier=MetadataClassifier(),
    consensus_caller=ConsensusCallerMetadata(),
    rule_factory=MetadataRuleFactory(),
):
    bam = pysam.AlignmentFile(bam_file, "rb")
    bam_indexed = pysam.IndexedReads(bam)
    bam_indexed.build()

    fastq_file = Path(result_fastq)

    # obtain all rules used for feature extraction used by classification
    rules = rule_factory.get_rules()

    full_metadata = {}
    with open(fastq_file, "w+", encoding="utf-8") as fastq:
        _quick_process_read(
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
        with open(metadata_json, "w", encoding="utf-8") as out_json:
            out_json.write(json.dumps(full_metadata))

    for read_type in set(classifier.result.values()):
        print(f"plotting {read_type}")
        _plot_specific_type(classifier, bam, read_type, amount=4, filters=filters)


def _process_reads_classification(
    read: str,
    bam_indexed: pysam.IndexedReads,
    filters: List[Filter],
    metadata_rules: List,
    structure_classifier: StructureClassifier,
    metadata_classifier: MetadataClassifier,
    config: ReadClassificationConfig,
    make_msa: bool = False,
) -> Iterator[tuple[dict, list]]:
    """
    Processes a single 3D read by generating alignment metadata and applying classification.

    Args:
        read: Read name to process.
        bam_indexed: Indexed BAM file object.
        filters: List of filters to apply to alignments.
        metadata_rules: Rules used to generate metadata.
        structure_classifier: Classifier for structural classification.
        metadata_classifier: Classifier for metadata-based classification.
        config: Configuration object containing strategy and thresholds.

    Returns:
        An iterator of tuples containing metadata dict and list of alignments.
    """

    consensus_caller = ConsensusCallerMetadata()

    # 1. Create annotated AlignmentGroup with config
    group = AlignmentGroup(read, bam_indexed, filters, config=config)
    original_structure = group.structure
    all_blocks = group.all_blocks

    # 2. Split reads based on relative ID
    split_groups = RelativeIdSplit().split(group)

    for subgroup in split_groups:
        # 3. Create metadata
        metadata = _create_metadata(subgroup, metadata_rules)
        metadata["original_structure"] = original_structure

        # Compute cluster metrics for the first cluster
        metrics = ClusterMetrics(
            cluster_id=subgroup.clusters[0], all_blocks=all_blocks
        ).calculate_metrics()

        metadata["cluster_metrics"] = metrics

        # 4. Classify structure based on the calculated metrics
        structure_type = structure_classifier.classify(subgroup.clusters, metrics)
        metadata["structure_classification"] = structure_type

        # 5. Classify metadata
        metadata["classification"] = metadata_classifier.classify(metadata)

        if make_msa:
            consensus_alignments = consensus_caller.create_consensus_alignment_objects(
                subgroup.alignments
            )
            consensus_caller.update_insert_locations(consensus_alignments)

            for aln in subgroup.alignments:
                aln.extended_seq = consensus_alignments[aln].seq
                aln.extended_qual = consensus_alignments[aln].qual
                aln.extended_cigar = consensus_alignments[aln].cigar
                aln.direction = consensus_alignments[aln].direction

        yield (metadata, subgroup.alignments)


def create_classifications(
    bam_file: str,
    output_bam: str,
    metadata_json: str,
    limit_calls: int = 0,
    config: Optional[ReadClassificationConfig] = None,
    filters: Optional[List[Filter]] = None,
    structure_classifier: Optional[StructureClassifier] = None,
    metadata_classifier: Optional[MetadataClassifier] = None,
    metadata_rule_factory: Optional[MetadataRuleFactory] = None,
    hdf5_file=None,
):
    """
    Classifies 3D reads from a BAM file and writes the resulting metadata to a JSON file.

    Args:
        bam_file: Path to the input BAM file.
        metadata_json: Path to write the resulting metadata JSON.
        limit_calls: Maximum number of reads to classify (0 means no limit).
        config: Optional ReadClassificationConfig configuration object.
        filters: List of filters to apply (default is filtering secondary alignments).
        structure_classifier: Optional custom structure classifier.
        metadata_classifier: Optional custom metadata classifier.
        metadata_rule_factory: Optional custom metadata rule factory.
    """
    # TODO: make filters user arguments
    filters = filters or [FilterSecondaryAlignments()]

    config = config or ReadClassificationConfig()
    structure_classifier = structure_classifier or StructureClassifier(config=config)
    metadata_classifier = metadata_classifier or MetadataClassifier()
    metadata_rule_factory = metadata_rule_factory or MetadataRuleFactory()
    metadata_rules = metadata_rule_factory.get_rules()

    hdf5_file = Hdf5Writer(hdf5_file) if hdf5_file else None

    # classification_counts = Counter()
    # count_processed = 0

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        bam_indexed = pysam.IndexedReads(bam)
        bam_indexed.build()

        full_metadata = {}
        processed_alignments = []

        for i, read in enumerate(tqdm(_get_read_names(bam))):
            processed = list(
                _process_reads_classification(
                    read,
                    bam_indexed,
                    filters,
                    metadata_rules,
                    structure_classifier,
                    metadata_classifier,
                    config=config,
                    make_msa=hdf5_file is not None,
                )
            )

            # metadata_blocks = []
            for metadata, alignments in processed:
                if hdf5_file:
                    hdf5_file.save_raw_blocks(metadata, alignments)

                full_metadata[metadata["readname"]] = metadata
                # metadata_blocks.append(metadata)
                processed_alignments.append(alignments)

            # full_metadata[read] = metadata_blocks

            if limit_calls > 0 and i >= limit_calls:
                logger.warning(f"Limit of {limit_calls} reads reached. Stopping.")
                break

        # Write BAM
        _write_split_bam(output_bam, bam, processed_alignments)

    # classification_counts["total"] = count_processed
    with open(metadata_json, "w", encoding="utf-8") as out_json:
        out_json.write(json.dumps(full_metadata))

    if hdf5_file:
        hdf5_file.close()

    # Clean up
    del bam_indexed
    gc.collect()


def _process_reads_consensus(
    read: str,
    bam_indexed: pysam.IndexedReads,
    filters: List[Filter],
    metadata_rules: List,
    structure_classifier: StructureClassifier,
    metadata_classifier: MetadataClassifier,
    consensus_caller: ConsensusCallerMetadata,
    calibrator: FractionCalibratedConsensus,
    config: ReadClassificationConfig,
) -> Iterator[dict]:
    """
    Processes a single 3D read by generating alignment metadata and applying classification.

    Args:
        read: Read name to process.
        bam_indexed: Indexed BAM file object.
        filters: List of filters to apply to alignments.
        metadata_rules: Rules used to generate metadata.
        structure_classifier: Classifier for structural classification.
        metadata_classifier: Classifier for metadata-based classification.

    Returns:
        A dictionary containing classification and metadata for the given read.
    """

    consensus_caller = ConsensusCallerMetadata()

    # 1. Create annotated AlignmentGroup with config
    group = AlignmentGroup(read, bam_indexed, filters, config=config)
    original_structure = group.structure
    all_blocks = group.all_blocks

    # 2. Split reads based on relative ID
    split_groups = RelativeIdSplit().split(group)

    for subgroup in split_groups:
        # 3. Create metadata
        metadata = _create_metadata(subgroup, metadata_rules)
        metadata["original_structure"] = original_structure

        # Compute cluster metrics for the first cluster
        metrics = ClusterMetrics(
            cluster_id=subgroup.clusters[0], all_blocks=all_blocks
        ).calculate_metrics()

        metadata["cluster_metrics"] = metrics

        # 4. Classify structure based on the calculated metrics
        structure_type = structure_classifier.classify(subgroup.clusters, metrics)
        metadata["structure_classification"] = structure_type

        if structure_type in ["too_few_inserts", "no_block_classification"]:
            yield (metadata, None, group.high_reject_ratio)
            continue

        # 5. Classify metadata
        metadata["classification"] = metadata_classifier.classify(metadata)

        # 6. Make MSA
        consensus_alignments = consensus_caller.create_consensus_alignment_objects(
            subgroup.alignments
        )
        consensus_caller.update_insert_locations(consensus_alignments)

        data = {
            "msa_seq": np.array(
                [list(aln.seq) for aln in consensus_alignments.values()]
            ),
            "msa_qual": np.array([aln.qual for aln in consensus_alignments.values()]),
            "msa_cigar": np.array(
                [list(aln.cigar) for aln in consensus_alignments.values()]
            ),
            "direction": [aln.direction for aln in consensus_alignments.values()],
            "metadata": metadata,
        }

        # TODO: this will copy the contents of data into itself...
        data["consensus_result"] = calibrator._make_consensus(data)
        data = calibrator.calibrate(data)

        consensus_read = _create_fastq_entry(
            f"{data['metadata']['readname']}_{data['metadata']['structure_classification']}",
            data["consensus_sequence"],
            data["consensus_quality_calibrated"],
        )

        yield (metadata, consensus_read, group.high_reject_ratio)


def call_consensus(
    bam_file: str,
    output_fastq: str,
    metadata_json: str,
    limit_calls: int = 0,
    config: Optional[ReadClassificationConfig] = None,
    filters: Optional[List[Filter]] = None,
    structure_classifier: Optional[StructureClassifier] = None,
    metadata_classifier: Optional[MetadataClassifier] = None,
    metadata_rule_factory: Optional[MetadataRuleFactory] = None,
    consensus_caller: Optional[ConsensusCallerMetadata] = None,
    calibrator: Optional[FractionCalibratedConsensus] = None,
):
    filters = filters or [FilterSecondaryAlignments()]

    config = config or ReadClassificationConfig()
    structure_classifier = structure_classifier or StructureClassifier(config=config)

    metadata_classifier = metadata_classifier or MetadataClassifier()
    metadata_rule_factory = metadata_rule_factory or MetadataRuleFactory()
    consensus_caller = consensus_caller or ConsensusCallerMetadata()

    calibrator = calibrator or FractionCalibratedConsensus()
    calibrator.load(config.calibration_model)

    metadata_rules = metadata_rule_factory.get_rules()
    full_metadata = {}

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        bam_indexed = pysam.IndexedReads(bam)
        bam_indexed.build()
        n_processed_reads = 0
        n_high_reject_reads = 0
        n_consensus_reads = 0

        with open(output_fastq, "w+") as fastq:
            for i, read in enumerate(tqdm(_get_read_names(bam))):
                reject_warnings = []
                processed = list(
                    _process_reads_consensus(
                        read,
                        bam_indexed,
                        filters,
                        metadata_rules,
                        structure_classifier,
                        metadata_classifier,
                        consensus_caller,
                        calibrator,
                        config=config,
                    )
                )

                n_processed_reads += 1

                for (
                    metadata,
                    consensus_read,
                    high_reject_warning,
                ) in processed:
                    reject_warnings.append(high_reject_warning)
                    full_metadata[metadata["readname"]] = metadata
                    if consensus_read:
                        n_consensus_reads += 1
                        fastq.write(consensus_read)

                if any(reject_warnings):
                    n_high_reject_reads += 1

                if limit_calls > 0 and i >= limit_calls:
                    logger.warning(f"Limit of {limit_calls} reads reached. Stopping.")
                    break

    logger.info(f"{n_processed_reads} reads processed.")
    logger.info(
        f"{n_high_reject_reads} reads had a high alignment rejection ratio (> 1.8)."
    )
    logger.info(f"{n_consensus_reads} consensus reads generated.")

    with open(metadata_json, "w", encoding="utf-8") as out_json:
        out_json.write(json.dumps(full_metadata))

    # Clean up
    del bam_indexed
    gc.collect()


def simulate_reads(*args, **kwargs):
    from src.read_simulation import simulate_reads

    return simulate_reads(*args, **kwargs)


def calibration_file_prep(*args, **kwargs):
    from src.calibrated_consensus.default import file_prep

    return file_prep(*args, **kwargs)


def calibration_fit_model(*args, **kwargs):
    from src.calibrated_consensus.default import fit_model

    return fit_model(*args, **kwargs)
