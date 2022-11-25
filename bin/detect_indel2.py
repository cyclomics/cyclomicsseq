#!/usr/bin/env python
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pysam
from vcf_tools import VCF_entry


@dataclass
class Indel:
    found: bool = False
    depth: int = 0
    depth_q: int = 0
    type: str = ""
    length: int = 0
    assemby: str = ""
    start_position: int = 0
    reference_nucleotide: str = ""
    variant_nucleotide: str = ""
    support_count: int = 0
    support_ratio: float = 0.0
    support_count_fwd: int = 0
    support_ratio_fwd: float = 0.0
    support_count_rev: int = 0
    support_ratio_rev: float = 0.0
    base_qualities: list = field(default_factory=list)
    alt_base_qualities: list = field(default_factory=list)


def check_indel(
    pu_column,
    variant_count_th=10,
):
    """Return most common indel above thresholds for a pileup column

    Checks pileupcolumn for indel evidence from
    pysam.PileupColumn.get_query_sequences,
    then chooses most common insert or deletion to call as a variant.
    """

    # TODO: Check homopolymer status
    # Initialize general statistics of this pileup position
    total = pu_column.nsegments
    total_q = pu_column.get_num_aligned()

    query_qualities = pu_column.get_query_qualities()

    variant = Indel(None)
    variant.depth = total
    variant.depth_q = total_q
    variant.assemby = pu_column.reference_name
    variant.start_position = pu_column.pos

    # Get positional sequences, from which to search for indels
    pu = pu_column.get_query_sequences(
        mark_matches=True, mark_ends=True, add_indels=True
    )

    quality_map = list(zip(pu, query_qualities))

    # The reference nucleotide should be the most common in this position
    variant.reference_nucleotide = Counter(pu).most_common(1)[0][0].upper()[0]

    # Exclude matches, SNPs or irrelevant characters
    # The remainder should be indels
    excluded_chars = [".", ",", "*", "^", "A", "C", "T", "G", "]", "$"]
    # pu_indel = [x for x in quality_map if x[0].upper() not in excluded_chars]
    pu_indel = [
        x for x in quality_map if (char not in x[0].upper() for char in excluded_chars)
    ]

    inserts = [(re.split(r"(\d+)", x[0])[-1], x[1]) for x in pu_indel if "+" in x[0]]
    # remove the directional indicator and make uppercase
    # eg convert A+3TTT' into TTT and A+1a to A
    inserts_fwd = [x[0].upper() for x in inserts if x[0][-1].isupper()]
    inserts_rev = [x[0].upper() for x in inserts if x[0][-1].islower()]

    deletions = [(re.split(r"(\d+)", x[0])[-1], x[1]) for x in pu_indel if "-" in x[0]]
    # remove the directional indicator and make uppercase
    # eg G-15NNNNNNNNNNNNNNN g-15nnnnnnnnnnnnnnn turn into
    # 15NN... and 15NN.. in their respective variable
    deletions_fwd = [x[0].upper() for x in deletions if x[0][-1] == "N"]
    deletions_rev = [x[0].upper() for x in deletions if x[0][-1] == "n"]

    # If no insert or deletion evidence was found in this position,
    # the 'inserts' list will be empty and the Counter will return
    # an IndexError, which can be understood as 'no evidence'.
    try:
        # new insert needs to be at least 50% of inserts_fwd
        new_insert_fwd = Counter(inserts_fwd).most_common(1)[0][0]
        new_insert_rev = Counter(inserts_rev).most_common(1)[0][0]
        new_insert_fwd_count = Counter(inserts_fwd).most_common(1)[0][1]
        new_insert_rev_count = Counter(inserts_rev).most_common(1)[0][1]
        new_insert_count = new_insert_fwd_count + new_insert_rev_count
        new_insert_orientation_ratio = min(
            new_insert_fwd_count / new_insert_rev_count,
            new_insert_rev_count / new_insert_fwd_count,
        )

    except IndexError:
        new_insert_count = 0
        new_insert_orientation_ratio = 0
        new_insert_fwd, new_insert_rev = "", ""

    try:
        new_deletion_fwd = Counter(deletions_fwd).most_common(1)[0][0]
        new_deletion_rev = Counter(deletions_rev).most_common(1)[0][0]
        new_deletion_fwd_count = Counter(deletions_fwd).most_common(1)[0][1]
        new_deletion_rev_count = Counter(deletions_rev).most_common(1)[0][1]
        new_deletion_count = new_deletion_fwd_count + new_deletion_rev_count
        new_deletion_orientation_ratio = min(
            new_deletion_fwd_count / new_deletion_rev_count,
            new_deletion_rev_count / new_deletion_fwd_count,
        )

    except IndexError:
        new_deletion_count = 0
        new_deletion_orientation_ratio = 0
        new_deletion_fwd, new_deletion_rev = "", ""

    # Decide if we have found either an Insert or a Deletion
    # TODO: inserts and deletions in separate functions, run both
    insert_indicator = False
    deletion_indicator = False
    try:
        if new_insert_count == 0 and new_deletion_count == 0:
            # No variant was found
            return variant

        elif (
            new_insert_fwd == new_insert_rev
            and new_insert_count > new_deletion_count
            and new_insert_fwd_count > variant_count_th
            and new_insert_rev_count > variant_count_th
        ):
            # Found an insertion
            insert_indicator = True

        elif (
            new_deletion_fwd == new_deletion_rev
            and new_deletion_count > new_insert_count
            and new_deletion_fwd_count > variant_count_th
            and new_deletion_rev_count > variant_count_th
        ):
            # Found a deletion
            deletion_indicator = True

        else:
            # No variant was found.
            # This include the case that there is indel evidence,
            # but it either didn't pass the threshold or there was
            # evidence found for both and indel and an insert with the
            # same number of supporting reads.
            return variant

    except IndexError:
        pass

    # Output variant found information and statistics
    if insert_indicator:
        variant.found = True
        variant.type = "insertion"
        variant.variant_nucleotide = new_insert_fwd
        variant.length = len(new_insert_fwd)
        variant.support_count = new_insert_count
        variant.support_ratio = new_insert_count / total_q
        variant.support_count_fwd = new_insert_fwd_count
        variant.support_ratio_fwd = new_insert_fwd_count / total_q
        variant.support_count_rev = new_insert_rev_count
        variant.support_ratio_rev = new_insert_rev_count / total_q
        variant.base_qualities = quality_map
        variant.alt_base_qualities = inserts

    elif deletion_indicator:
        variant.found = True
        variant.type = "deletion"
        variant.variant_nucleotide = new_deletion_fwd
        variant.length = len(new_deletion_fwd)
        variant.support_count = new_deletion_count
        variant.support_ratio = new_deletion_count / total_q
        variant.support_count_fwd = new_deletion_fwd_count
        variant.support_ratio_fwd = new_deletion_fwd_count / total_q
        variant.support_count_rev = new_deletion_rev_count
        variant.support_ratio_rev = new_deletion_rev_count / total_q
        variant.base_qualities = quality_map
        variant.alt_base_qualities = deletions

    return variant


def extract_indel_evidence(
    pileupcolumn: pysam.PileupColumn,
    assembly: str,
    reference: pysam.FastaFile,
    pos: int,
    high_base_quality_cutoff=80,
    end_of_amplicon=False,
):
    """Returns found variants in a given pileup position"""

    if not end_of_amplicon:
        indel = check_indel(pileupcolumn)
        if indel.found:
            if indel.type == "deletion":
                # Adjust reference allele to include deleted seq
                ref_seq = reference.fetch(
                    reference=assembly, start=pos, end=pos + indel.length + 1
                )
                alleles = (ref_seq, ref_seq[0])
            else:
                alleles = (
                    indel.reference_nucleotide,
                    str(indel.reference_nucleotide + indel.variant_nucleotide),
                )

            # Update variant statistics with the found indel
            indel_type = indel.type
            depth = indel.depth
            depth_q = indel.depth_q
            total_count = indel.support_count
            total_ratio = indel.support_ratio
            fwd_count = indel.support_count_fwd
            fwd_ratio = indel.support_ratio_fwd
            rev_count = indel.support_count_rev
            rev_ratio = indel.support_ratio_rev

            abq = np.mean(
                [
                    x[1]
                    for x in indel.alt_base_qualities
                    if indel.variant_nucleotide in x[0]
                ]
            )

            obq = np.mean([x[1] for x in indel.base_qualities])

            high_quality_bases = [
                x for x in indel.base_qualities if x[1] > high_base_quality_cutoff
            ]
            high_quality_alt_bases = [
                x
                for x in indel.alt_base_qualities
                if x[1] > high_base_quality_cutoff and x[0] == indel.variant_nucleotide
            ]
            hcr = len(high_quality_alt_bases) / len(high_quality_bases)

    # Return statistics for this position,
    # whether a variant was found or not
    vcf_entry = VCF_entry(
        DP=depth,
        DPQ=depth_q,
        FREQ=total_ratio,
        VAF=total_ratio,
        FWDC=fwd_count,
        FWDR=fwd_ratio,
        REVC=rev_count,
        REVR=rev_ratio,
        TOTC=total_count,
        TOTR=total_ratio,
        SAME=1,
        OBSR=total_ratio,
        ABQ=abq,
        OBQ=obq,
        HCR=hcr,
    )

    return (assembly, alleles, vcf_entry)


if __name__ == "__main__":
    import argparse

    dev = True
    if not dev:
        parser = argparse.ArgumentParser(
            description=("Detect indels in BAM alignment file.")
        )

        parser.add_argument("fasta", type=Path)
        parser.add_argument("bed", type=Path)
        parser.add_argument("bam", type=Path)
        parser.add_argument("vcf_out", type=Path)
        args = parser.parse_args()

        extract_indel_evidence(args.bam, args.bed, args.fasta, args.vcf_out)

    if dev:
        # EGFR
        fasta = Path(
            "/data/references/Homo_sapiens/GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna"
        )
        bed = Path("pos_EGFR.bed")
        bam = Path(
            "/scratch/projects/ROD_0908_63_variantcalling/results/PR_test/consensus_aligned/FAS12641.taged.bam"
        )
        vcf_out = Path("testindel_EGFR.vcf")

        extract_indel_evidence(bam, bed, fasta, vcf_out)
