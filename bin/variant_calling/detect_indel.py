#!/usr/bin/env python
import re
import time
from collections import Counter
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pysam
from variant_calling.vcf_tools import VCF_entry, create_bed_lines, create_ref_mapper


@dataclass
class Indel:
    """
    Object to store Indel evidence information.
    """

    found: bool = False
    depth: int = 0
    depth_q: int = 0
    type: str = ""
    length: int = 0
    assemby: str = ""
    start_position: int = 0
    # reference_nucleotide: str = ""
    variant_nucleotide: str = ""
    support_count: int = 0
    support_ratio: float = 0.0
    support_count_fwd: int = 0
    support_ratio_fwd: float = 0.0
    support_count_rev: int = 0
    support_ratio_rev: float = 0.0
    base_qualities: list = field(default_factory=list)
    alt_base_qualities: list = field(default_factory=list)


def check_indel(pu_column: pysam.PileupColumn, variant_count_th: int = 10) -> Indel:
    """
    Choose most common Indel evidence in a given position.

    Args:
        pu_column: Pysam pileup column at a given position.
        variant_count_th: Minimum count of an indel for it to be considered.
            (Default = 10)


    Returns:
        Indel object encoding information at this position.
    """

    # TODO: Check homopolymer status
    # Initialize general statistics of this pileup position
    total = pu_column.nsegments
    total_q = pu_column.get_num_aligned()

    fwd_count = 0
    rev_count = 0
    for pu_read in pu_column.pileups:
        if pu_read.alignment.is_forward:
            fwd_count += 1
        elif pu_read.alignment.is_reverse:
            rev_count += 1

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
    # variant.reference_nucleotide = Counter(pu).most_common(1)[0][0].upper()

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

    except IndexError:
        new_insert_count = 0
        new_insert_fwd, new_insert_rev = "", ""

    try:
        new_deletion_fwd = Counter(deletions_fwd).most_common(1)[0][0]
        new_deletion_rev = Counter(deletions_rev).most_common(1)[0][0]
        new_deletion_fwd_count = Counter(deletions_fwd).most_common(1)[0][1]
        new_deletion_rev_count = Counter(deletions_rev).most_common(1)[0][1]
        new_deletion_count = new_deletion_fwd_count + new_deletion_rev_count

    except IndexError:
        new_deletion_count = 0
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
        variant.type = "ins"
        variant.variant_nucleotide = new_insert_fwd
        variant.length = len(new_insert_fwd)
        variant.support_count = new_insert_count
        variant.support_ratio = new_insert_count / total_q
        variant.support_count_fwd = new_insert_fwd_count
        variant.support_ratio_fwd = new_insert_fwd_count / fwd_count
        variant.support_count_rev = new_insert_rev_count
        variant.support_ratio_rev = new_insert_rev_count / rev_count
        variant.base_qualities = quality_map
        variant.alt_base_qualities = inserts

    elif deletion_indicator:
        variant.found = True
        variant.type = "del"
        variant.variant_nucleotide = new_deletion_fwd
        variant.length = len(new_deletion_fwd)
        variant.support_count = new_deletion_count
        variant.support_ratio = new_deletion_count / total_q
        variant.support_count_fwd = new_deletion_fwd_count
        variant.support_ratio_fwd = new_deletion_fwd_count / fwd_count
        variant.support_count_rev = new_deletion_rev_count
        variant.support_ratio_rev = new_deletion_rev_count / rev_count
        variant.base_qualities = quality_map
        variant.alt_base_qualities = deletions

    return variant


def extract_indel_evidence(
    pileupcolumn: pysam.PileupColumn,
    assembly: str,
    ref_nt: str,
    reference_mapper: Dict[int, str],
    pos: int,
    amplicon_end: int,
    high_base_quality_cutoff: int = 80,
    edge_of_amplicon: bool = False,
) -> Tuple[str, Tuple[str, str], VCF_entry]:
    """
    Find Indel in a given pileup position.

    Args:
        pileupcolumn: Pysam pileup column at a given position.
        assembly: Reference/contig name.
        ref_nt: Reference nucleotide at a given position.
        reference: Reference genome or sequence, FASTA.
        pos: Position in the alignment pileup to check for variants.
        high_base_quality_cutoff: Cutoff to calculate ratios on high base
            quality nucleotides only (Default = 80).
        edge_of_amplicon: Flag to determine if we are in either the start
            or end of the amplicon.

    Returns:
        A tuple with reference name, reference and alternative alleles,
        and a VCF_entry object with variant information.
    """

    vcf_entry = VCF_entry()
    alleles = (".", ".")
    info = {
        "TYPE": "snp",
        "DP": vcf_entry.DP,
        "QA": vcf_entry.ABQ * vcf_entry.TOTC,
        "AO": vcf_entry.TOTC,
        "SAF": vcf_entry.FWDC,
        "SAR": vcf_entry.REVC,
    }

    if not edge_of_amplicon:
        indel = check_indel(pileupcolumn)
        if indel.found:
            if indel.type == "del":
                indel_stop = pos + indel.length + 1
                if indel_stop <= amplicon_end:
                    # Adjust reference allele to include deleted seq
                    deletion_seq = [
                        reference_mapper[x].upper() for x in range(pos, indel_stop)
                    ]
                    deletion_seq = "".join(deletion_seq).upper()
                    ref_seq = deletion_seq[0]
                    alleles = (deletion_seq, ref_seq)
                else:
                    return (assembly, alleles, vcf_entry, info)
            else:
                alleles = (ref_nt, ref_nt + indel.variant_nucleotide)

            # Update variant statistics with the found indel
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
            vcf_entry.DP = depth
            vcf_entry.DPQ = depth_q
            vcf_entry.FREQ = total_ratio
            vcf_entry.VAF = total_ratio
            vcf_entry.FWDC = fwd_count
            vcf_entry.FWDR = fwd_ratio
            vcf_entry.REVC = rev_count
            vcf_entry.REVR = rev_ratio
            vcf_entry.TOTC = total_count
            vcf_entry.TOTR = total_ratio
            vcf_entry.SAME = 1
            vcf_entry.OBSR = total_ratio
            vcf_entry.ABQ = abq
            vcf_entry.OBQ = obq
            vcf_entry.HCR = hcr

            info = {
                "TYPE": indel.type,
                "DP": vcf_entry.DP,
                "QA": vcf_entry.ABQ * vcf_entry.TOTC,
                "AO": vcf_entry.TOTC,
                "SAF": vcf_entry.FWDC,
                "SAR": vcf_entry.REVC,
            }

    return (assembly, alleles, vcf_entry, info)


def main(
    bam: Path,
    bed: Path,
    fasta: Path,
    output_path: Path,
    pileup_depth: int = 1_000_000,
    min_edge_distance: int = 4,
):
    """
    Run Indel detection over given positions.

    Args:
        bam: Read alignments, BAM.
        bed: Genomic locations, BED.
        fasta: Reference genome or sequence, FASTA.
        output_path: Output path for Indels, VCF.
        pileup_depth: Maximum pileup depth.
        min_edge_distance: Minimum distance to the edge of the amplicon.
            If position distance is lower, do not call variants on it.
    """

    from vcf_tools import initialize_output_vcf

    # logging.debug("started main")
    # Open input files and create empty output VCF
    bam_af = pysam.AlignmentFile(bam, "rb")
    reference = pysam.FastaFile(fasta)
    vcf = initialize_output_vcf(output_path, bam_af.references)

    for bed_file_line in create_bed_lines(bed):
        contig = bed_file_line[0]
        contig_region_start = int(bed_file_line[1])
        contig_region_stop = int(bed_file_line[2])
        contig_reference_mapper = create_ref_mapper(
            fasta, contig, contig_region_start, contig_region_stop
        )

        for pos in range(contig_region_start, contig_region_stop):
            close_to_edge = (
                # close to start
                pos - min_edge_distance <= int(contig_region_start)
                or
                # close to end
                pos + min_edge_distance >= int(contig_region_stop)
            )

            positional_pileup = bam_af.pileup(
                contig,
                pos,
                pos + 1,
                truncate=True,
                max_depth=pileup_depth,
                min_base_quality=10,
                stepper="all",
            )

            ref_seq = reference.fetch(reference=contig, start=pos, end=pos + 1).upper()

            for pileupcolumn in positional_pileup:
                result = extract_indel_evidence(
                    pileupcolumn=pileupcolumn,
                    assembly=contig,
                    ref_nt=ref_seq,
                    reference_mapper=contig_reference_mapper,
                    pos=pos,
                    edge_of_amplicon=close_to_edge,
                    amplicon_end=contig_region_stop,
                )

                pos_alleles = result[1]
                ref_allele = pos_alleles[0]
                variant_format = result[2]
                variant_info = result[3]

                if ref_allele != ".":
                    # Reference allele is not '.', then a variant was found
                    # The 'start' value is 0-based, 'stop' is 1-based
                    r = vcf.new_record(
                        contig=contig, start=pos, alleles=pos_alleles, filter="PASS"
                    )

                    # Write found variant as a new entry to VCF output
                    for fld in fields(variant_format):
                        fld_value = getattr(variant_format, fld.name)
                        if isinstance(fld_value, (float, np.float64, np.float32)):
                            fld_entry = str(f"{getattr(variant_format, fld.name):.6f}")
                        elif isinstance(fld_value, int):
                            fld_entry = str(f"{getattr(variant_format, fld.name)}")
                        else:
                            fld_entry = str(fld_value)

                        r.samples["SAMPLE1"][fld.name] = fld_entry

                    # Write info tag values to VCF output
                    # These values will be used to filter the VCF
                    for tag_id, value in variant_info.items():
                        r.info[tag_id] = value

                    vcf.write(r)

                else:
                    # Reference allele is '.', then no variant was found
                    # Don't write anthing to VCF output file
                    continue

    time.sleep(0.5)
    vcf.close()


if __name__ == "__main__":
    import argparse

    dev = False
    if not dev:
        parser = argparse.ArgumentParser(
            description=("Detect indels in BAM alignment file.")
        )

        parser.add_argument("fasta", type=Path)
        parser.add_argument("bed", type=Path)
        parser.add_argument("bam", type=Path)
        parser.add_argument("vcf_out", type=Path)
        args = parser.parse_args()

        main(args.bam, args.bed, args.fasta, args.vcf_out)

    if dev:
        fasta = Path(
            "/data/projects/ROD_tmp/2f/2da9d6bfb181300787503ab54f79de/GRCh38_renamed_BBCS.fasta"
        )
        bed = Path(
            "/data/projects/ROD_tmp/2f/2da9d6bfb181300787503ab54f79de/FAW08675_roi.bed"
        )
        bam = Path(
            "/data/projects/ROD_tmp/2f/2da9d6bfb181300787503ab54f79de/FAW08675.YM_gt_3.bam"
        )
        vcf_out = Path("./test.vcf")

        main(bam, bed, fasta, vcf_out)
