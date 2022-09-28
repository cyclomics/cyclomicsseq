#!/usr/bin/env python

from dataclasses import dataclass
import logging
from collections import Counter, namedtuple
from pathlib import Path
import time
import re

import pysam
from tqdm import tqdm
import numpy as np


def is_intable(value):
    """Test if the input can be converted to an int"""

    try:
        float(value)
        return True
    except ValueError:
        return False


def create_bed_positions(bed_file, end_warning_length=4):
    """Returns positions on which to pileup reads and find variants"""

    with open(bed_file) as f:
        for line in f:
            # bed should be tab delimited
            L = line.strip().split("\t")
            if len(L) == 1:
                L = L[0].split(" ")
                L = [x for x in L if x]

            # if not is_intable(L[2]):
            # logging.critical("error in bed file")

            for pos in range(int(L[1]), int(L[2])):
                close_to_end = pos + end_warning_length >= int(L[2])
                yield L[0], pos, close_to_end


def initialize_output_vcf(vcf_path, contigs):
    """Returns pre-formatted VCF file on which to output variants"""

    # Create a VCF header
    vcfh = pysam.VariantHeader()
    # Add a sample named "ahstram" to our VCF header
    vcfh.add_sample("Sample1")
    # Add a "FILTER" value, other than "PASS"
    vcfh.add_meta(
        "FILTER",
        items=[("ID", "RF"), ("Description", "Variant failed filter due to low RF")],
    )
    # Add a contig (chromosome 1) to our VCF header
    for contig in contigs:
        vcfh.add_meta("contig", items=[("ID", contig)])

    # Add GT to FORMAT in our VCF header
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "DP"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Depth at position"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "DPQ"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Depth at position above Q threshold"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "FREQ"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Frequency, raw frequency of non majority allele"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "VAF"),
            ("Number", 1),
            ("Type", "String"),
            (
                "Description",
                "Variant Allele Frequency, filtered frequency on base quaility of non majority allele",
            ),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "FWDC"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Forward count after Q filtering"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "FWDR"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Forward ratio after Q filtering"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "REVC"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Reverse count after Q filtering"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "REVR"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Reverse ratio after Q filtering"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "TOTC"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Total count after Q filtering"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "TOTR"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Total ratio after Q filtering"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "SAME"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Same base found forward and reverse"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "OBSR"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "observational ratio of alternative allele"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "ABQ"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "alternative base quality (mean)"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "OBQ"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "overall base quality (mean)"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "HCR"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "high quality ratio"),
        ],
    )

    vcf = pysam.VariantFile(vcf_path, "w", header=vcfh)
    return vcf


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


def check_indel(pu_column, variant_th=0.003, variant_count_th=10):
    """Return most common indel above thresholds for a pileup column

    Checks pileupcolumn for indel evidence from
    pysam.PileupColumn.get_query_sequences,
    then chooses most common insert or deletion to call as a variant.
    """

    # TODO: Check homopolymer status
    # Initialize general statistics of this pileup position
    total = pu_column.nsegments
    total_q = pu_column.get_num_aligned()
    variant = Indel(None)
    variant.depth = total
    variant.depth_q = total_q
    variant.assemby = pu_column.reference_name
    variant.start_position = pu_column.pos

    # Get positional sequences, from which to search for indels
    pu = pu_column.get_query_sequences(
        mark_matches=True, mark_ends=True, add_indels=True
    )

    # The reference nucleotide should be the most common in this position
    variant.reference_nucleotide = Counter(pu).most_common(1)[0][0].upper()

    # Exclude matches, SNPs or irrelevant characters
    # The remainder should be indels
    excluded_chars = [".", ",", "*", "^", "A", "C", "T", "G", "]"]
    pu_indel = [x for x in pu if x.upper() not in excluded_chars]

    inserts = [x for x in pu_indel if x[1] == "+"]
    # remove the directional indicator and make uppercase
    # eg convert A+3TTT' into TTT and A+1a to A
    inserts_fwd = [x[3:].upper() for x in inserts if x[-1].isupper()]
    inserts_rev = [x[3:].upper() for x in inserts if x[-1].islower()]

    deletions = [x for x in pu_indel if x[1] == "-"]
    # remove the directional indicator and make uppercase
    # eg G-15NNNNNNNNNNNNNNN g-15nnnnnnnnnnnnnnn turn into
    # 15NN... and 15NN.. in their respective variable
    deletions_fwd = [x[2:].upper() for x in deletions if x[-1] == "N"]
    deletions_rev = [x[2:].upper() for x in deletions if x[-1] == "n"]

    # If no insert or deletion evidence was found in this position,
    # the 'inserts' list will be empty and the Counter will return
    # an IndexError, which can be understood as 'no evidence'.
    try:
        new_insert_fwd = Counter(inserts_fwd).most_common(1)[0][0]
        new_insert_rev = Counter(inserts_rev).most_common(1)[0][0]
        new_insert_fwd_count = Counter(inserts_fwd).most_common(1)[0][1]
        new_insert_rev_count = Counter(inserts_rev).most_common(1)[0][1]
        new_insert_count = new_insert_fwd_count + new_insert_rev_count

    except IndexError:
        new_insert_count = 0
        new_insert_fwd, new_insert_rev = "", ""

    try:
        new_deletion_fwd = Counter(deletions_fwd).most_common(1)[0][0][1:]
        new_deletion_rev = Counter(deletions_rev).most_common(1)[0][0][1:]
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
            and new_insert_fwd_count / total_q > variant_th
            and new_insert_fwd_count > variant_count_th
            and new_insert_rev_count / total_q > variant_th
            and new_insert_rev_count > variant_count_th
        ):
            # Found an insertion
            insert_indicator = True

        elif (
            new_deletion_fwd == new_deletion_rev
            and new_deletion_count > new_insert_count
            and new_deletion_fwd_count / total_q > variant_th
            and new_deletion_fwd_count > variant_count_th
            and new_deletion_rev_count / total_q > variant_th
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
        variant.variant_nucleotide = variant.reference_nucleotide + new_insert_fwd
        variant.length = len(new_insert_fwd)
        variant.support_count = new_insert_count
        variant.support_ratio = new_insert_count / total_q
        variant.support_count_fwd = new_insert_fwd_count
        variant.support_ratio_fwd = new_insert_fwd_count / total_q
        variant.support_count_rev = new_insert_rev_count
        variant.support_ratio_rev = new_insert_rev_count / total_q

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

    return variant


def extract_nucleotide_count(
    bam: pysam.AlignmentFile,
    assembly: str,
    reference: pysam.FastaFile,
    pos: int,
    pileup_depth: int,
    minimum_base_quality=10,
    high_base_quality_cutoff=80,
    end_of_amplicon=False,
):
    """Returns found variants in a given pileup position"""

    # Initialize VCF tags
    tags = [
        "DP",
        "DPQ",
        "FREQ",
        "VAF",
        "FWDC",
        "FWDR",
        "REVC",
        "REVR",
        "TOTC",
        "TOTR",
        "SAME",
        "OBSR",
        "ABQ",
        "OBQ",
        "HCR",
    ]

    # Initialize the position with no variant
    Variant = namedtuple("Variant", tags)
    alleles = (".", ".")
    indel_type = ""
    depth = 0
    depth_q = 0
    fwd_count = 0
    fwd_ratio = 0
    rev_count = 0
    rev_ratio = 0
    total_count = 0
    total_ratio = 0

    # Create a pileup for a given position in the search space
    positional_pileup = bam.pileup(
        assembly,
        pos,
        pos + 1,
        truncate=True,
        max_depth=pileup_depth,
        min_base_quality=minimum_base_quality,
        stepper="all",
    )

    # Each positional pileup will have 1 pileup column
    for pileupcolumn in positional_pileup:
        # Check for Indels
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
                    alleles = (indel.reference_nucleotide, indel.variant_nucleotide)

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

    # Return statistics for this position,
    # whether a variant was found or not
    var = Variant(
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
        OBSR="N/A",
        ABQ="N/A",
        OBQ="N/A",
        HCR="N/A",
    )

    return (assembly, alleles, var, indel_type)


def main(bam: Path, bed: Path, fasta: Path, output_path: Path, pileup_depth=1_000_000):
    """Run indel detection

    This will output a separate VCF file only with detected indels.
    Requires as input:
    - BAM file with read alignments,
    - BED file with genomic positions over which to detect indels,
    - FASTA file with the reference genome used for read alignment,
    - Output path to which a VCF result file will be written,
    - Maximum pileup depth (Default=1_000_000)
    """

    # logging.debug("started main")
    # Open input files and create empty output VCF
    bam_af = pysam.AlignmentFile(bam, "rb")
    reference = pysam.FastaFile(fasta)
    vcf = initialize_output_vcf(output_path, bam_af.references)

    # Iterate over positions in search space indicated in BED file
    for assembly, pos, amplicon_ending in tqdm(create_bed_positions(bed)):
        # Check statistics for this position,
        # potentially finding a new variant

        # TODO: results is a list
        result = extract_nucleotide_count(
            bam=bam_af,
            assembly=assembly,
            reference=reference,
            pos=pos,
            pileup_depth=pileup_depth,
            end_of_amplicon=amplicon_ending,
        )

        # TODO: for result in results, write to VCF
        if result[1][0] != ".":
            # Reference allele is not '.', then a variant was found
            # The 'start' value is 0-based, 'stop' is 1-based
            # if result[-1] == 'deletion':
            #    pos = pos + len(result[1][0]) - 2

            r = vcf.new_record(
                contig=assembly, start=pos + 1, alleles=result[1], filter="PASS"
            )

            # Write found variant as a new entry to VCF output
            for fld in result[2]._fields:
                fld_value = getattr(result[2], fld)
                if type(fld_value) in [float, np.float64, np.float32]:
                    fld_entry = str(f"{getattr(result[2], fld):.6f}")
                elif type(fld_value) == int:
                    fld_entry = str(f"{getattr(result[2], fld)}")
                else:
                    fld_entry = str(fld_value)

                r.samples["Sample1"][fld] = fld_entry

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
        fasta = Path("tmp/PNK_01_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna")
        # PNK_01
        bed = Path("tmp/PNK_Rob_custom_GRCh38.bed")
        bam = Path("tmp/PNK_01_GRCh38.p14/FAS12641.taged.bam")
        vcf_out = Path("tmp/PNK_01_GRCh38.p14/PNK_01_GRCh38.p14.indel.vcf")

        # Cyclomics TP53
        # bed = Path("tmp/TP53_custom.bed")
        # bam = Path("tmp/0.7.0/FAU48563.taged.bam")
        # vcf_out = Path("tmp/0.7.0/TP53_000025_indel.vcf")

        main(bam, bed, fasta, vcf_out)
