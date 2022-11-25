#!/usr/bin/env python
import logging
import time
from pathlib import Path

import pysam
from detect_indel import *
from detect_snp import *
from tqdm import tqdm
from vcf_tools import create_bed_positions, initialize_output_vcf, write_vcf_entry


def write_variants(
    snp_vcf,
    indel_vcf,
    bam,
    contig,
    reference,
    pos,
    pileup_depth,
    minimum_base_quality=10,
    add_dels=False,
    high_base_quality_cutoff=80,
    end_of_amplicon=False,
):
    """
    Returns found variants in a given pileup position
    """
    positional_pileup = bam.pileup(
        contig,
        pos,
        pos + 1,
        truncate=True,
        max_depth=pileup_depth,
        min_base_quality=minimum_base_quality,
        stepper="all",
    )

    # list of jobs to do: contig:position
    # based on contig:position, do snp + indel
    # multiprocessing.pool, run all positions in parallel

    # This whole thing becomes a single nxf process

    for pileupcolumn in positional_pileup:
        # positional_pileup only has 1 element
        snv_evidence = extract_snv_evidence(pileupcolumn)
        # filter here, write immediately to vcf
        # later merge with indels
        write_vcf_entry(snp_vcf, contig, pos, snv_evidence)

        indel_evidence = extract_indel_evidence(pileupcolumn, depth=5000)
        if indel_evidence[1][0] != ".":
            # if reference allele is not '.', then an indel was found
            write_vcf_entry(indel_vcf, contig, pos, indel_evidence)


def main(
    bam: Path,
    variants: Path,
    snp_output_path: Path,
    indel_output_path: Path,
    pileup_depth: int = 1_000_000,
):
    """ """

    logging.debug("started main")
    bam_af = pysam.AlignmentFile(bam, "rb")

    snp_vcf = initialize_output_vcf(snp_output_path, bam_af.references)
    indel_vcf = initialize_output_vcf(indel_output_path, bam_af.references)

    for contig, pos, amplicon_ending in tqdm(create_bed_positions(variants)):
        result = write_variants(
            snp_vcf,
            indel_vcf,
            bam_af,
            contig,
            pos,
            pileup_depth,
            end_of_amplicon=amplicon_ending,
        )

    time.sleep(0.5)
    snp_vcf.close()
    indel_vcf.close()


if __name__ == "__main__":
    import argparse

    dev = True
    if not dev:
        parser = argparse.ArgumentParser(
            description="Process the information in the sequencing summary and add it to the bam."
        )

        parser.add_argument("variant_bed", type=Path)
        parser.add_argument("bam", type=Path)
        parser.add_argument("snp_file_out", type=Path)
        parser.add_argument("indel_file_out", type=Path)
        args = parser.parse_args()
        logging.info(args)

        main(args.bam, args.variant_bed, args.snp_file_out, args.indel_file_out)

    if dev:
        # PNK
        # bed = Path("dilution_series_expected_mutations.bed")
        bed = Path("/home/dami/projects/variantcalling/depth/PNK_region.bed")
        bam = Path("/home/dami/Data/dev/ROB_pnk/FAS12639.YM_gt_3.bam")

        main(bam, bed, "tmp.vcf")
