#!/usr/bin/env python
import logging
import time
from pathlib import Path

import pysam
from detect_indel import extract_indel_evidence
from detect_snp import extract_snp_evidence
from tqdm import tqdm
from vcf_tools import create_bed_positions, initialize_output_vcf, write_vcf_entry


def write_variants(
    reference: pysam.FastaFile,
    contig: str,
    pos: int,
    bam_af: pysam.AlignmentFile,
    snp_vcf: Path,
    indel_vcf: Path,
    pileup_depth: int,
    minimum_base_quality: int = 10,
    high_base_quality_cutoff: int = 80,
    end_of_amplicon: bool = False,
):
    """
    Identify SNPs and Indels in a given pileup position.
    Results are written to VCF output files.

    Args:
        reference: Reference genome or sequence, FASTA.
        contig: Reference sequence name.
        pos: Position in the alignment pileup to check for variants.
        bam_af: Read alignments, BAM.
        snp_vcf: Output path for SNPs, VCF.
        indel_vcf: Output path for Indels, VCF.
        pileup_depth: Maximum pileup depth, integer (Default=1_000_000).
        minimum_base_quality: Minimum base quality to consider for
            variant calling.
        high_base_quality_cutoff: Cutoff to calculate ratios on high base
            quality nucleotides only (Default = 80).
        end_of_amplicon: (Not used) Flag to determine if we are reaching
            the end of the amplicon, based on positions in BED file.
    """

    positional_pileup = bam_af.pileup(
        contig,
        pos,
        pos + 1,
        truncate=True,
        max_depth=pileup_depth,
        min_base_quality=minimum_base_quality,
        stepper="all",
    )

    # list of jobs TODO: contig:position
    # based on contig:position, do snp + indel
    # multiprocessing.pool, run all positions in parallel
    for pileupcolumn in positional_pileup:
        # positional_pileup only has 1 element
        snv_evidence = extract_snp_evidence(
            pileupcolumn, contig, high_base_quality_cutoff
        )
        # TODO: filter here, write immediately to vcf
        # later merge with indels
        write_vcf_entry(snp_vcf, contig, pos, snv_evidence)

        indel_evidence = extract_indel_evidence(
            pileupcolumn,
            contig,
            reference,
            pos,
            high_base_quality_cutoff,
            end_of_amplicon,
        )
        if indel_evidence[1][0] != ".":
            # if reference allele is not '.', then an indel was found
            write_vcf_entry(indel_vcf, contig, pos, indel_evidence)


def main(
    fasta: Path,
    bed: Path,
    bam: Path,
    snp_output_path: Path,
    indel_output_path: Path,
    pileup_depth: int = 1_000_000,
):
    """
    Run variant detection over given positions.

    Args:
        fasta: Reference genome or sequence, FASTA.
        bed: Genomic locations, BED.
        bam: Read alignments, BAM.
        snp_output_path: Output path for SNPs, VCF.
        indel_output_path: Output path for Indels, VCF.
        pileup_depth: Maximum pileup depth, integer (DEFAULT=1_000_000).
    """

    logging.debug("started main")
    bam_af = pysam.AlignmentFile(bam, "rb")
    reference = pysam.FastaFile(fasta)

    snp_vcf = initialize_output_vcf(snp_output_path, bam_af.references)
    indel_vcf = initialize_output_vcf(indel_output_path, bam_af.references)

    for contig, pos, amplicon_ending in tqdm(
        create_bed_positions(bed, end_warning_length=4)
    ):
        # TODO: split detection and writing, so that we can multiprocess the detection
        write_variants(
            reference,
            contig,
            pos,
            bam_af,
            snp_vcf,
            indel_vcf,
            pileup_depth,
            end_of_amplicon=amplicon_ending,
        )

    time.sleep(0.5)
    snp_vcf.close()
    indel_vcf.close()


if __name__ == "__main__":
    import argparse

    dev = False
    if not dev:
        parser = argparse.ArgumentParser(description="")

        parser.add_argument("fasta", type=Path)
        parser.add_argument("bed", type=Path)
        parser.add_argument("bam", type=Path)
        parser.add_argument("snp_vcf_out", type=Path)
        parser.add_argument("indel_vcf_out", type=Path)
        args = parser.parse_args()
        logging.info(args)

        main(args.fasta, args.bed, args.bam, args.snp_vcf_out, args.indel_vcf_out)

    if dev:
        # EGFR
        fasta = Path(
            "/data/references/Homo_sapiens/GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna"
        )
        bed = Path("/data/projects/ROD_1125_variant_improvements/EGFR.bed")
        bam = Path(
            "/data/projects/ROD_1125_variant_improvements/ONT_20221121_EGFR/consensus_aligned/284.taged.bam"
        )
        snp_vcf_out = Path(
            "/data/projects/ROD_1125_variant_improvements/testsnp_EGFR.vcf"
        )
        indel_vcf_out = Path(
            "/data/projects/ROD_1125_variant_improvements/testindel_EGFR.vcf"
        )
        main(fasta, bed, bam, snp_vcf_out, indel_vcf_out)
