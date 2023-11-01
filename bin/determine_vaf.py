#!/usr/bin/env python
import logging
import time
from concurrent.futures import ALL_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path
from threading import Lock

import pysam
from variant_calling.detect_indel import extract_indel_evidence
from variant_calling.detect_snp import extract_snp_evidence
from variant_calling.vcf_tools import (
    create_bed_lines,
    initialize_output_vcf,
    write_vcf_entry,
)


def process_pileup_column(
    contig: str,
    pos: int,
    bam: str,
    reference_fasta: str,
    amplicon_ending: bool = False,
    pileup_depth: int = 1_000_000,
    minimum_base_quality: int = 10,
):
    """
    Function to process a single location in a genome alignment. Used to multiprocess the variant calling

    Args:
        contig: Reference genome contig name.
        pos: Genomic locations within the contig.
        bam: path to the bam file, BAM.
        fasta: path to the reference genome, Fasta.
        amplicon_ending: Bolean to indicate if the amplicon is about to end, used in indel calling.
        pileup_depth: Maximum pileup depth, integer (DEFAULT=1_000_000).
        minimum_base_quality: minimum base quality at a given position in a read to be considered, integer.
    """
    logging.debug(f"process column {contig} | {pos}")
    bam_af = pysam.AlignmentFile(bam, "r")
    # Only have one process access the reference genome at a time, multi access is not supported on all filesystems
    with Lock() as lock:
        # Double try to read the reference genome, with a sleep time for the retry
        try:
            reference = pysam.FastaFile(reference_fasta)
        except OSError:
            try:
                time.sleep(0.2)
                reference = pysam.FastaFile(reference_fasta)
            except:
                raise OSError(f"Reference genome {reference_fasta} could not be read.")

        print(reference.references)
        try:
            ref_nuc = str(
                reference.fetch(reference=contig, start=pos, end=pos + 1)
            ).upper()
        except KeyError:
            try:
                ref_nuc = str(
                    reference.fetch(reference=contig, start=pos, end=pos + 1)
                ).upper()
            except KeyError:
                ref_nuc = "."
                logging.warning(f"Unable to fetch reference genome at {contig} - {pos}")

        del reference

    iteration = 0
    # create enteties for when no forloop event happens
    snv_evidence = None
    indel_evidence = None

    for pu_column in bam_af.pileup(
        contig,
        pos,
        pos + 1,
        truncate=True,
        max_depth=pileup_depth,
        min_base_quality=minimum_base_quality,
        stepper="all",
    ):
        snv_evidence = extract_snp_evidence(
            pu_column, contig, ref_nuc, minimum_base_quality
        )
        indel_evidence = extract_indel_evidence(
            pu_column,
            contig,
            ref_nuc,
            reference_fasta,
            pos,
            minimum_base_quality,
            amplicon_ending,
        )
        if iteration > 0:
            logging.critical("Single thread worker recieved multiple locations.")

        if indel_evidence[1][0] == "." and indel_evidence[1][1] == ".":
            indel_evidence = None

        logging.debug(f"done process column {contig} | {pos}")
        return (contig, pos, snv_evidence, indel_evidence)


def main(
    fasta: Path,
    bed_path: Path,
    bam_path: Path,
    snp_output_path: Path,
    indel_output_path: Path,
    pileup_depth: int = 1_000_000,
    minimum_base_quality: int = 10,
    end_of_amplicon_warn_limit: int = 4,
    threads: int = 16,
):
    """
    Run variant detection over given positions.

    Args:
        fasta: Reference genome or sequence, FASTA.
        bed_file_path: Genomic locations, BED.
        bam_path: Read alignments, BAM.
        snp_output_path: Output path for SNPs, VCF.
        indel_output_path: Output path for Indels, VCF.
        pileup_depth: Maximum pileup depth, integer (DEFAULT=1_000_000).
    """

    logging.debug("started main")
    pileup_parameters = []
    promisses = []

    with ProcessPoolExecutor(max_workers=threads) as executor:
        logging.debug("creating location promisses in ProcessPoolExecutor")
        for bed_file_line in create_bed_lines(bed_path):
            print(bed_file_line)
            contig = bed_file_line[0]
            contig_region_start = int(bed_file_line[1])
            contig_region_stop = int(bed_file_line[2])

            for pos in range(contig_region_start, contig_region_stop):
                amplicon_ending = (
                    True
                    if contig_region_stop - pos > end_of_amplicon_warn_limit
                    else False
                )
                logging.debug("appending to ProcessPoolExecutor")
                promisses.append(
                    executor.submit(
                        process_pileup_column,
                        contig,
                        pos,
                        bam_path,
                        fasta,
                        amplicon_ending,
                        pileup_depth,
                        minimum_base_quality,
                    )
                )
        print("done with submitting")
        logging.debug("Asking for results from ProcessPoolExecutor")
        done, not_done = wait(promisses, return_when=ALL_COMPLETED)
        print("done with collecting")
        results = [x.result() for x in promisses]
        print("done with obtaining results")

    bam_af = pysam.AlignmentFile(bam_path, "r")

    snp_vcf = initialize_output_vcf(snp_output_path, bam_af.references)
    indel_vcf = initialize_output_vcf(indel_output_path, bam_af.references)

    print("Started writing vcf files")
    for result in results:
        if not result:
            continue
        contig = result[0]
        position = result[1]
        snp_results = result[2]
        indel_results = result[3]

        if snp_results:
            write_vcf_entry(snp_vcf, contig, position, snp_results)

        if indel_results:
            write_vcf_entry(indel_vcf, contig, position, indel_results)

    print("finishing writing and closing files.")
    # Wait for all the threads to finish writing
    delay_for_pysam_variantfile = 0.5
    time.sleep(delay_for_pysam_variantfile)
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
        parser.add_argument("-t", "--threads", default=16, type=int)
        args = parser.parse_args()
        logging.info(args)

        main(args.fasta, args.bed, args.bam, args.snp_vcf_out, args.indel_vcf_out)

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
        snp_vcf_out = Path("./test_snp.vcf")
        indel_vcf_out = Path("./test_indel.vcf")

        main(fasta, bed, bam, snp_vcf_out, indel_vcf_out)
