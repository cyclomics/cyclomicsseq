#!/usr/bin/env python
import logging
import time
from concurrent.futures import ALL_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path
from threading import Lock
from typing import Dict

import pysam
from variant_calling.detect_indel import extract_indel_evidence
from variant_calling.detect_snp import extract_snp_evidence
from variant_calling.vcf_tools import (create_bed_lines, initialize_output_vcf,
                                       write_vcf_entry)


def process_pileup_column(
    contig: str,
    pos: int,
    bam: str,
    reference_mapper: Dict[int,str],
    amplicon_edge: bool = False,
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
        amplicon_edge: Bolean to indicate if current position is at end of amplicon, used in indel calling.
        pileup_depth: Maximum pileup depth, integer (DEFAULT=1_000_000).
        minimum_base_quality: minimum base quality at a given position in a read to be considered, integer.
    """
    logging.debug(f"process column {contig} | {pos}")
    
    # obtain the base for the snp
    ref_nuc = reference_mapper[pos]
    
    # create enteties for when no forloop event happens
    iteration = 0
    snv_evidence = None
    indel_evidence = None

    bam_af = pysam.AlignmentFile(bam, "r")

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
        # indel needs the full mapper
        indel_evidence = extract_indel_evidence(
            pu_column,
            contig,
            ref_nuc,
            reference_mapper,
            pos,
            minimum_base_quality,
            amplicon_edge,
        )
        if iteration > 0:
            logging.critical("Single thread worker recieved multiple locations.")

        if indel_evidence[1][0] == "." and indel_evidence[1][1] == ".":
            indel_evidence = None

        logging.debug(f"done process column {contig} | {pos}")
        return (contig, pos, snv_evidence, indel_evidence)


def create_ref_mapper(reference_fasta, contig,start,stop)-> Dict[int,str]:
    """
    Code to create a dict with pos -> base for a given contig.
    """
    try:
        reference = pysam.FastaFile(reference_fasta)
    except OSError:
        try:
            time.sleep(0.2)
            reference = pysam.FastaFile(reference_fasta)
        except:
            raise OSError(f"Reference genome {reference_fasta} could not be read.")

    ref_mapper = {}
    for pos in range(start,stop+1):
        ref_nuc = str(
                reference.fetch(reference=contig, start=pos, end=pos + 1)
            ).upper()
        ref_mapper[pos] = ref_nuc
    
    return ref_mapper


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
    promisses = []

    #  Loop over all lines in the bed file and start a process.
    with ProcessPoolExecutor(max_workers=threads) as executor:
        logging.debug("creating location promisses in ProcessPoolExecutor")
        for bed_file_line in create_bed_lines(bed_path):
            logging.debug(f"processing {bed_file_line}")
            contig = bed_file_line[0]
            contig_region_start = int(bed_file_line[1])
            contig_region_stop = int(bed_file_line[2])
            contig_reference_mapper = create_ref_mapper(fasta, contig, contig_region_start, contig_region_stop)
            for pos in range(contig_region_start, contig_region_stop):
                amplicon_edge = (
                    # amplicon start edge
                    pos - contig_region_start <= end_of_amplicon_warn_limit
                    or
                    # amplicon stop edge
                    contig_region_stop - pos <= end_of_amplicon_warn_limit
                )
                logging.debug("appending to ProcessPoolExecutor")
                promisses.append(
                    executor.submit(
                        process_pileup_column,
                        contig,
                        pos,
                        bam_path,
                        contig_reference_mapper,
                        amplicon_edge,
                        pileup_depth,
                        minimum_base_quality,
                    )
                )
        logging.info("done with submitting")
        logging.debug("Asking for results from ProcessPoolExecutor")
        done, not_done = wait(promisses, return_when=ALL_COMPLETED)
        logging.info("done with collecting")
        results = [x.result() for x in promisses]
        logging.info("done with obtaining results")

    bam_af = pysam.AlignmentFile(bam_path, "r")

    snp_vcf = initialize_output_vcf(snp_output_path, bam_af.references)
    indel_vcf = initialize_output_vcf(indel_output_path, bam_af.references)

    logging.info("Started writing vcf files")
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

    logging.info("finishing writing and closing files.")
    # Wait for all the threads to finish writing
    delay_for_pysam_variantfile = 0.5
    time.sleep(delay_for_pysam_variantfile)
    snp_vcf.close()
    indel_vcf.close()
    logging.info(f"done, results are in :\n SNP: {snp_output_path}\n INDEL: {indel_output_path}")


if __name__ == "__main__":
    import argparse

    dev = True
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
            "/home/dami/Software/cyclomicsseq/dev_files/chm13v2_BBCS.fasta"
        )
        bed = Path(
            "/home/dami/Software/cyclomicsseq/dev_files/fastq_roi.bed"
        )
        bam = Path(
            "/home/dami/Software/cyclomicsseq/dev_files/fastq.YM_gt_3.bam"
        )
        snp_vcf_out = Path("./test_snp.vcf")
        indel_vcf_out = Path("./test_indel.vcf")

        main(fasta, bed, bam, snp_vcf_out, indel_vcf_out)
