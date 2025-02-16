#!/usr/bin/env python

import argparse
from pathlib import Path
from typing import Dict, Union

import pysam

FilePath = Union[str, Path]


def get_amplicon_lengths(in_file: str) -> Dict[str, int]:
    """Create a dictionary from a TSV file with amplicon names and expected lengths.

    Args:
        in_file (FilePath): Path to input BED file

    Returns:
        Dict[str, int]: Dictionary where keys are amplicon names,
            values are expected amplicon lengths
    """
    amplicon_lengths = {}
    with open(in_file, "r") as file:
        for line in file:
            fields = line.strip().split("\t")
            amplicon_length = int(fields[2]) - int(fields[1])
            amplicon_name = fields[3]

            amplicon_lengths[amplicon_name] = amplicon_length

    return amplicon_lengths


def filter_alignment_rate(
    in_file: str, out_file: str, amplicons_file: str, min_rate: float
):
    """Filter alignments in input BAM based on contig-specific minimum alignment rates.

    Args:
        in_file (FilePath): Path to input BAM file to be filtered by alignment rate.
        out_file (FilePath): Path to output BAM file after filtering by alignment rate.
        amplicons_file (FilePath): Path to input TSV file containing contig-specific minimum
            alignment rates (column 1: amplicon name; column 2: expected amplicon length).
        min_rate (float): Minimum alignment rate for a read to be kept. Alignment rate is calculated
            as (alignment end in reference - alignment start in reference) / expected amplicon length
    """
    input_bam = pysam.AlignmentFile(in_file, "rb")
    output_bam = pysam.AlignmentFile(out_file, "wb", template=input_bam)
    amplicon_lengths = get_amplicon_lengths(amplicons_file)

    for read in input_bam:
        if read.reference_name in amplicon_lengths.keys():
            # If reference contig name is found, then filter based on minimum alignment rate
            align_rate = (read.reference_end - read.reference_start) / amplicon_lengths[
                read.reference_name
            ]
            if align_rate > min_rate:
                output_bam.write(read)
        else:
            # If reference contig name is not found, write read to output anyway
            output_bam.write(read)

    input_bam.close()
    output_bam.close()


if __name__ == "__main__":
    dev = False
    if not dev:
        parser = argparse.ArgumentParser(
            description="Filter alignments based on contig-specific minimum alignment rates."
        )
        parser.add_argument("--input_bam", type=str)
        parser.add_argument("--output_bam", type=str)
        parser.add_argument("--amplicons_bed", type=str)
        parser.add_argument("--min_align_rate", type=float)
        args = parser.parse_args()

        filter_alignment_rate(
            args.input_bam, args.output_bam, args.amplicons_bed, args.min_align_rate
        )

    else:
        input_bam = "PAW56644.YM_gt_3.bam"
        output_bam = "test.bam"
        amplicons_bed = "bin/variant_calling/egfr_panel.bed"
        min_align_rate = 0.8

        filter_alignment_rate(input_bam, output_bam, amplicons_bed, min_align_rate)
