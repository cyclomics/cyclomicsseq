#!/usr/bin/env python

import argparse
import logging
from collections import defaultdict
from email.policy import default
from pathlib import Path
from typing import Dict, Tuple

import pysam

HEADER_TO_TAG = defaultdict(lambda: -1)

HEADER_TO_TAG["adapter_duration"] = "XA"
HEADER_TO_TAG["channel"] = "XC"
HEADER_TO_TAG["duration"] = "XD"
HEADER_TO_TAG["end_reason"] = "XE"
HEADER_TO_TAG["sequence_length_template"] = "XL"
HEADER_TO_TAG["mad_template"] = "XM"
HEADER_TO_TAG["median_template"] = "XP"
HEADER_TO_TAG["start_time"] = "XT"
HEADER_TO_TAG["mean_qscore_template"] = "XQ"
HEADER_TO_TAG["mux"] = "XX"


USED_TAGS = list(HEADER_TO_TAG.keys())


def load_seqsum(
    seqsum_path: str, readname_index=3
) -> Tuple[Dict[str, str], Dict[str, int]]:
    """
    Load in the sequnecing summary file.
    """
    seqsum_indexes = defaultdict(lambda: -1)

    logging.info("Loading sequencing summary")
    with open(seqsum_path, "r") as seqsum:
        header = seqsum.readline()
        for i, name in enumerate(header.split()):
            seqsum_indexes[name] = i
            logging.info("\t", i, name)

        lines = [line.split() for line in seqsum.readlines()]
        readname_info = {x[readname_index]: x for x in lines}
    logging.info("- Done")
    return readname_info, seqsum_indexes


def tag_aln(aln, query_name, seqsum, seqsum_indexes, items=USED_TAGS):
    """
    for a set of information in `items`, look up the index and the tag and add it to the tags list of aln
    If the tag is not in the summary, or an item in the summary is not in tags, dont do anything.
    """
    readsum = seqsum[query_name]
    new_tags = []
    for i in items:
        tag = HEADER_TO_TAG[i]
        header_index = seqsum_indexes[i]

        # tagless column
        if tag == -1:
            continue

        # something not in the sequencing summary
        if header_index == -1:
            continue

        value = readsum[header_index]

        new_tags.append((tag, value))

    aln.tags = aln.tags + new_tags
    return aln


def process_query_name(name, split_queryname=True, splitter="_"):
    if not split_queryname:
        return name

    return name.split(splitter)[0]


def main(seqsum_path, in_bam_path, out_bam_path, split_queryname=True):
    # Initialize file handles
    in_bam_file = pysam.AlignmentFile(in_bam_path, "rb")
    seqsum, seqsum_indexes = load_seqsum(seqsum_path)
    out_bam_file = pysam.AlignmentFile(out_bam_path, "wb", header=in_bam_file.header)

    for aln in in_bam_file.fetch():
        query_name = process_query_name(aln.query_name, split_queryname, "_")

        if query_name in seqsum:
            aln = tag_aln(aln, query_name, seqsum, seqsum_indexes)
        out_bam_file.write(aln)
        # exit()

    in_bam_file.close()
    out_bam_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process the information in the sequencing summary and add it to the bam."
    )

    parser.add_argument("file_seqsum", type=Path)
    parser.add_argument("file_bam", type=Path)
    parser.add_argument("file_out", type=Path)
    args = parser.parse_args()
    logging.info(args)

    intermediate_bam = args.file_bam.with_suffix(".unsorted.bam")
    logging.info(intermediate_bam)

    pore_time_aln = main(args.file_seqsum, args.file_bam, intermediate_bam)

    # pysam requires pure strings
    pysam.sort("-o", str(args.file_out), str(intermediate_bam))
    pysam.index(str(args.file_out))

    # file_seqsum = Path('/media/dami/cyclomics_003/ont/001_accuracy_testing/SS_220209_cyclomics/002/20220209_1609_X5_FAS04073_a6645565/sequencing_summary_FAS04073_77432eb4.txt')
    # file_bam = Path('/media/dami/cyclomics_003/tmp/annotationtest/Minimap2Align/SamToBam/FAS04073_pass_77432eb4_140.fastq.bam')
    # file_out = Path('testing.bam')
    # intermediate_bam = file_out.with_suffix(".unsorted.bam")

    # pore_time_aln = main(file_seqsum, file_bam, intermediate_bam)
