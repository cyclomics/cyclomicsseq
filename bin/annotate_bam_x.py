#!/usr/bin/env python

from typing import Dict, Tuple
import argparse
import pathlib

import pysam
import logging

header_to_tag = {
    "channel": "XC",
    "mux": "XM",
    "channel_mux": "Xc",
    "start_time": "XT",
    "duration": "XD",
    "adapter_duration": "XA",
    "sequence_length_template": "XL",
    "mean_qscore_template": "XQ",
    "median_template": "XT",
    "mad_template": "XM",
    "end_reason": "XE",
}


def load_seqsum(
    seqsum_path: str, readname_index=3
) -> Tuple[Dict[str, str], Dict[str, int]]:
    """
    Load in the sequnecing summary file.
    """
    seqsum_indexes = {}
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


def tag_aln(
    aln,
    query_name,
    seqsum,
    seqsum_indexes,
    items=[
        "channel",
        "mux",
        "channel_mux",
        "start_time",
        "duration",
        "adapter_duration",
        "sequence_length_template",
        "mean_qscore_template",
        "median_template",
        "mad_template",
        "end_reason",
    ],
):
    """
    for a set of information in `items`, look up the index and the tag and add it to the tags list of aln
    """
    readsum = seqsum[query_name]
    new_tags = []
    for i in items:
        tag = header_to_tag[i]
        value = int(readsum[seqsum_indexes[i]])

        new_tags.append((tag, value ))

    aln.tags = aln.tags + new_tags
    return aln


def main(seqsum_path, in_bam_path, out_bam_path, split_queryname=True):
    # Initialize file handles
    in_bam_file = pysam.AlignmentFile(in_bam_path, "rb")
    seqsum, seqsum_indexes = load_seqsum(seqsum_path)
    out_bam_file = pysam.AlignmentFile(out_bam_path, "wb", header=in_bam_file.header)

    for aln in in_bam_file.fetch():
        query_name = aln.query_name
        if split_queryname:
            query_name = query_name.split("_")[0]
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

    parser.add_argument("file_seqsum", type=pathlib.Path)
    parser.add_argument("file_bam", type=pathlib.Path)
    parser.add_argument("file_out", type=pathlib.Path)
    args = parser.parse_args()
    logging.info(args)

    intermediate_bam = args.file_bam.with_suffix(".unsorted.bam")
    logging.info(intermediate_bam)

    pore_time_aln = main(args.file_seqsum, args.file_bam, intermediate_bam)

    # pysam requires pure strings
    pysam.sort("-o", str(args.file_out), str(intermediate_bam))
    pysam.index(str(args.file_out))
