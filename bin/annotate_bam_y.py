#!/usr/bin/env python

import sys

import pysam


def load_json(seqsum_path, readname_index=3):
    seqsum_indexes = {}
    print("Loading sequencing summary")
    with open(seqsum_path, "r") as seqsum:
        header = seqsum.readline()
        for i, name in enumerate(header.split()):
            seqsum_indexes[name] = i
            print("\t", i, name)
        lines = [line.split() for line in seqsum.readlines()]
        readnames = [x[readname_index] for x in lines]
    print("- Done")
    return header, dict(zip(readnames, lines)), seqsum_indexes


def add_y_tags(seqsum_path, in_bam_path, out_bam_path, split_queryname=True):
    # Initialize file handles
    in_bam_file = pysam.AlignmentFile(in_bam_path, "rb")
    hss, seqsum, seqsum_indexes = load_json(seqsum_path)
    out_bam_file = pysam.AlignmentFile(
        out_bam_path + ".bam", "wb", header=in_bam_file.header
    )

    for aln in in_bam_file.fetch():
        query_name = aln.query_name
        if split_queryname:
            query_name = query_name.split("_")[0]
        if query_name in seqsum:
            readsum = seqsum[query_name]
            # TODO: make this extenable
            aln.tags = aln.tags + [
                (header_to_tag["channel"], int(readsum[seqsum_indexes["channel"]])),
                (header_to_tag["mux"], int(readsum[seqsum_indexes["mux"]])),
                (
                    header_to_tag["channel_mux"],
                    readsum[seqsum_indexes["channel"]]
                    + "."
                    + readsum[seqsum_indexes["mux"]],
                ),
                (
                    header_to_tag["start_time"],
                    float(readsum[seqsum_indexes["start_time"]]),
                ),
                (header_to_tag["duration"], float(readsum[seqsum_indexes["duration"]])),
                (
                    header_to_tag["adapter_duration"],
                    float(readsum[seqsum_indexes["duration"]])
                    - float(readsum[seqsum_indexes["start_time"]]),
                ),
                (
                    header_to_tag["sequence_length_template"],
                    int(readsum[seqsum_indexes["sequence_length_template"]]),
                ),
                (
                    header_to_tag["mean_qscore_template"],
                    float(readsum[seqsum_indexes["mean_qscore_template"]]),
                ),
                (
                    header_to_tag["median_template"],
                    float(readsum[seqsum_indexes["median_template"]]),
                ),
                (
                    header_to_tag["mad_template"],
                    float(readsum[seqsum_indexes["mad_template"]]),
                ),
                (header_to_tag["end_reason"], readsum[seqsum_indexes["end_reason"]]),
            ]
            out_bam_file.write(aln)
        # exit()

    in_bam_file.close()
    out_bam_file.close()


file_json = sys.argv[1]
# file_seqsum = '/media/dami/cyclomics_003/ont/001_accuracy_testing/SS_220209_cyclomics/002/20220209_1609_X5_FAS04073_a6645565/sequencing_summary_FAS04073_77432eb4.txt'
file_bam = sys.argv[2]
# file_bam = '/media/dami/cyclomics_003/tmp/annotationtest/Minimap2Align/SamToBam/FAS04073_pass_77432eb4_140.fastq.bam'
file_out = sys.argv[3]
# file_out = 'testing'

add_y_tags(file_json, file_bam, file_out)

pysam.sort("-o", file_out + ".sort.bam", file_out + ".bam")
pysam.index(file_out + ".sort.bam")
