#!/usr/bin/env python

import argparse
import logging
from collections import defaultdict
from email.policy import default
from pathlib import Path
from typing import Dict, Tuple

import pysam

HEADER_TO_TAG = defaultdict(lambda: -1)

HEADER_TO_TAG["channel"] = "XC"
HEADER_TO_TAG["mux"] = "XM"
HEADER_TO_TAG["start_time"] = "XT"
HEADER_TO_TAG["duration"] = "XD"
HEADER_TO_TAG["adapter_duration"] = "XA"
HEADER_TO_TAG["sequence_length_template"] = "XL"
HEADER_TO_TAG["mean_qscore_template"] = "XQ"
HEADER_TO_TAG["median_template"] = "XT"
HEADER_TO_TAG["mad_template"] = "XM"
HEADER_TO_TAG["end_reason"] = "XE"

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


import json

def create_metadata(json_path):
    with open(json_path) as d:
        dict_data = json.load(d)
    return dict_data

def apply_metadata_tags(bam:str, metadata:Dict):
    pass

def extract_segment_from_meta(general_meta, full_name):
    c_reads = general_meta['consensus_reads']
    segment_data = [(k,v) for k,v in c_reads.items()]
    result = None
    for i in segment_data:

        # first check presence to prevent keyerror
        if 'readname_unique' in i[1] and i[1]['readname_unique'] == full_name:
            result = i

    return result

def extract_barcode(meta):
    return "NNNN"

def extract_partner_locations(general_meta, joiner = "|"):
    locations = []
    for k,v in general_meta['consensus_reads'].items():
        locations.append(v['alignment_position'])
    
    print(locations)
    # convert the list to a string 
    return joiner.join(locations)

def update_tags(aln, items):
    new_tags = []
    for tag, value in items.items():
        new_tags.append((tag, value))
    aln.tags = aln.tags + new_tags
    return aln

def main(metadata_json, in_bam_path, out_bam_path, split_queryname=True):
    metadata = create_metadata(metadata_json)
    in_bam = pysam.AlignmentFile(in_bam_path, "rb")

    # apply_metadata_tags(in_bam_path, metadata)
    for aln in in_bam.fetch():
        full_name = aln.query_name
        if full_name == '0e2f0d9f-2e6d-4d37-b107-93b3092eb39e_I_0_chr17:7580007:7580156':
            print('holdup')
        query_name = process_query_name(full_name, split_queryname, "_")
        gen_meta = metadata[query_name]
        seg_id, seg_meta = extract_segment_from_meta(gen_meta, full_name)

        tags = {}

        tags['YC'] = gen_meta['classification']
        tags['YT'] = gen_meta['raw_length']
        tags['YB'] = extract_barcode(gen_meta)
        tags['YP'] = extract_partner_locations(gen_meta)
        # add more info if we find the segment data
        if seg_id:
            tags['YI'] = seg_id
    
        if seg_meta:
            tags['Yt'] = seg_meta['aligned_bases_before_consensus']
            tags['YL'] = seg_meta['len']
            tags['YA'] = seg_meta['alignment_position']
            tags['YR'] = seg_meta['alignment_orientation']
            tags['YM'] = seg_meta['alignment_count']

        aln = update_tags(aln,tags)
        print(full_name)
        print(tags)
        # print(query_name)
        # print(seg_id)
        print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process the information in the sequencing summary and add it to the bam."
    )

    # parser.add_argument("file_metadata", type=Path)
    # parser.add_argument("file_bam", type=Path)
    # parser.add_argument("file_out", type=Path)
    # args = parser.parse_args()
    # # logging.info(args)
    # intermediate_bam = args.file_bam.with_suffix(".unsorted.bam")
    # # logging.info(intermediate_bam)

    # pore_time_aln = main(args.file_seqsum, args.file_bam, intermediate_bam)

    # pysam requires pure strings
    # pysam.sort("-o", str(args.file_out), str(intermediate_bam))
    # pysam.index(str(args.file_out))

    test_metadata = '/home/dami/Software/cycloseq/FAT55666_pass_8a93c5bd_123_filtered.metadata.json'
    test_bam = '/home/dami/Software/cycloseq/FAT55666_pass_8a93c5bd_123_filtered.bam'

    main(test_metadata,test_bam, 'annotatetest.bam')
    # file_seqsum = Path('/media/dami/cyclomics_003/ont/001_accuracy_testing/SS_220209_cyclomics/002/20220209_1609_X5_FAS04073_a6645565/sequencing_summary_FAS04073_77432eb4.txt')
    # file_bam = Path('/media/dami/cyclomics_003/tmp/annotationtest/Minimap2Align/SamToBam/FAS04073_pass_77432eb4_140.fastq.bam')
    # file_out = Path('testing.bam')
    # intermediate_bam = file_out.with_suffix(".unsorted.bam")

    # pore_time_aln = main(file_seqsum, file_bam, intermediate_bam)
