#!/usr/bin/env python

import argparse
from pathlib import Path
from datetime import datetime
from typing import Any, Dict
import json

from tqdm import tqdm
import pysam


def process_query_name(name, split_queryname=True, splitter="_"):
    """
    Split a string called name if split_queryname is True based on the splitter provided
    eg:

    abc_def with splitter _ will become abc
    abc_def_efg with splitter _ will become abc
    """
    if not split_queryname:
        return name

    return name.split(splitter)[0]


def create_metadata(json_path) -> Dict:
    """
    Open the json path and load into a dict
    """
    with open(json_path) as d:
        dict_data = json.load(d)
    return dict_data


def extract_segment_from_meta(general_meta, full_name) -> Dict:
    """
    Extract the data relevant to the specific read from the all the data available about the raw read.
    """
    c_reads = general_meta["consensus_reads"]
    segment_data = [(k, v) for k, v in c_reads.items()]
    result = None
    for i in segment_data:

        # first check presence to prevent keyerror
        if "readname_unique" in i[1] and i[1]["readname_unique"] == full_name:
            result = i

    return result


def extract_barcode(general_meta) -> str:
    """
    Extract the relevant barcode from the overall metadata
    """
    c_reads = general_meta["consensus_reads"]
    segment_data = [(k, v) for k, v in c_reads.items()]
    barcodes_full = [x[1]["barcode"] for x in segment_data]
    barcodes_filtered = [x for x in barcodes_full if x != "NNNN"]

    if len(barcodes_full) == 1:
        barcode = barcodes_full[0]
    elif len(barcodes_filtered) == 1:
        barcode = barcodes_filtered[0]
    elif len(barcodes_filtered) > 2:
        barcode = "|".join(barcodes_filtered)
    else:
        barcode = "NNNN"

    return barcode


def extract_partner_locations(general_meta, joiner="|"):
    """
    Find the location of all subreads for a given raw read
    """
    locations = []
    for k, v in general_meta["consensus_reads"].items():
        locations.append(v["alignment_position"])

    # convert the list to a string
    return joiner.join(locations)


def update_tags(
    aln: pysam.AlignedSegment, items: Dict[str, Any]
) -> pysam.AlignedSegment:
    """
    Update the BAM tags on a read.
    """
    new_tags = []
    for tag, value in items.items():
        new_tags.append((tag, value))
    aln.tags = aln.tags + new_tags
    return aln


def make_tags(gen_meta, seg_meta, seg_id) -> Dict:
    """
    Given the metadata, create a dict with the tags as key and the corresponding value as value.
    """
    tags = {}

    tags["YC"] = gen_meta["classification"]
    tags["YT"] = gen_meta["raw_length"]
    tags["YB"] = extract_barcode(gen_meta)
    tags["YP"] = extract_partner_locations(gen_meta)
    # add more info if we find the segment data
    if seg_id:
        tags["YI"] = seg_id

    if seg_meta:
        tags["YE"] = seg_meta["aligned_bases_before_consensus"]
        tags["YL"] = seg_meta["len"]
        tags["YA"] = seg_meta["alignment_position"]
        tags["YR"] = seg_meta["alignment_orientation"]
        tags["YM"] = seg_meta["alignment_count"]

    return tags


def main(metadata_json, in_bam_path, out_bam_path, split_queryname=True):
    """
    Look up the Y tags in the metadata for all reads in a bam.
    """
    start_time = datetime.now()

    metadata = create_metadata(metadata_json)
    in_bam = pysam.AlignmentFile(in_bam_path, "rb")
    out_bam_file = pysam.AlignmentFile(out_bam_path, "wb", header=in_bam.header)
    aln_count = 0
    metadata_count = 0
    for aln in tqdm(in_bam.fetch()):
        aln_count += 1
        full_name = aln.query_name
        query_name = process_query_name(full_name, split_queryname, "_")

        if query_name in metadata:
            metadata_count += 1
            gen_meta = metadata[query_name]
            seg_id, seg_meta = extract_segment_from_meta(gen_meta, full_name)
            tags = make_tags(gen_meta, seg_meta, seg_id)
            aln = update_tags(aln, tags)
        out_bam_file.write(aln)

    out_bam_file.close()
    print(
        f"Added metadata to {metadata_count} alignments out of {aln_count} in {datetime.now() - start_time}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process the information in the metadata and add it to the bam."
    )

    parser.add_argument("file_metadata", type=Path)
    parser.add_argument("file_bam", type=Path)
    parser.add_argument("file_out", type=Path)
    args = parser.parse_args()
    intermediate_bam = args.file_bam.with_suffix(".unsorted.bam")

    main(args.file_metadata, args.file_bam, intermediate_bam)

    # pysam requires pure strings
    pysam.sort("-o", str(args.file_out), str(intermediate_bam))
    pysam.index(str(args.file_out))

    # test_metadata = (
    #     "/home/dami/Software/cycloseq/FAT55692_pass_43d236d5_123_filtered.metadata.json"
    # )
    # test_bam = "/home/dami/Software/cycloseq/FAT55692_pass_43d236d5_123_filtered.annotated.bam"

    # main(test_metadata, test_bam, "annotatetest.bam")
