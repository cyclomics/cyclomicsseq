#!/usr/bin/env python

import json
import sys
from collections import Counter
from glob import glob

processed_reads = {}


def process_metadata_json(json_files, output):
    print(f"\t recieved the pattern: {json_files}")
    for json_file in glob(json_files):
        print(f"\t Processing {json_file}")
        # read information from file
        with open(json_file, "r") as f:
            json_data = json.load(f)
        # parse info
        for i in json_data:
            if i["id"] in processed_reads.keys():
                assert i["classification"] == processed_reads[i["id"]]

            else:
                processed_reads[i["id"]] = i["classification"]

    result = [f"{k}, {v}" for k, v in Counter(processed_reads.values()).items()]
    result.sort()
    print("writing output file")
    with open(output, "w") as f:
        [f.write(f"{x}\n") for x in result]


json_files = sys.argv[1]
output_file = file_seqsum = sys.argv[2]

# json_files = "/media/dami/cyclomics_003/results/Cyclomics/000008_onco_panel/CycasConsensus/Cycas/*.metadata.json"
# output_file = "test_008.txt"

process_metadata_json(json_files, output_file)
