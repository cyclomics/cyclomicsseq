#!/usr/bin/env python3

import argparse
import json
from collections import defaultdict
from pathlib import Path


def bed_to_json(bed):
    result = defaultdict(list)
    for reg in bed.readlines():
        chrom, pos, depth = reg.strip().split("\t")
        result[chrom].append({"position": pos, "depth": depth})
    return result


def main(args):

    bed_file = Path(args.bed)
    json_file = Path(bed_file.stem).with_suffix(".json")

    with bed_file.open(mode="r") as bed:
        result = bed_to_json(bed)

    with json_file.open(mode="w") as json_out:
        json_out.write(json.dumps(result))


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--bed", type=str)
    args = argparser.parse_args()
    main(args)
