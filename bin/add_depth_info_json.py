#!/usr/bin/env python
import argparse
import json


def main(global_json, depth_json, output):
    # global_json = "/home/dami/Software/cycloseq/work/22/eedc59dc73f26558e7de53fdda10d1/global.json"
    # depth_json = "/home/dami/Software/cycloseq/work/d5/9b57d8f8623d6c9546d03199098189/AIG157_pass_dab2ad05_8.json"
    # output = 'test.json'

    with open(global_json) as json_file:
        global_information = json.load(json_file)
    with open(depth_json) as json_file:
        depth = json.load(json_file)

    # add the data
    global_information["depth"] = depth
    with open(output, "w") as json_file:
        json.dump(global_information, json_file)

    print("done")


if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--global_json", type=str)
    argparser.add_argument("--depth_json", type=str)
    argparser.add_argument("--output", type=str)
    args = argparser.parse_args()
    main(args.global_json, args.depth_json, args.output)
