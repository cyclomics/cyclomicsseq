#!/usr/bin/env python

import glob
import json
import numpy as np

from bokeh.plotting import figure, show
from bokeh.layouts import row, column
from bokeh.io import save, output_file



def read_jsons_into_plots(json_folder, plot_file):
    alignment_ratio = []
    repeat_data = []
    raw_lens = []
    segments = []

    for test_json in glob.glob(f"{json_folder}/*.json"):
        with open(test_json) as d:
            dict_data = json.load(d)
            for read in dict_data.keys():
                data = dict_data[read]
                raw_len = data['raw_length']
                aln_len = sum([v['aligned_bases_before_consensus'] for k,v in data['consensus_reads'].items()])
                alignment_ratio.append((raw_len, aln_len))


                raw_len = data['raw_length']
                segment = data ['alignment_count'] 
                repeat_data.append((raw_len, segment))
                raw_lens.append(raw_len)
                segments.append(segment)


    X = [x[0] for x in alignment_ratio]
    Y = [x[1] for x in alignment_ratio]
    Y_n = [(x[1])/x[0] for x in alignment_ratio]
    hist, edges = np.histogram(Y_n, density=True, bins=100)

    p1 = figure(title="Normalized Mappability")
    p1.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], line_color="white")

    p2 = figure()
    p2.scatter(repeat_data )

    output_file(plot_file, title="metadata plots")
    save(column([p1,p2]))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Create hist plot from a regex for fastq and fastq.gz files."
    )

    parser.add_argument("json_glob_path")
    parser.add_argument("plot_file")
    args = parser.parse_args()

    read_jsons_into_plots(args.json_glob_path, args.plot_file)
