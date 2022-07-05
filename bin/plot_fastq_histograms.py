#!/usr/bin/env python

from collections import Counter
from glob import glob
from tqdm import tqdm

from pyfastx import Fastq
import pandas as pd
import numpy as np

from bokeh.io import save, output_file, show
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource
from bokeh.layouts import row, column


def process_fastqs(fastqs_path):
    overall = Counter()
    lengths = []
    for fq in tqdm(glob(fastqs_path)):
        read_file = Fastq(fq, build_index=False)
        for read in read_file:
            read_q_count = Counter(read[2])
            overall.update(read_q_count)
            lengths.append(len(read[1]))
    return overall, lengths


def plot_overall_Q_hist(overall_Q, my_title):
    df = pd.DataFrame(overall_Q, columns=["Q", "relative_count", "count"])

    source = ColumnDataSource(data=df)

    p = figure(
        plot_height=500,
        plot_width=1000,
        title=my_title,
        tooltips="@Q: @count (@relative_count{%0.2f})",
    )
    p.vbar(x="Q", top="relative_count", width=0.8, source=source)
    return p


def plot_length_hist(lengths, my_title_len):
    hist, edges = np.histogram(lengths, bins=90)
    print("hi")
    p = figure(plot_height=500, plot_width=1000, title=my_title_len)
    p.quad(
        top=hist,
        bottom=0,
        width=0.8,
        left=edges[:-1],
        right=edges[1:],
        line_color="white",
    )

    p.y_range.start = 0
    # p.legend.location = "center_right"
    # p.legend.background_fill_color = "#fefefe"
    # p.xaxis.axis_label = 'x'
    # p.yaxis.axis_label = 'Pr(x)'
    # p.grid.grid_line_color="white"
    return p


def main(file_extention, output_file_name, my_title_q, my_title_len):

    Q_scores, lengths = process_fastqs(f"*.{file_extention}")
    total = sum(Q_scores.values())
    overall_Q = [[ord(k) - 33, v / total, v] for k, v in Q_scores.items()]

    print(overall_Q)
    output_file(filename=output_file_name, title="fastq information")

    q_hist = plot_overall_Q_hist(overall_Q, my_title_q)

    len_hist = plot_length_hist(lengths, my_title_len)

    save(column(q_hist, len_hist))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Create hist plot from a regex for fastq and fastq.gz files."
    )

    parser.add_argument("fastq_regex_suffix")
    parser.add_argument("plot_file")
    args = parser.parse_args()

    my_title_q = f"Q scores relative abundance "
    my_title_len = f"Length distribution"
    main(args.fastq_regex_suffix, args.plot_file, my_title_q, my_title_len)

    # fastq_file_name = "/media/dami/a2bc89fb-be6b-4e23-912a-0c7137cd69ad/raw_data/Cyclomics/000014/HC01_CG_001/20220630_1612_MN40283_FAT55621_6a2c8b02/fastq_pass/*13.fastq.gz"
    # my_title_q = f"Q scores in all Fastq files found using {fastq_file_name}"
    # my_title_len = f"Length distribution"
    # output_file_name = "qscore.html"
    # main(fastq_file_name, output_file_name, my_title_q,  my_title_len)
