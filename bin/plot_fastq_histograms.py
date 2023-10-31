#!/usr/bin/env python

from collections import Counter
from glob import glob
from tqdm import tqdm
import json
from pathlib import Path
import os

from pyfastx import Fastq
import pandas as pd
import numpy as np

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource
from bokeh.layouts import column
from bokeh.embed import components

from plotting_defaults import cyclomics_defaults

TAB_PRIORITY = 6


def process_fastqs(fastqs_path):
    overall = Counter()
    lengths = []
    for fq in tqdm(glob(fastqs_path)):
        if os.stat(fq).st_size == 0:
            print(f"fastq is empty: {fq}")
            continue
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
        plot_width=cyclomics_defaults.width,
        title=my_title,
        tooltips="@Q: @count (@relative_count{%0.2f})",
    )
    p.vbar(x="Q", top="relative_count", width=0.8, source=source)
    p.y_range.start = 0
    p.x_range.start = 0

    p.title.text_font_size = "18pt"
    p.xaxis.axis_label = "Base Q"
    p.xaxis.axis_label_text_font_size = "16pt"
    p.yaxis.axis_label = "Relative abundance"
    p.yaxis.axis_label_text_font_size = "16pt"
    p.xaxis.major_label_text_font_size = "12pt"
    p.yaxis.major_label_text_font_size = "12pt"
    return p


def plot_length_hist(lengths, my_title_len):
    hist, edges = np.histogram(lengths, bins=90)
    p = figure(plot_height=500, plot_width=cyclomics_defaults.width, title=my_title_len)
    p.quad(
        top=hist,
        bottom=0,
        width=0.8,
        left=edges[:-1],
        right=edges[1:],
        line_color="white",
    )

    p.y_range.start = 0
    p.x_range.start = 0

    p.title.text_font_size = "18pt"
    p.xaxis.axis_label = "read length"
    p.xaxis.axis_label_text_font_size = "16pt"
    p.yaxis.axis_label = "Count"
    p.yaxis.axis_label_text_font_size = "16pt"
    p.xaxis.major_label_text_font_size = "12pt"
    p.yaxis.major_label_text_font_size = "12pt"

    return p


def main(
    file_extention,
    output_file_name,
    my_title_q,
    my_title_len,
    tab_name,
    priority_limit: int,
):
    json_obj = {}
    json_obj[tab_name] = {}
    json_obj[tab_name]["name"] = tab_name
    json_obj[tab_name]["priority"] = TAB_PRIORITY

    if TAB_PRIORITY < priority_limit:
        Q_scores, lengths = process_fastqs(f"*.{file_extention}")
        total = sum(Q_scores.values())
        overall_Q = [[ord(k) - 33, v / total, v] for k, v in Q_scores.items()]

        q_hist = plot_overall_Q_hist(overall_Q, my_title_q)

        len_hist = plot_length_hist(lengths, my_title_len)

        final_plot = column(q_hist, len_hist)

        json_obj[tab_name]["script"], json_obj[tab_name]["div"] = components(final_plot)
        json_obj["additional_info"] = {f"reads{tab_name}": len(lengths)}

    with open(Path(output_file_name).with_suffix(".json"), "w") as f:
        f.write(json.dumps(json_obj))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Create hist plot from a regex for fastq and fastq.gz files."
    )

    parser.add_argument("fastq_regex_suffix")
    parser.add_argument("plot_file")
    parser.add_argument("tab_name", default="fastq information")
    parser.add_argument("priority_limit", type=int, default=89)
    args = parser.parse_args()

    my_title_q = f"Q scores relative abundance"
    my_title_len = f"Length distribution"
    main(
        args.fastq_regex_suffix,
        args.plot_file,
        my_title_q,
        my_title_len,
        args.tab_name,
        args.priority_limit,
    )
