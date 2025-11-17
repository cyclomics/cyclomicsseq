#!/usr/bin/env python
import json
import os
from collections import Counter
from glob import glob
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotting_defaults import plotly_components
from pyfastx import Fastq
from tqdm import tqdm

TAB_PRIORITY = 6


def process_fastqs(fastqs_path):
    overall = Counter()
    lengths = []
    for fq in tqdm(glob(fastqs_path)):
        if os.stat(fq).st_size == 0:
            continue
        read_file = Fastq(fq, build_index=False)
        for read in read_file:
            overall.update(Counter(read[2]))
            lengths.append(len(read[1]))
    return overall, lengths


def plot_overall_Q_hist(overall_Q, figtitle: str):
    df = pd.DataFrame(overall_Q, columns=["Q", "relative_count", "count"])
    fig = go.Figure(
        go.Bar(
            x=df["Q"],
            y=df["relative_count"],
            text=[
                f"count={c}<br>rel={r:.2f}"
                for c, r in zip(df["count"], df["relative_count"])
            ],
            hoverinfo="text",
            marker=dict(color="mediumseagreen"),
        )
    )

    fig.update_layout(
        title_text=figtitle,
        xaxis_title="Base Q",
        yaxis_title="Relative abundance",
        width=600,
        height=400,
        margin=dict(l=50, r=20, t=60, b=40),
    )
    return fig


def plot_length_hist(lengths, figtitle: str):
    df = pd.DataFrame({"length": lengths})
    fig = go.Figure(
        go.Histogram(
            x=df["length"],
            nbinsx=90,
            marker=dict(color="steelblue"),
            hovertemplate="Length: %{x}<br>Count: %{y}<extra></extra>",
        )
    )

    fig.update_layout(
        title_text=figtitle,
        xaxis_title="Read length",
        yaxis_title="Count",
        width=600,
        height=400,
        margin=dict(l=40, r=20, t=60, b=40),
    )
    return fig


def main(
    file_extention,
    plot_file,
    tab_name,
    priority_limit: int,
):
    json_obj = {tab_name: {"name": tab_name, "priority": TAB_PRIORITY}}

    if TAB_PRIORITY < priority_limit:
        Q_scores, lengths = process_fastqs(f"*.{file_extention}")

        if len(lengths) > 0:
            total = sum(Q_scores.values())
            overall_Q = [[ord(k) - 33, v / total, v] for k, v in Q_scores.items()]

            p1 = plot_overall_Q_hist(overall_Q, "Q scores relative abundance")
            p2 = plot_length_hist(lengths, "Length distribution")

            scripts, divs = plotly_components([p1, p2])
            json_obj[tab_name]["script"] = scripts
            json_obj[tab_name]["div"] = divs

        else:
            json_obj[tab_name]["script"] = []
            json_obj[tab_name]["div"] = ["<h1>No reads found.</h1>"]

        json_obj["additional_info"] = {f"reads{tab_name}": len(lengths)}

    with open(Path(plot_file).with_suffix(".json"), "w", encoding="utf-8") as f:
        json.dump(json_obj, f, ensure_ascii=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Create hist plot from fastq files.")
    parser.add_argument("fastq_regex_suffix")
    parser.add_argument("plot_file")
    parser.add_argument("tab_name", default="fastq information")
    parser.add_argument("priority_limit", type=int, default=89)
    args = parser.parse_args()

    main(
        args.fastq_regex_suffix,
        args.plot_file,
        args.tab_name,
        args.priority_limit,
    )
