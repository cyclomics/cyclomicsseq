#!/usr/bin/env python

import glob
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotting_defaults import plotly_components

TAB_PRIORITY = 92

concat_type_colors = {
    "BB-I": "DodgerBlue",
    "BB-only": "Crimson",
    "I-only": "MediumSlateBlue",
    "Unmapped": "Gray",
    "mBB-only": "Tomato",
    "mI-only": "DarkOrchid",
    "BB-mI": "RoyalBlue",
    "mBB-I": "Navy",
    "mBB-mI": "SkyBlue",
    "Unknown": "Gold",
}

cycas_class_mapper = {
    "BackboneInsert": concat_type_colors["BB-I"],
    "SingleBackbone": concat_type_colors["BB-only"],
    "Unknown": concat_type_colors["Unknown"],
    "SingleInsert": concat_type_colors["I-only"],
    "SingleInsertUncertain": "Indigo",
    "DoubleInsertUncertain": "MediumOrchid",
    "BackboneDoubleInsert": concat_type_colors["BB-mI"],
    "MessyAlignment": concat_type_colors["Unknown"],
    "SingleInsertUnalignedGaps": "SlateBlue",
    "SingleBackboneUnalignedGaps": "Red",
    "LowAlignmentCount": "LightGray",
}

cycas_structure_mapper = {
    "Unknown": "Grey",
    "too_few_inserts": "LightGray",
    "1D": "DodgerBlue",
    "2D": "Gold",
    "3D": "MediumSeaGreen",
}


def sum_overlaps(structure_str: str | None) -> int:
    """Sum absolute values of overlap lengths (values preceding 'O') in structure string."""
    if not isinstance(structure_str, str) or not structure_str:
        return 0

    total = 0
    parts = structure_str.split(",")
    for part in parts:
        match = re.match(r"(-?\d+):O", part)
        if match:
            total += abs(int(match.group(1)))
    return total


def _plot_donut(classif_data, figtitle: str):
    df = pd.DataFrame.from_dict(Counter(classif_data), orient="index").reset_index()
    df.columns = ["type", "count"]
    df["color"] = df["type"].map(lambda x: cycas_class_mapper.get(x, "Grey"))

    fig = go.Figure(
        go.Pie(
            labels=df["type"],
            values=df["count"],
            marker=dict(colors=df["color"]),
            hole=0.5,
            textinfo="percent",
        )
    )

    fig.update_layout(
        title_text=figtitle,
        width=500,
        height=400,
        margin=dict(l=20, r=20, t=60, b=20),
        legend=dict(
            orientation="v",
            yanchor="middle",
            y=0.5,
            xanchor="left",
            x=1.02,
        ),
        uniformtext_minsize=10,
        uniformtext_mode="hide",
    )

    return fig


def _plot_mappability(aln_len_data, raw_len_data, figtitle: str):
    Y_n = np.divide(
        aln_len_data,
        raw_len_data,
        out=np.zeros_like(aln_len_data, dtype=float),
        where=raw_len_data != 0,
    )
    fig = go.Figure(go.Histogram(x=Y_n, nbinsx=100))
    fig.update_layout(
        title_text=figtitle,
        xaxis_title="bases/bases mapped",
        yaxis_title="Occurrence",
        width=500,
        height=400,
    )
    return fig


def _plot_lengthsegments(raw_len_data, segment_data, figtitle: str):
    fig = go.Figure(
        go.Scatter(
            x=raw_len_data,
            y=segment_data,
            mode="markers",
            marker=dict(size=8, opacity=0.5),
        )
    )
    fig.update_layout(
        title_text=figtitle,
        xaxis_title="Read length",
        yaxis_title="Aligned segments",
        width=500,
        height=400,
    )
    return fig


def _plot_segmentdist(segment_data, figtitle: str):
    fig = go.Figure(go.Histogram(x=segment_data, nbinsx=100))
    fig.update_layout(
        title_text=figtitle,
        xaxis_title="# Segments",
        yaxis_title="Occurrence",
        width=500,
        height=400,
    )
    return fig


def _plot_structures_per_read(structures, figtitle: str):
    structure_order = ["Unknown", "too_few_inserts", "1D", "2D", "3D"]

    df = (
        pd.DataFrame.from_dict(structures, orient="index")
        .rename_axis("readname")
        .reset_index()
    )
    reads_count = df["structure"].value_counts().to_dict()
    bases_count = df.groupby("structure")["raw_length"].sum().to_dict()

    total_reads = len(df)
    total_bases = df["raw_length"].sum()

    reads_pct = {k: v / total_reads * 100 for k, v in reads_count.items()}
    bases_pct = {k: v / total_bases * 100 for k, v in bases_count.items()}

    fig = go.Figure()
    for struct in structure_order:
        reads_n = reads_count.get(struct, 0)
        bases_n = bases_count.get(struct, 0)
        reads_p = reads_pct.get(struct, 0)
        bases_p = bases_pct.get(struct, 0)

        customdata = [
            [reads_n, "reads"],
            [bases_n, "bases"],
        ]

        fig.add_trace(
            go.Bar(
                name=struct,
                x=["Reads (%)", "Bases (%)"],
                y=[reads_p, bases_p],
                text=[f"{reads_p:.1f}%", f"{bases_p:.1f}%"],
                textposition="inside",
                marker_color=cycas_structure_mapper.get(struct, "Grey"),
                customdata=customdata,
                hovertemplate=(
                    "<b>%{x}</b><br>"
                    "Structure: %{fullData.name}<br>"
                    "%{customdata[0]:,} %{customdata[1]}<br>"
                    "(%{y:.1f}%)<extra></extra>"
                )
                .replace("%{x_label}", "reads")
                .replace("%{x_label}", "bases"),
            )
        )

    fig.update_layout(
        barmode="stack",
        template="simple_white",
        title_text=figtitle,
        uniformtext_minsize=10,
        uniformtext_mode="hide",
        width=500,
        height=400,
        legend_title_text="Structure",
        legend_traceorder="normal",
    )
    return fig


def _plot_structures_per_aln_count(structures, figtitle: str):
    structure_order = ["Unknown", "too_few_inserts", "1D", "2D", "3D"]
    hidden_by_default = ["Unknown", "too_few_inserts"]

    df = (
        pd.DataFrame.from_dict(structures, orient="index")
        .rename_axis("readname")
        .reset_index()
    )
    df["structure"] = pd.Categorical(
        df["structure"], categories=structure_order, ordered=True
    )

    fig = go.Figure()

    for struct in structure_order:
        subset = df[df["structure"] == struct]
        counts = subset["alignment_count"].value_counts().sort_index()
        x_vals = counts.index.values
        y_vals = counts.values / counts.sum() * 100

        fig.add_trace(
            go.Bar(
                x=x_vals,
                y=y_vals,
                name=struct,
                marker_color=cycas_structure_mapper.get(struct, "Grey"),
                customdata=counts.values,
                hovertemplate=(
                    "Alignment count: %{x}<br>"
                    "Structure: %{fullData.name}<br>"
                    "Reads: %{customdata:,} (%{y:.1f}%)<extra></extra>"
                ),
                visible="legendonly" if struct in hidden_by_default else True,
            )
        )

    fig.update_layout(
        barmode="group",
        template="simple_white",
        title_text=figtitle,
        xaxis=dict(range=[0, 30]),
        xaxis_title="Number of segments",
        yaxis_title="Percent of reads per type",
        width=500,
        height=400,
        legend_title_text="Structure",
    )

    return fig


def _plot_structures_per_length(structures, figtitle: str, bin_size=1000):
    structure_order = ["Unknown", "too_few_inserts", "1D", "2D", "3D"]
    hidden_by_default = ["Unknown", "too_few_inserts"]

    df = (
        pd.DataFrame.from_dict(structures, orient="index")
        .rename_axis("readname")
        .reset_index()
    )
    df["structure"] = pd.Categorical(
        df["structure"], categories=structure_order, ordered=True
    )

    fig = go.Figure()

    # Define bins
    max_len = df["raw_length"].max() if not df.empty else 20_000
    bins = np.arange(1, max_len + bin_size, bin_size)

    for struct in structure_order:
        subset = df[df["structure"] == struct]
        if subset.empty:
            continue

        # Compute counts per bin
        counts, edges = np.histogram(subset["raw_length"], bins=bins)
        y_vals = counts / counts.sum() * 100
        x_vals = edges[:-1]
        bin_upper = edges[1:] - 1

        fig.add_trace(
            go.Bar(
                x=x_vals,
                y=y_vals,
                width=bin_size,
                name=struct,
                marker_color=cycas_structure_mapper.get(struct, "Grey"),
                customdata=np.stack([counts, bin_upper], axis=-1),
                opacity=0.7,
                hovertemplate=(
                    "Read length: %{x} - %{customdata[1]}<br>"
                    "Structure: %{fullData.name}<br>"
                    "Reads: %{customdata[0]:,} (%{y:.1f}%)<extra></extra>"
                ),
                visible="legendonly" if struct in hidden_by_default else True,
            )
        )

    fig.update_layout(
        barmode="overlay",
        template="simple_white",
        title_text=figtitle,
        xaxis_title="Read length",
        yaxis_title="Percent of reads per type",
        width=500,
        height=400,
        legend_title_text="Structure",
    )

    return fig


def main(json_folder, plot_file, priority_limit: int, subsample_size: int = 10000):
    tab_name = "Metadata"
    json_obj = {tab_name: {"name": tab_name, "priority": TAB_PRIORITY}}

    nr_reads = 0
    raw_lens = []
    aln_lens = []
    segments = []
    classifs = []
    structures = {}

    # Load JSON metadata
    for data_json in glob.glob(f"{json_folder}/*.json"):
        with open(str(data_json)) as d:
            dict_data = json.load(d)
        for read, data in dict_data.items():
            raw_lens.append(data.get("raw_length", 0))
            total_aligned_bp = data.get("total_aligned_bases", 0)
            overlap_length = sum_overlaps(data.get("original_structure", None))
            aln_lens.append(total_aligned_bp - overlap_length)
            segments.append(data.get("alignment_count", 0))
            classifs.append(data.get("classification", "Unknown"))
            structures[read] = {
                "raw_length": data.get("raw_length", 0),
                "structure": data.get("structure_classification", "Unknown"),
                "alignment_count": data.get("alignment_count", 0),
            }
            nr_reads += 1
            if subsample_size != 0 and nr_reads >= subsample_size:
                break
        if subsample_size != 0 and nr_reads >= subsample_size:
            break

    if nr_reads == 0:
        json_obj[tab_name]["script"] = []
        json_obj[tab_name]["div"] = ["<h1>No metadata found.</h1>"]

    else:
        p1 = _plot_structures_per_read(structures, f"Read structures (n={nr_reads})")
        p2 = _plot_structures_per_aln_count(
            structures, f"Read structures (n={nr_reads})"
        )
        p3 = _plot_structures_per_length(structures, f"Read structures (n={nr_reads})")
        p4 = _plot_donut(classifs, f"Per read classification (n={nr_reads})")
        p5 = _plot_mappability(
            np.array(aln_lens),
            np.array(raw_lens),
            f"Normalized Mappability (n={nr_reads})",
        )
        p6 = _plot_lengthsegments(
            np.array(raw_lens), np.array(segments), f"Length vs segments (n={nr_reads})"
        )

        scripts, divs = plotly_components([p1, p2, p3, p4, p5, p6])
        json_obj[tab_name]["script"] = scripts
        json_obj[tab_name]["div"] = divs

    with open(Path(plot_file).with_suffix(".json"), "w", encoding="utf-8") as f:
        json.dump(json_obj, f, ensure_ascii=False)


if __name__ == "__main__":
    import argparse

    dev = True
    if dev:
        main(
            "./plot_metadata/",
            "./test_report_data/metadata",
            priority_limit=9999,
            subsample_size=1000,
        )
    else:
        parser = argparse.ArgumentParser(description="Plot read metadata statistics.")
        parser.add_argument("json_glob_path")
        parser.add_argument("plot_file")
        parser.add_argument("priority_limit", type=int, default=89)
        parser.add_argument("subsample_size", type=int, default=10000)
        args = parser.parse_args()
        main(
            args.json_glob_path,
            args.plot_file,
            args.priority_limit,
            args.subsample_size,
        )
