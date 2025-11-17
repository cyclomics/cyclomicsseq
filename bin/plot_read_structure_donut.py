#!/usr/bin/env python

import json
from collections import Counter
from itertools import chain
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import plotly.graph_objects as go
import pysam
from plotting_defaults import plotly_components

TAB_PRIORITY_CONTIG_COUNT = 90
TAB_PRIORITY_DONUT = 1


def determine_read_type(
    read_info_list: List[pysam.AlignedSegment],
    per_read_counter: Counter,
    per_base_counter: Counter,
) -> None:
    if not read_info_list:
        return

    mapped = list(
        set([x.reference_name if x.is_mapped else "*" for x in read_info_list])
    )
    raw_read_length = read_info_list[0].infer_read_length()
    if not raw_read_length:
        raw_read_length = len(read_info_list[0].seq)

    mapped_bb = [x for x in mapped if x.startswith("BB")]
    mapped_ins = [x for x in mapped if not x.startswith("BB")]

    if len(mapped) == 1:
        if mapped == ["*"]:
            per_read_counter.update({"Unmapped"})
            per_base_counter.update({"Unmapped": raw_read_length})
        elif mapped[0] in mapped_bb:
            per_read_counter.update({"BB-only"})
            per_base_counter.update({"BB-only": raw_read_length})
        else:
            per_read_counter.update({"I-only"})
            per_base_counter.update({"I-only": raw_read_length})
    elif len(mapped) == 2:
        if mapped_bb and mapped_ins:
            per_read_counter.update({"BB-I"})
            per_base_counter.update({"BB-I": raw_read_length})
        elif mapped_bb and not mapped_ins:
            per_read_counter.update({"BB-only"})
            per_base_counter.update({"BB-only": raw_read_length})
        else:
            per_read_counter.update({"I-only"})
            per_base_counter.update({"I-only": raw_read_length})
    else:
        if not mapped_ins and len(mapped_bb) > 1:
            per_read_counter.update({"mBB-only"})
            per_base_counter.update({"mBB-only": raw_read_length})
        elif not mapped_bb and len(mapped_ins) > 1:
            per_read_counter.update({"mI-only"})
            per_base_counter.update({"mI-only": raw_read_length})
        elif len(mapped_ins) > 1 and len(mapped_bb) == 1:
            per_read_counter.update({"BB-mI"})
            per_base_counter.update({"BB-mI": raw_read_length})
        elif len(mapped_bb) > 1 and len(mapped_ins) == 1:
            per_read_counter.update({"mBB-I"})
            per_base_counter.update({"mBB-I": raw_read_length})
        elif len(mapped_bb) > 1 and len(mapped_ins) > 1:
            per_read_counter.update({"mBB-mI"})
            per_base_counter.update({"mBB-mI": raw_read_length})
        else:
            per_read_counter.update({"Unknown"})
            per_base_counter.update({"Unknown": raw_read_length})


def update_chromosome_counts(
    reads: List[pysam.AlignedSegment], chromosome_counts: Counter
) -> None:
    mapped = list(set([x.reference_name if x.reference_name else "*" for x in reads]))
    for map in mapped:
        chromosome_counts[map] = chromosome_counts.get(map, 0) + 1


def initialize_counters():
    concat_types = {
        "BB-I": 0,
        "BB-only": 0,
        "I-only": 0,
        "Unmapped": 0,
        "mBB-only": 0,
        "mI-only": 0,
        "BB-mI": 0,
        "mBB-I": 0,
        "mBB-mI": 0,
        "Unknown": 0,
    }
    concat_type_stats = Counter(concat_types)
    concat_type_stats_by_bases = Counter(concat_types)
    return concat_type_stats, concat_type_stats_by_bases


def direct_to_Couter(bam) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, int]]:
    concat_type_stats, concat_type_stats_by_bases = initialize_counters()
    chromosome_counts = {}
    read_info = []
    readname = None

    for read in bam:
        if read.qname == readname:
            read_info.append(read)
        else:
            if read_info:
                determine_read_type(
                    read_info, concat_type_stats, concat_type_stats_by_bases
                )
                update_chromosome_counts(read_info, chromosome_counts)
            readname = read.qname
            read_info = [read]

    if read_info:
        determine_read_type(read_info, concat_type_stats, concat_type_stats_by_bases)
        update_chromosome_counts(read_info, chromosome_counts)

    return concat_type_stats, concat_type_stats_by_bases, chromosome_counts


def create_assembly_count_plot(chromosome_counts, my_title, priority_limit: int):
    tab_name = "Alignment"
    json_obj = {
        tab_name: {"name": tab_name, "priority": TAB_PRIORITY_CONTIG_COUNT},
        "additional_info": {},
    }

    if TAB_PRIORITY_CONTIG_COUNT < priority_limit:
        df = pd.DataFrame(list(chromosome_counts.items()), columns=["CHROM", "count"])
        fig = go.Figure(
            go.Bar(
                x=df["CHROM"],
                y=df["count"],
                marker_color="steelblue",
                hovertemplate="CHROM=%{x}<br>Count=%{y}<extra></extra>",
            )
        )
        fig.update_layout(
            title=my_title,
            xaxis_title="Chromosome",
            yaxis_title="Count",
            xaxis_tickangle=-45,
            width=1000,
            height=500,
        )

        script, div = plotly_components([fig])
        json_obj[tab_name]["script"], json_obj[tab_name]["div"] = script, div

    return json_obj


def create_readtype_donuts(
    concat_type_stats,
    concat_type_stats_by_bases,
    plot_title,
    priority_limit: int,
):
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

    def _plot_donut(stats, subtitle):
        df = pd.DataFrame.from_dict(stats, orient="index").reset_index()
        df.columns = ["type", "count"]
        df["color"] = df["type"].map(lambda x: concat_type_colors.get(x, "gray"))
        total_count = df["count"].sum()

        fig = go.Figure(
            go.Pie(
                labels=df["type"],
                values=df["count"],
                hole=0.5,
                marker=dict(colors=df["color"]),
                textinfo="label+percent",
                insidetextorientation="radial",
                textposition="inside",
                hovertemplate="%{label}: %{value} (%{percent})<extra></extra>",
                # Hide labels for slices <5%
                text=[
                    f"{row['type']} {row['count']} ({row['count'] / total_count * 100:.1f}%)"
                    if row["count"] / total_count > 0.05
                    else ""
                    for _, row in df.iterrows()
                ],
            )
        )
        fig.update_layout(
            title_text=subtitle,
            width=400,
            height=400,
            margin=dict(l=20, r=20, t=60, b=20),
            legend=dict(
                orientation="v",
                yanchor="top",
                y=1.0,
                xanchor="left",
                x=1.05,
                font=dict(size=12),
            ),
        )
        return fig

    tab_name = "Read structure"
    json_obj = {tab_name: {"name": tab_name, "priority": TAB_PRIORITY_DONUT}}

    if TAB_PRIORITY_DONUT < priority_limit:
        add_info = {}
        total = sum(concat_type_stats.values())
        add_info["read_struc_prec_bbi"] = (
            concat_type_stats.get("BB-I", 0) / total * 100 if total else 0
        )

        # Generate two independent figures
        fig_per_read = _plot_donut(concat_type_stats, "Per read")
        fig_per_base = _plot_donut(concat_type_stats_by_bases, "Per base")

        # Pass both figures as independent plots to plotly_components
        script, div = plotly_components([fig_per_read, fig_per_base])
        json_obj[tab_name]["script"], json_obj[tab_name]["div"] = script, div
        json_obj["additional_info"] = add_info

    return json_obj


def main(
    bam_path,
    assembly_plot_title,
    plot_file,
    donut_plot_title,
    priority_limit: int,
):
    aln_file = pysam.AlignmentFile(bam_path, "rb")
    read_counts, base_counts, chromosome_counts = direct_to_Couter(aln_file)

    assembly_info = create_assembly_count_plot(
        chromosome_counts, assembly_plot_title, priority_limit
    )
    donut_info = create_readtype_donuts(
        read_counts, base_counts, donut_plot_title, priority_limit
    )

    json_obj = dict(chain(assembly_info.items(), donut_info.items()))
    json_obj["additional_info"] = assembly_info["additional_info"]
    json_obj["additional_info"].update(donut_info.get("additional_info", {}))

    with open(Path(plot_file).with_suffix(".json"), "w") as f:
        json.dump(json_obj, f, indent=2)


if __name__ == "__main__":
    import argparse

    dev = False
    if dev:
        import os

        os.chdir("test_read_structure")
        main(
            "ALK_multiplex_panel_AMP2_lSQDW0i_barcode06.merged.bam",
            "Segment count per assembly in reference",
            "ALK_multiplex_panel_AMP2_lSQDW0i_barcode06_read_structure.json",
            "Read structure based on chromosomal presence per readname",
            999,
        )
    else:
        parser = argparse.ArgumentParser(
            description="Create plots for a split read bam."
        )
        parser.add_argument("bam_file")
        parser.add_argument("plot_file")
        parser.add_argument("priority_limit", type=int, default=89)
        args = parser.parse_args()

        main(
            args.bam_file,
            "Segment count per assembly in reference",
            args.plot_file,
            "Read structure based on chromosomal presence per readname",
            args.priority_limit,
        )
