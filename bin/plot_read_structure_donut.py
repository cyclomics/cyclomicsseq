#!/usr/bin/env python

import json
from collections import Counter
from itertools import chain
from math import pi
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import pysam
from bokeh.embed import components
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Div, LabelSet
from bokeh.plotting import figure
from bokeh.transform import cumsum

TAB_PRIORITY_CONTIG_COUNT = 90
TAB_PRIORITY_DONUT = 1


def determine_read_type(
    read_info_list: List[pysam.AlignedSegment],
    per_read_counter: Counter,
    per_base_counter: Counter,
) -> None:
    """
    Given a list of pysam objects, update the per_read_counter and per_base_counter based on their infered read type from the present chromosomes.

    """
    # ignore empty reads. should not occur since the read should be unmapped.
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

    # two elements mapped
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

    # complex mapping
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
    """
    Given a list of pysam objects, update the per_read_counter and per_base_counter based on their infered read type from the present chromosomes.
    """

    mapped = list(set([x.reference_name if x.reference_name else "*" for x in reads]))

    for map in mapped:
        if map in chromosome_counts.keys():
            chromosome_counts[map] += 1
        else:
            chromosome_counts[map] = 1


def direct_to_Couter(bam) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, int]]:
    """
    Given a coordinate sorted and indexed bam, create counters for read types and chromosome counts.

    """

    concat_type_stats, concat_type_stats_by_bases = initialize_counters()
    chromosome_counts = {}
    read_info = []
    readname = None
    for i, read in enumerate(bam):
        if read.qname == readname:
            read_info.append(read)
        else:
            if not readname:
                readname = read.qname
                read_info = [read]

            determine_read_type(
                read_info, concat_type_stats, concat_type_stats_by_bases
            )
            update_chromosome_counts(read_info, chromosome_counts)

            readname = read.qname
            read_info = [read]

    return concat_type_stats, concat_type_stats_by_bases, chromosome_counts


def initialize_counters():
    concat_types = {
        "BB-I": 0,
        "BB-only": 0,
        "I-only": 0,
        "Unmapped": 0,
        "mBB-only": 0,  # multiple BB, no I
        "mI-only": 0,  # multiple I, no BB
        "BB-mI": 0,  # single BB, multiple I
        "mBB-I": 0,  # multiple BB, single I
        "mBB-mI": 0,  # multiple BB, multiple I
        "Unknown": 0,  # anything else
    }
    concat_type_stats = Counter(concat_types)
    concat_type_stats_by_bases = Counter(concat_types)
    return concat_type_stats, concat_type_stats_by_bases


def create_assembly_count_plot(chromosome_counts, my_title, priority_limit: int):
    """
    Create a barplot for the assembly statistics.
    """
    tab_name = "Alignment"
    add_info = {}
    json_obj = {}
    json_obj[tab_name] = {}
    json_obj[tab_name]["name"] = tab_name
    json_obj[tab_name]["priority"] = TAB_PRIORITY_CONTIG_COUNT

    if TAB_PRIORITY_CONTIG_COUNT < priority_limit:
        ##Segments count
        # _d = Counter(split_table.RNAME.values)
        segments = (
            pd.DataFrame.from_dict(chromosome_counts, orient="index")
            .reset_index()
            .rename(columns={"index": "CHROM", 0: "count"})
        )
        print(segments)
        source = ColumnDataSource(data=segments)

        p = figure(
            plot_height=500,
            plot_width=1000,
            x_range=segments.CHROM,
            title=my_title,
            tooltips="@CHROM: @count",
        )
        p.vbar(x="CHROM", top="count", width=0.8, source=source)
        p.xaxis.major_label_orientation = "vertical"

        json_obj[tab_name]["script"], json_obj[tab_name]["div"] = components(p)

    json_obj["additional_info"] = add_info

    return json_obj


def create_readtype_donuts(
    concat_type_stats,
    concat_type_stats_by_bases,
    plot_title,
    priority_limit: int,
):
    """
    Create a donut plot for the read structures.
    """

    concat_type_colors = {
        "BB-I": "DodgerBlue",  # perfect
        "BB-only": "Crimson",  # waste
        "I-only": "MediumSlateBlue",  # informing but missing barcode
        "Unmapped": "Gray",  # waste
        "mBB-only": "Tomato",  # waste
        "mI-only": "DarkOrchid",  # informing but missing barcode
        "BB-mI": "RoyalBlue",  # OK
        "mBB-I": "Navy",  # OK
        "mBB-mI": "SkyBlue",  # OK
        "Unknown": "Gold",  # Waste
    }

    def _calculate_angle_and_color(stats, concat_type_colors):
        _df1 = pd.DataFrame.from_dict(stats, orient="index").reset_index()
        _df1.columns = ["type", "count"]
        _df1["angle"] = _df1["count"] / _df1["count"].sum() * 2 * pi
        _df1["percentage"] = _df1["count"] / _df1["count"].sum() * 100
        _df1["color"] = _df1.type.map(lambda x: concat_type_colors[x])

        return _df1

    def _plot_donut(
        data,
        donut_plot_height,
        donut_plot_width,
        donut_plot_x_range,
        donut_plot_y_range,
        subtitle,
    ):
        # data = ColumnDataSource(data)

        # create empty plot object
        donut_plot = figure(
            plot_height=donut_plot_height,
            plot_width=donut_plot_width,
            title=subtitle,
            tools="hover",
            tooltips=[("type:", "@type"), ("count:", "@count"), ("%", "@percentage")],
            x_range=donut_plot_x_range,
            y_range=donut_plot_y_range,
            toolbar_location=None,
        )
        # add wedges per type
        donut_plot.annular_wedge(
            x=0,
            y=1,
            inner_radius=0.2,
            outer_radius=0.45,
            direction="anticlock",
            start_angle=cumsum("angle", include_zero=True),
            end_angle=cumsum("angle"),
            line_color="white",
            fill_color="color",
            legend_group="type",
            source=data,
        )
        # mask all % below 2%
        data["percentage"] = data["percentage"].mask(data["percentage"] < 2)
        # format to 1 decimal
        data["percentage"] = data["percentage"].map("{:,.1f}%".format)
        # remove nan string
        data["percentage"] = data["percentage"].replace(["nan%"], "")
        # padd to appear in right place
        data["percentage"] = data["percentage"].str.pad(18, side="left")
        source = ColumnDataSource(data)

        labels = LabelSet(
            x=0,
            y=1,
            text="percentage",
            angle=cumsum("angle", include_zero=True),
            source=source,
            render_mode="canvas",
        )

        donut_plot.add_layout(labels)

        # remove chart elements
        donut_plot.axis.axis_label = None
        donut_plot.axis.visible = False
        donut_plot.grid.grid_line_color = None
        donut_plot.title.text_font_size = "16pt"
        return donut_plot

    tab_name = "Read structure"
    json_obj = {}
    json_obj[tab_name] = {}
    json_obj[tab_name]["name"] = tab_name
    json_obj[tab_name]["priority"] = TAB_PRIORITY_DONUT

    if TAB_PRIORITY_DONUT < priority_limit:
        _df1 = _calculate_angle_and_color(concat_type_stats, concat_type_colors)
        _df2 = _calculate_angle_and_color(
            concat_type_stats_by_bases, concat_type_colors
        )

        add_info = {}
        add_info["read_struc_prec_bbi"] = _df1[_df1.type == "BB-I"].percentage[0]

        donut_plot_height = 400
        donut_plot_width = 600
        donut_plot_x_range = (-0.6, 1.4)
        donut_plot_y_range = (0, 2)
        p1 = _plot_donut(
            _df1,
            donut_plot_height,
            donut_plot_width,
            donut_plot_x_range,
            donut_plot_y_range,
            "per read",
        )
        p2 = _plot_donut(
            _df2,
            donut_plot_height,
            donut_plot_width,
            donut_plot_x_range,
            donut_plot_y_range,
            "per base",
        )

        donut_row = row(p1, p2)
        donut_plot = column(Div(text=f"<h1>{plot_title}</h1>"), donut_row)

        json_obj[tab_name]["script"], json_obj[tab_name]["div"] = components(donut_plot)
        json_obj["additional_info"] = add_info

    return json_obj


def main(
    bam_path,
    assembly_plot_title,
    plot_file,
    donut_plot_title,
    priority_limit: int,
):
    """
    Extract data from bam and plot it using two functions that further select data
    """

    aln_file = pysam.AlignmentFile(bam_path, "rb")

    read_counts, base_counts, chromosome_counts = direct_to_Couter(aln_file)

    assembly_info = create_assembly_count_plot(
        chromosome_counts, assembly_plot_title, priority_limit
    )
    donut_info = create_readtype_donuts(
        read_counts, base_counts, donut_plot_title, priority_limit
    )

    json_obj = chain(assembly_info.items(), donut_info.items())
    json_obj = dict(json_obj)
    # correct the add info by merging them
    json_obj["additional_info"] = assembly_info["additional_info"]
    json_obj["additional_info"].update(donut_info["additional_info"])

    with open(Path(plot_file).with_suffix(".json"), "w") as f:
        f.write(json.dumps(json_obj))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Create plots for a split read bam.")

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

    # split_bam = (
    #     "/scratch/nxf_work/dami/bf/b4e997120a32723705a5ffd59f7ba1/tmp_readname_sorted_splibams_merged.bam"
    # )

    # main(aln_file, "assem.html", "assem count", "Donuts.html", "Donut plot 123")
