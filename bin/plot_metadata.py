#!/usr/bin/env python

import glob
import json
from collections import Counter
from math import pi
from pathlib import Path
import random

import numpy as np
import pandas as pd

from bokeh.plotting import figure
from bokeh.layouts import column
from bokeh.io import save, output_file
from bokeh.transform import cumsum
from bokeh.models import LabelSet, ColumnDataSource, HoverTool
from bokeh.embed import components


from plotting_defaults import cyclomics_defaults

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
cycas_class_mapper = {
    "BackboneInsert": concat_type_colors["BB-I"],
    "SingleBackbone": concat_type_colors["BB-only"],
    "Unknown": concat_type_colors["Unknown"],
    "SingleInsert": concat_type_colors["I-only"],
    "BackboneDoubleInsert": concat_type_colors["BB-mI"],
    "MessyAlignment": concat_type_colors["Unknown"],
    "SingleInsertUnalignedGaps": "SlateBlue",
    "SingleBackboneUnalignedGaps": "Red",
    "LowAlignmentCount": "LightGray",
}


def _calculate_angle_and_color(stats, concat_type_colors):
    print(stats)
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


def parse_Cycas_metadata(dict_data):
    print("cycas parsing")

    alignment_ratio, repeat_data, raw_lens, segments, classifications = (
        [],
        [],
        [],
        [],
        [],
    )

    for read in dict_data.keys():
        data = dict_data[read]
        raw_len = data["raw_length"]
        aln_len = sum(
            [
                v["aligned_bases_before_consensus"]
                for k, v in data["consensus_reads"].items()
            ]
        )
        alignment_ratio.append((raw_len, aln_len))

        raw_len = data["raw_length"]
        segment = data["alignment_count"]
        classification = data["classification"]

        repeat_data.append((raw_len, segment))
        raw_lens.append(raw_len)
        segments.append(segment)
        classifications.append(classification)

    return alignment_ratio, repeat_data, raw_lens, segments, classifications


def parse_Tidehunter_metadata(dict_data):
    print("tidehunter parsing")
    alignment_ratio, repeat_data, raw_lens, segments, classifications = (
        [],
        [],
        [],
        [],
        [],
    )

    for i, data in enumerate(dict_data):
        if i == 0:
            continue
        print(data)
        aln_ratio = float(data["baseunit_copies"]) * float(data["baseunit_length"])
        alignment_ratio.append((int(data["raw_length"]), aln_ratio))
        repeat_data.append((int(data["raw_length"]), float(data["baseunit_copies"])))
        raw_lens.append(int(data["raw_length"]))
        segments.append(float(data["baseunit_copies"]))
        classifications.append("Unknown")
    print(dict_data)
    return alignment_ratio, repeat_data, raw_lens, segments, classifications


def read_jsons_into_plots(json_folder, plot_file):
    print(json_folder)
    dict_data = None
    for test_json in glob.glob(f"{json_folder}/*.json"):
        with open(test_json) as d:
            print(test_json)
            dict_data_json = json.load(d)
            # If cycas
            if type(dict_data_json) == dict:
                if not dict_data:
                    dict_data = {}
                dict_data = {**dict_data, **dict_data_json}
            # Tidehunter
            else:
                if not dict_data:
                    dict_data = []
                dict_data += dict_data_json

    if type(dict_data) == dict:
        (
            alignment_ratio,
            repeat_data,
            raw_lens,
            segments,
            classifications,
        ) = parse_Cycas_metadata(dict_data)
    else:
        (
            alignment_ratio,
            repeat_data,
            raw_lens,
            segments,
            classifications,
        ) = parse_Tidehunter_metadata(dict_data)

    X = [x[0] for x in alignment_ratio]
    Y = [x[1] for x in alignment_ratio]
    Y_n = [(x[1]) / x[0] for x in alignment_ratio]
    hist, edges = np.histogram(Y_n, density=True, bins=100)

    p1 = figure(
        plot_height=500,
        plot_width=cyclomics_defaults.width,
        title="Normalized Mappability",
    )
    p1.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], line_color="white")

    p1.title.text_font_size = "18pt"
    p1.xaxis.axis_label = "bases/bases mapped"
    p1.xaxis.axis_label_text_font_size = "16pt"
    p1.yaxis.axis_label = "Occurence"
    p1.yaxis.axis_label_text_font_size = "16pt"
    p1.xaxis.major_label_text_font_size = "12pt"
    p1.yaxis.major_label_text_font_size = "12pt"
    hover = HoverTool(tooltips=[("From", "@left"), ("Until", "@right")])
    p1.add_tools(hover)

    p2 = figure(
        plot_height=500,
        plot_width=cyclomics_defaults.width,
        title="Length vs segments identified",
    )
    if len(raw_lens) > 10_000:
        full = list(zip(raw_lens, segments))
        subset = random.sample(full, 10_000)
        raw_lens_subset = [x[0] for x in subset]
        segments_subset = [x[1] for x in subset]
    else:
        raw_lens_subset = raw_lens
        segments_subset = segments

    p2.scatter(x=raw_lens_subset, y=segments_subset)

    p2.title.text_font_size = "18pt"
    p2.xaxis.axis_label = "read length"
    p2.xaxis.axis_label_text_font_size = "16pt"
    p2.yaxis.axis_label = "alinged segments"
    p2.yaxis.axis_label_text_font_size = "16pt"
    p2.xaxis.major_label_text_font_size = "12pt"
    p2.yaxis.major_label_text_font_size = "12pt"
    hover = HoverTool(tooltips=[("Segments", "@y"), ("length", "@x")])
    p2.add_tools(hover)

    classification_count = Counter(classifications)
    _df1 = _calculate_angle_and_color(classification_count, cycas_class_mapper)

    donut_plot_height = 500
    donut_plot_width = cyclomics_defaults.width
    donut_plot_x_range = (-0.6, 1.4)
    donut_plot_y_range = (0, 2)
    donut = _plot_donut(
        _df1,
        donut_plot_height,
        donut_plot_width,
        donut_plot_x_range,
        donut_plot_y_range,
        "Per read classification based on consensus caller ",
    )

    my_title_segment_hist = "Distribution of Segments identified"

    int_centered_bins = np.arange(0, min(100, max(segments)) + 1.5) - 0.5
    hist, edges = np.histogram(segments, density=True, bins=int_centered_bins)
    density_plot = figure(
        plot_height=500,
        plot_width=cyclomics_defaults.width,
        title=my_title_segment_hist,
    )
    density_plot.quad(
        top=hist, bottom=0, left=edges[:-1], right=edges[1:], line_color="white"
    )

    density_plot.y_range.start = 0
    density_plot.x_range.start = 0
    density_plot.title.text_font_size = "18pt"
    density_plot.xaxis.axis_label = "# Segments"
    density_plot.xaxis.axis_label_text_font_size = "16pt"
    density_plot.yaxis.axis_label = "Occurence"
    density_plot.yaxis.axis_label_text_font_size = "16pt"
    density_plot.xaxis.major_label_text_font_size = "12pt"
    density_plot.yaxis.major_label_text_font_size = "12pt"
    hover = HoverTool(
        tooltips=[("Value", "@top{%0.4f}"), ("From", "@left"), ("Until", "@right")]
    )
    density_plot.add_tools(hover)

    output_file(plot_file, title="metadata plots")
    final_plot = column([donut, p1, p2, density_plot])

    tab_name = "metadata"
    with open(Path(plot_file).with_suffix(".json"), "w") as f:
        json_obj = {}
        json_obj[tab_name] = {}
        json_obj[tab_name]["name"] = tab_name
        json_obj[tab_name]["script"], json_obj[tab_name]["div"] = components(final_plot)
        f.write(json.dumps(json_obj))

    save(final_plot)


if __name__ == "__main__":
    import argparse

    dev = False

    if dev:
        read_jsons_into_plots("metadata_dev2", "metadata2.html")
    else:
        parser = argparse.ArgumentParser(
            description="Create hist plot from a regex for fastq and fastq.gz files."
        )

        parser.add_argument("json_glob_path")
        parser.add_argument("plot_file")
        args = parser.parse_args()

        read_jsons_into_plots(args.json_glob_path, args.plot_file)
    # read_jsons_into_plots('dummy_json', 'tmp.html')
