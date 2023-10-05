#!/usr/bin/env python

import logging
import glob
import json
import random
from collections import Counter
from math import pi
from pathlib import Path

import numpy as np
import pandas as pd
from bokeh.embed import components
from bokeh.io import output_file, save
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, HoverTool, LabelSet
from bokeh.plotting import figure
from bokeh.transform import cumsum
from plotting_defaults import cyclomics_defaults
from concurrent.futures import ProcessPoolExecutor, wait, ALL_COMPLETED
from threading import Lock

raw_len_dt = np.int32
aln_len_dt = np.int32
segment_dt = np.int16
class_dt = "U50"

TAB_PRIORITY = 92

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
    "SingleInsertUncertain": "Indigo",
    "DoubleInsertUncertain": "MediumOrchid",
    "BackboneDoubleInsert": concat_type_colors["BB-mI"],
    "MessyAlignment": concat_type_colors["Unknown"],
    "SingleInsertUnalignedGaps": "SlateBlue",
    "SingleBackboneUnalignedGaps": "Red",
    "LowAlignmentCount": "LightGray",
}


def _calculate_angle_and_color(stats, concat_type_colors):
    # print(stats)
    _df1 = pd.DataFrame.from_dict(stats, orient="index").reset_index()
    _df1.columns = ["type", "count"]
    _df1["angle"] = _df1["count"] / _df1["count"].sum() * 2 * pi
    _df1["percentage"] = _df1["count"] / _df1["count"].sum() * 100
    _df1["color"] = _df1.type.map(lambda x: concat_type_colors[x])

    return _df1


def _plot_donut(
    classif_data, figtitle: str = "Per read classification based on consensus caller"
):
    # Plot classifications donut
    classification_count = Counter(classif_data)
    data = _calculate_angle_and_color(classification_count, cycas_class_mapper)

    donut_plot_height = 500
    donut_plot_width = cyclomics_defaults.width
    donut_plot_x_range = (-0.6, 1.4)
    donut_plot_y_range = (0, 2)

    # create empty plot object
    plt = figure(
        plot_height=donut_plot_height,
        plot_width=donut_plot_width,
        title=figtitle,
        tools="hover",
        tooltips=[("type:", "@type"), ("count:", "@count"), ("%", "@percentage")],
        x_range=donut_plot_x_range,
        y_range=donut_plot_y_range,
        toolbar_location=None,
    )
    # add wedges per type
    plt.annular_wedge(
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

    plt.add_layout(labels)
    # remove chart elements
    plt.axis.axis_label = None
    plt.axis.visible = False
    plt.grid.grid_line_color = None
    plt.title.text_font_size = "16pt"

    return plt


def _plot_mappability(
    aln_len_data, raw_len_data, figtitle: str = "Normalized Mappability"
):
    # Plot mappability
    Y_n = np.divide(aln_len_data, raw_len_data)
    hist, edges = np.histogram(Y_n, density=True, bins=100)

    plt = figure(
        plot_height=500,
        plot_width=cyclomics_defaults.width,
        title=figtitle,
    )
    plt.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], line_color="white")

    plt.title.text_font_size = "18pt"
    plt.xaxis.axis_label = "bases/bases mapped"
    plt.xaxis.axis_label_text_font_size = "16pt"
    plt.yaxis.axis_label = "Occurence"
    plt.yaxis.axis_label_text_font_size = "16pt"
    plt.xaxis.major_label_text_font_size = "12pt"
    plt.yaxis.major_label_text_font_size = "12pt"
    hover = HoverTool(tooltips=[("From", "@left"), ("Until", "@right")])
    plt.add_tools(hover)

    return plt


def _plot_lengthsegments(
    raw_len_data, segment_data, figtitle: str = "Length vs segments identified"
):
    # Plot lengths vs nr segments
    plt = figure(
        plot_height=500,
        plot_width=cyclomics_defaults.width,
        title=figtitle,
    )

    plt.scatter(x=raw_len_data, y=segment_data)

    plt.title.text_font_size = "18pt"
    plt.xaxis.axis_label = "read length"
    plt.xaxis.axis_label_text_font_size = "16pt"
    plt.yaxis.axis_label = "alinged segments"
    plt.yaxis.axis_label_text_font_size = "16pt"
    plt.xaxis.major_label_text_font_size = "12pt"
    plt.yaxis.major_label_text_font_size = "12pt"
    hover = HoverTool(tooltips=[("Segments", "@y"), ("length", "@x")])
    plt.add_tools(hover)

    return plt


def _plot_segmentdist(
    segment_data, figtitle: str = "Distribution of Segments identified"
):
    # Plot distribution of segments
    my_title_segment_hist = figtitle

    int_centered_bins = np.arange(0, min(100, max(segment_data)) + 1.5) - 0.5
    hist, edges = np.histogram(segment_data, density=True, bins=int_centered_bins)
    plt = figure(
        plot_height=500,
        plot_width=cyclomics_defaults.width,
        title=my_title_segment_hist,
    )
    plt.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], line_color="white")

    plt.y_range.start = 0
    plt.x_range.start = 0
    plt.title.text_font_size = "18pt"
    plt.xaxis.axis_label = "# Segments"
    plt.xaxis.axis_label_text_font_size = "16pt"
    plt.yaxis.axis_label = "Occurence"
    plt.yaxis.axis_label_text_font_size = "16pt"
    plt.xaxis.major_label_text_font_size = "12pt"
    plt.yaxis.major_label_text_font_size = "12pt"
    hover = HoverTool(
        tooltips=[("Value", "@top{%0.4f}"), ("From", "@left"), ("Until", "@right")]
    )
    plt.add_tools(hover)

    return plt


def parse_Tidehunter_metadata(
    dict_data,
    raw_len_arr=np.array([], dtype=raw_len_dt),
    aln_len_arr=np.array([], dtype=aln_len_dt),
    segment_arr=np.array([], dtype=segment_dt),
    classif_arr=np.array([], dtype=class_dt),
):
    for i, data in enumerate(dict_data):
        if i == 0:
            continue
        raw_len = data["raw_length"]
        aln_len = float(data["baseunit_copies"]) * float(data["baseunit_length"])
        segment = data["baseunit_copies"]
        classification = "Unknown"

        raw_len_arr = np.append(raw_len_arr, np.array([raw_len], dtype=raw_len_dt))
        aln_len_arr = np.append(aln_len_arr, np.array([aln_len], dtype=aln_len_dt))
        segment_arr = np.append(segment_arr, np.array([segment], dtype=segment_dt))
        classif_arr = np.append(classif_arr, np.array([classification], dtype=class_dt))

    return raw_len_arr, aln_len_arr, segment_arr, classif_arr


def parse_Cycas_metadata(
    dict_data,
    raw_len_arr=np.array([], dtype=raw_len_dt),
    aln_len_arr=np.array([], dtype=aln_len_dt),
    segment_arr=np.array([], dtype=segment_dt),
    classif_arr=np.array([], dtype=class_dt),
):
    for read in dict_data.keys():
        data = dict_data[read]
        raw_len = data["raw_length"]
        aln_len = sum(
            [
                v["aligned_bases_before_consensus"]
                for k, v in data["consensus_reads"].items()
            ]
        )
        segment = data["alignment_count"]
        classification = data["classification"]

        raw_len_arr = np.append(raw_len_arr, np.array([raw_len], dtype=raw_len_dt))
        aln_len_arr = np.append(aln_len_arr, np.array([aln_len], dtype=aln_len_dt))
        segment_arr = np.append(segment_arr, np.array([segment], dtype=segment_dt))
        classif_arr = np.append(classif_arr, np.array([classification], dtype=class_dt))

    return raw_len_arr, aln_len_arr, segment_arr, classif_arr


def read_jsons_into_plots(
    json_folder,
    plot_file,
    subsample_files: bool = False,
    file_subset_size: int = 200,
    subsample_reads: bool = True,
    read_subset_size: int = 10_000,
    seed: int = 42,
    threads: int = 16,
):
    promises = []
    rng = np.random.default_rng(seed)

    # Parse JSON metadata with a process pool
    with ProcessPoolExecutor(max_workers=threads) as executor:
        logging.debug("creating location promisses in ProcessPoolExecutor")
        json_files = glob.glob(f"{json_folder}/*.json")
        if subsample_files:
            json_files = rng.choice(json_files, file_subset_size)

        for test_json in json_files:
            with open(str(test_json)) as d:
                dict_data_json = json.load(d)
                logging.debug("appending to ProcessPoolExecutor")
                if type(dict_data_json) == dict:
                    # Parse Cycas
                    promises.append(
                        executor.submit(parse_Cycas_metadata, dict_data_json)
                    )
                else:
                    # Parse Tidehunter
                    promises.append(
                        executor.submit(parse_Tidehunter_metadata, dict_data_json)
                    )

        logging.debug("Asking for results from ProcessPoolExecutor")
        done, not_done = wait(promises, return_when=ALL_COMPLETED)
        results = np.hstack([x.result() for x in promises])

        # Separate results into lists of correct np.dtype
        raw_lens = results[0].astype(raw_len_dt)
        aln_lens = results[1].astype(aln_len_dt)
        segments = results[2].astype(segment_dt)
        classifs = results[3].astype(class_dt)

    # Initialize HTML tab
    tab_name = "Metadata"
    json_obj = {}
    json_obj[tab_name] = {}
    json_obj[tab_name]["name"] = tab_name
    json_obj[tab_name]["priority"] = TAB_PRIORITY
    # Write out empty file if there are no metadata
    if raw_lens.size == 0:
        f = open(plot_file, "w")
        f.write("<h1>No metadata found.</h1>")
        f.close()
        json_obj[tab_name]["script"], json_obj[tab_name]["div"] = (
            "",
            "<h1>No metadata found: metadata JSON files were empty.</h1>",
        )
        with open(Path(plot_file).with_suffix(".json"), "w") as f:
            f.write(json.dumps(json_obj))
        return

    # Subsample results for lighter plotting
    n_reads = len(raw_lens)
    if n_reads > read_subset_size:
        idx = rng.choice(np.arange(n_reads), read_subset_size)
        if subsample_reads:
            raw_lens = raw_lens[idx]
            aln_lens = aln_lens[idx]
            segments = segments[idx]
            classifs = classifs[idx]

            p1 = _plot_donut(
                classifs,
                figtitle=f"Per read classification based on consensus caller (n = {read_subset_size:,})",
            )
            p2 = _plot_mappability(
                aln_lens,
                raw_lens,
                figtitle=f"Normalized Mappability (n = {read_subset_size:,})",
            )
            p3 = _plot_lengthsegments(
                raw_lens,
                segments,
                figtitle=f"Length vs segments identified (n = {read_subset_size:,})",
            )
            p4 = _plot_segmentdist(
                segments,
                figtitle=f"Distribution of Segments identified (n = {read_subset_size:,})",
            )

        else:
            subset_raw_lens = raw_lens[idx]
            subset_segments = segments[idx]

            p1 = _plot_donut(classifs)
            p2 = _plot_mappability(aln_lens, raw_lens)
            p3 = _plot_lengthsegments(
                subset_raw_lens,
                subset_segments,
                figtitle=f"Length vs segments identified (n = {read_subset_size:,})",
            )
            p4 = _plot_segmentdist(segments)

    output_file(plot_file, title="metadata plots")
    final_plot = column([p1, p2, p3, p4])

    with open(Path(plot_file).with_suffix(".json"), "w") as f:
        json_obj[tab_name]["script"], json_obj[tab_name]["div"] = components(final_plot)
        f.write(json.dumps(json_obj))

    save(final_plot)


if __name__ == "__main__":
    import argparse

    # import tracemalloc

    # tracemalloc.start()

    dev = False

    if dev:
        read_jsons_into_plots(
            "/data/projects/ROD_1002_cycseq-093/plot_metadata/jsons/",
            "./metadata.html",
            subsample_files=False,
            file_subset_size=200,
            subsample_reads=True,
            read_subset_size=10_000,
            seed=42,
            threads=16,
        )

    else:
        parser = argparse.ArgumentParser(
            description="Create hist plot from a regex for fastq and fastq.gz files."
        )

        parser.add_argument("json_glob_path")
        parser.add_argument("plot_file")
        parser.add_argument("--subsample_files", default=False)
        parser.add_argument("--file_subset_size", default=200)
        parser.add_argument("--subsample_reads", default=True)
        parser.add_argument("--read_subset_size", default=10_000)
        parser.add_argument("--seed", default=42)
        parser.add_argument("--threads", default=16)
        args = parser.parse_args()

        read_jsons_into_plots(
            args.json_glob_path,
            args.plot_file,
            args.subsample_files,
            int(args.file_subset_size),
            args.subsample_reads,
            int(args.read_subset_size),
            int(args.seed),
            int(args.threads),
        )

    # snapshot = tracemalloc.take_snapshot()
    # top_stats = snapshot.statistics("lineno")

    # print("[ Top 10 ]")
    # for stat in top_stats[:10]:
    # print(stat)
