#!/usr/bin/env python

import io
from typing import List, Tuple, Any
import pandas as pd
import numpy as np
from pathlib import Path
import json

from bokeh.plotting import figure
from bokeh.layouts import column
from bokeh.models import HoverTool
from bokeh.embed import components

from plotting_defaults import cyclomics_defaults

chromosomal_region = Tuple[str, int, int]
TAB_PRIORITY = 3


def relaxed_float(x: Any) -> float:
    """Return a float, with value error catch"""
    try:
        my_float = float(x)
    except ValueError:
        my_float = x
    return my_float


def get_roi_pileup_df(
    df: pd.DataFrame, distance: int = 100
) -> List[chromosomal_region]:
    """
    Given a dataframe with CHROM and POS columns, find all streches in chromosomes that are less than `distance` appart.
    reports

    """

    def get_continous_strech(array: np.array, distance: int):
        """
        https://stackoverflow.com/questions/47183828/pandas-how-to-find-continuous-values-in-a-series-whose-differences-are-within-a
        """
        m = np.concatenate(([True], array[1:] > array[:-1] + distance, [True]))
        idx = np.flatnonzero(m)
        l = array.tolist()
        return [l[i:j] for i, j in zip(idx[:-1], idx[1:])]

    assemblies = df.CHROM.unique()
    roi = []
    for chrom in assemblies:
        df_tmp = df[df.CHROM == chrom]
        streches = get_continous_strech(df_tmp.POS.to_numpy(), distance)
        for strech in streches:
            roi.append((chrom, strech[0], strech[-1]))

    return roi


def read_vcf(path):
    with open(path, "r") as f:
        lines = [l for l in f if not l.startswith("##")]
    df = (
        pd.read_csv(
            io.StringIO("".join(lines)),
            dtype={
                "#CHROM": str,
                "POS": int,
                "ID": str,
                "REF": str,
                "ALT": str,
                "QUAL": str,
                "FILTER": str,
                "INFO": str,
            },
            sep="\t",
        )
        .rename(columns={"#CHROM": "CHROM"})
        .rename(columns=str.upper)
    )

    if df.empty:
        return df

    formats = df.FORMAT[0].split(":")
    for i, fmt in enumerate(formats):
        df[fmt] = df.SAMPLE1.apply(
            lambda x: relaxed_float(x.split(":")[i]) if (x.split(":")[i]) else 0
        )

    del df["FORMAT"]
    del df["SAMPLE1"]
    return df


def update_vcf_with_truth(true_vcf, df):
    true_df = read_vcf(true_vcf)
    true_df = true_df[["CHROM", "POS"]]
    true_df["Variant"] = 1
    df = pd.merge(
        df, true_df, "left", left_on=["CHROM", "POS"], right_on=["CHROM", "POS"]
    )
    df["Variant"].fillna(0, inplace=True)

    return df


def make_scatter_plots(data, roi):
    hover = HoverTool(
        tooltips=[
            ("position", "@CHROM : @POS"),
            ("alleles", "@REF -> @ALT"),
            ("VAF", "@VAF"),
            ("Frequency", "@FREQ"),
            ("forward ratio", "@FWDR"),
            ("reverse ratio", "@REVR"),
            ("Same base", "@SAME"),
            ("Depth", "@DP"),
            ("Depth after Qfilter", "@DPQ"),
        ]
    )

    plots = []
    for region in roi:
        chrom = region[0]
        start = region[1]
        stop = region[2]

        data_relevant = data.copy()
        data_relevant = data_relevant[data_relevant.CHROM == chrom]
        data_relevant = data_relevant[data_relevant.POS >= start]
        data_relevant = data_relevant[data_relevant.POS <= stop]

        data_snp = data_relevant[
            data_relevant.REF.str.len() == data_relevant.ALT.str.len()
        ]
        data_indel = data_relevant[
            data_relevant.REF.str.len() != data_relevant.ALT.str.len()
        ]

        # Pos scatter vaf
        p_vaf = figure(
            title=f"Positional VAF {chrom}:{start}-{stop}",
            width=cyclomics_defaults.width,
        )
        p_vaf.xaxis.axis_label = "Position"
        p_vaf.yaxis.axis_label = "VAF"
        p_vaf.xaxis.formatter.use_scientific = False

        p_vaf.scatter("POS", "VAF", source=data_snp, legend_label="SNP")
        p_vaf.scatter(
            "POS", "VAF", source=data_indel, color="red", legend_label="InDel"
        )

        p_vaf.add_tools(hover)
        p_vaf.title.text_font_size = "16pt"
        p_vaf.xaxis.axis_label_text_font_size = "12pt"
        p_vaf.yaxis.axis_label_text_font_size = "12pt"

        # fwd ratio vs reverse ratio
        p_ratio = figure(
            title=f"Support ratio's per direction {chrom}:{start}-{stop}",
            width=cyclomics_defaults.width,
        )
        p_ratio.xaxis.axis_label = "forward  ratio"
        p_ratio.yaxis.axis_label = "reverse ratio"

        p_ratio.scatter("FWDR", "REVR", source=data_snp, legend_label="SNP")
        p_ratio.scatter(
            "FWDR", "REVR", source=data_indel, color="red", legend_label="InDel"
        )

        p_ratio.add_tools(hover)
        p_ratio.title.text_font_size = "18pt"
        p_ratio.xaxis.axis_label_text_font_size = "16pt"
        p_ratio.yaxis.axis_label_text_font_size = "16pt"
        p_ratio.xaxis.major_label_text_font_size = "12pt"
        p_ratio.yaxis.major_label_text_font_size = "12pt"

        # Depth
        p_depth = figure(
            title=f"Depth pre and post base-quality filtering {chrom}:{start}-{stop}",
            width=cyclomics_defaults.width,
        )
        p_depth.xaxis.axis_label = "Position"
        p_depth.yaxis.axis_label = "Depth"
        p_depth.xaxis.formatter.use_scientific = False

        p_depth.scatter(
            "POS", "DP", source=data_relevant, color="orange", legend_label="depth"
        )
        p_depth.scatter(
            "POS", "DPQ", source=data_relevant, legend_label="depth after Q filtering"
        )
        p_depth.add_tools(hover)
        p_depth.legend.location = "bottom_center"
        p_depth.title.text_font_size = "18pt"
        p_depth.xaxis.axis_label_text_font_size = "16pt"
        p_depth.yaxis.axis_label_text_font_size = "16pt"
        p_depth.xaxis.major_label_text_font_size = "12pt"
        p_depth.yaxis.major_label_text_font_size = "12pt"

        plots.append(p_vaf)
        plots.append(p_ratio)
        plots.append(p_depth)

    return column(*plots)


def main(vcf_file, plot_file, priority_limit: int):
    data = read_vcf(vcf_file)
    tab_name = "Variants"
    json_obj = {}
    json_obj[tab_name] = {}
    json_obj[tab_name]["name"] = tab_name
    json_obj[tab_name]["priority"] = TAB_PRIORITY

    if TAB_PRIORITY < priority_limit:
        if data.empty:
            json_obj[tab_name]["script"], json_obj[tab_name]["div"] = (
                "",
                "<h1>One of the pileups was not deep enough.</h1>",
            )
            with open(Path(plot_file).with_suffix(".json"), "w") as f:
                f.write(json.dumps(json_obj))
            return

        roi = get_roi_pileup_df(data)
        # print(data)
        plot = make_scatter_plots(data, roi)

        json_obj[tab_name]["script"], json_obj[tab_name]["div"] = components(plot)

    with open(Path(plot_file).with_suffix(".json"), "w") as f:
        f.write(json.dumps(json_obj))


if __name__ == "__main__":
    dev = False
    if not dev:
        import argparse

        parser = argparse.ArgumentParser(
            description="Create hist plot from a regex for fastq and fastq.gz files."
        )

        parser.add_argument("vcf_file")
        parser.add_argument("plot_file")
        parser.add_argument("priority_limit", type=int, default=89)
        args = parser.parse_args()

        main(args.vcf_file, args.plot_file, args.priority_limit)

    else:
        vcf_file = "ABZ922.noisy_merged.vcf"
        plot_file = "ABZ922.noisy_merged.html"
        priority_limit = 9999

        main(vcf_file, plot_file, priority_limit)
