#!/usr/bin/env python

import io
from typing import List, Tuple
import pandas as pd
import numpy as np

from bokeh.io import save, output_file
from bokeh.plotting import figure, show
from bokeh.layouts import row, column
from bokeh.models import HoverTool

chromosomal_region = Tuple[str, int, int]


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
    df = pd.read_csv(
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
    ).rename(columns={"#CHROM": "CHROM"})

    formats = df.FORMAT[0].split(":")
    for i, fmt in enumerate(formats):
        df[fmt] = df.Sample1.apply(
            lambda x: float((x.split(":")[i] if (x.split(":")[i]) else 0))
        )
    del df["FORMAT"]
    del df["Sample1"]
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


        # Pos scatter vaf
        p_vaf = figure(title=f"Positional VAF {chrom}:{start}-{stop}")
        p_vaf.xaxis.axis_label = "Position"
        p_vaf.yaxis.axis_label = "VAF ratio"
        p_vaf.scatter("POS", "VAF", source=data_relevant)
        p_vaf.add_tools(hover)
        # show(p_vaf)

        # fwd ratio vs reverse ratio
        p_ratio = figure(title=f"Support ratio's per direction {chrom}:{start}-{stop}")
        p_ratio.xaxis.axis_label = "forward  ratio"
        p_ratio.yaxis.axis_label = "reverse ratio"
        p_ratio.scatter("FWDR", "REVR", source=data_relevant)
        p_ratio.add_tools(hover)

        # Depth
        p_depth = figure(title=f"Depth pre and post base-quality filtering {chrom}:{start}-{stop}")
        p_depth.xaxis.axis_label = "Position"
        p_depth.yaxis.axis_label = "Depth"
        p_depth.scatter("POS", "DP", source=data_relevant, color="orange", legend_label="depth")
        p_depth.scatter("POS", "DPQ", source=data_relevant, legend_label="depth after Q filtering")
        p_depth.add_tools(hover)
        p_depth.legend.location = "bottom_center"
        plots.append(row(p_vaf, p_ratio, p_depth))

    return column(*plots)


def main(vcf_file, plot_file):
    data = read_vcf(args.vcf_file)
    roi = get_roi_pileup_df(data)
    # print(data)
    plot = make_scatter_plots(data, roi)
    output_file(args.plot_file, title="variant plots")
    save(plot)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Create hist plot from a regex for fastq and fastq.gz files."
    )

    parser.add_argument("vcf_file")
    parser.add_argument("plot_file")
    args = parser.parse_args()
    # vcf = '/home/dami/projects/variantcalling/depth/datasets/000010_v4_bb41/2000/potential_vars_1000.vcf'
    # true_vcf = '/home/dami/projects/variantcalling/depth/datasets_cyclomics/000010_v4_bb41/Native/real_variants.vcf'

    main(args.vcf_file,args.plot_file)