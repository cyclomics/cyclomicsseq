#!/usr/bin/env python

import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotting_defaults import plotly_components

chromosomal_region = Tuple[str, int, int]
SIDE_DIST_PLOT_SIZE = 100
TAB_PRIORITY = 91


def get_roi_pileup_df(
    df: pd.DataFrame, distance: int = 100
) -> List[chromosomal_region]:
    print("get_roi_pileup_df)")
    """
    Given a dataframe with CHROM and POS columns, find all streches in chromosomes that are less than `distance` appart.
    reports 

    """

    def get_continous_strech(array: np.array, distance: int):
        print("get_continous_strech)")
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


def acc_to_Q(acc: float, max_Q=50) -> int:
    """
    Convert an accuracy (eg 0.99) to a Q score
    """
    # perfect acc is not Q score convertable
    if acc == 1:
        return max_Q
    else:
        return -10 * math.log10(1 - acc)


def _parse_pileup(row):
    print("_parse_pileup)")
    """
    . (dot) means a base that matched the reference on the forward strand
    , (comma) means a base that matched the reference on the reverse strand
    </> (less-/greater-than sign) denotes a reference skip. This occurs, for example, if a base in the reference genome is intronic and a read maps to two flanking exons. If quality scores are given in a sixth column, they refer to the quality of the read and not the specific base.
    AGTCN (upper case) denotes a base that did not match the reference on the forward strand
    agtcn (lower case) denotes a base that did not match the reference on the reverse strand
    A sequence matching the regular expression \+[0-9]+[ACGTNacgtn]+ denotes an insertion of one or more bases starting from the next position. For example, +2AG means insertion of AG in the forward strand
    A sequence matching the regular expression \-[0-9]+[ACGTNacgtn]+ denotes a deletion of one or more bases starting from the next position. For example, -2ct means deletion of CT in the reverse strand
    ^ (caret) marks the start of a read segment and the ASCII of the character following `^' minus 33 gives the mapping quality
    $ (dollar) marks the end of a read segment
    * (asterisk) is a placeholder for a deleted base in a multiple basepair deletion that was mentioned in a previous line by the -[0-9]+[ACGTNacgtn]+ notation
    """
    pileup = row.PILEUP[:]
    starts = []
    while True:
        try:
            i = pileup.index("^")
            starts.append(pileup[i : i + 3])
            pileup = pileup[:i] + pileup[i + 3 :]
        except ValueError:
            break

    ends = []
    while True:
        try:
            i = pileup.index("$")
            ends.append(pileup[i - 1 : i + 1])
            pileup = pileup[: i - 1] + pileup[i + 1 :]
        except ValueError:
            break

    skips = []
    while True:
        try:
            i = pileup.index(">")
            skips.append(pileup[i])
            pileup = pileup[:i] + pileup[i + 1 :]
        except ValueError:
            break
    while True:
        try:
            i = pileup.index("<")
            skips.append(pileup[i])
            pileup = pileup[:i] + pileup[i + 1 :]
        except ValueError:
            break

    forw_ins = []
    p = re.compile("\+[0-9]+[ACGTN]+")
    matches = list(p.finditer(pileup))
    spans = reversed([m.span() for m in matches])
    forw_ins += [m.group() for m in matches]
    for i, j in spans:
        pileup = pileup[:i] + pileup[j:]

    rev_ins = []
    p = re.compile("\+[0-9]+[acgtn]+")
    matches = list(p.finditer(pileup))
    spans = reversed([m.span() for m in matches])
    rev_ins += [m.group() for m in matches]
    for i, j in spans:
        pileup = pileup[:i] + pileup[j:]

    forw_del = []
    p = re.compile("\-[0-9]+[ACGTN]+")
    matches = list(p.finditer(pileup))
    spans = reversed([m.span() for m in matches])
    forw_del += [m.group() for m in matches]
    for i, j in spans:
        pileup = pileup[:i] + pileup[j:]

    while True:
        try:
            i = pileup.index("*")
            forw_del.append(pileup[i])
            pileup = pileup[:i] + pileup[i + 1 :]
        except ValueError:
            break

    rev_del = []
    p = re.compile("\-[0-9]+[acgtn]+")
    matches = list(p.finditer(pileup))
    spans = reversed([m.span() for m in matches])
    rev_del += [m.group() for m in matches]
    for i, j in spans:
        pileup = pileup[:i] + pileup[j:]

    while True:
        try:
            i = pileup.index("#")  # del in reverse strand
            rev_del.append(pileup[i])
            pileup = pileup[:i] + pileup[i + 1 :]
        except ValueError:
            break

    # Put back start and end bases
    for seq in starts:
        pileup += seq[-1]
    for seq in ends:
        pileup += seq[0]

    # Count matches/mismatches
    forw_matches = pileup.count(".")
    rev_matches = pileup.count(",")
    pileup = pileup.replace(".", "").replace(",", "")
    snps = Counter(pileup)

    calls = Counter(
        {
            "A": 0,
            "T": 0,
            "C": 0,
            "G": 0,
            "N": 0,
            "a": 0,
            "t": 0,
            "c": 0,
            "g": 0,
            "n": 0,
            "fINS": 0,
            "fDEL": 0,
            "rINS": 0,
            "rDEL": 0,
            "SKIP": 0,
            "FMB": 0,
            "LMB": 0,
        }
    )

    calls.update(snps)
    calls.update({"FMB": len(starts)})  # first mapping base
    calls.update({"LMB": len(ends)})  # last mapping base
    calls.update({"fDEL": len(forw_del)})
    calls.update({"fINS": len(forw_ins)})
    calls.update({"rDEL": len(forw_del)})
    calls.update({"rINS": len(forw_ins)})
    calls.update({row.REF: forw_matches})
    calls.update({row.REF.lower(): rev_matches})
    return pd.Series(calls)


def pileup_to_df(pileup_path: str) -> pd.DataFrame:
    print("pileup_to_df)")

    with open(pileup_path, "r") as pu:
        lines = pu.readlines()
        lines = [x.split("\t") for x in lines]
        df1 = pd.DataFrame(
            data=lines,
            columns=["CHROM", "POS", "REF", "COV", "PILEUP", "BASEQ", "MAPQ", "READ"],
        )
        df1["POS"] = df1.POS.astype(int)
        df1["COV"] = df1.COV.astype(int)
        df1["CONCAT"] = [len(set(v.split(","))) for v in df1.READ.values]  # concatemers
        df1 = df1.join(df1.apply(_parse_pileup, axis=1), lsuffix="", rsuffix="_right")
        df1["ACCURACY"] = df1.apply(
            compute_accuracy, insertions=True, deletions=True, axis=1
        )
        df1["Q"] = df1.ACCURACY.apply(acc_to_Q)

        return df1


def perbase_table_to_df(table_path):
    """
    https://github.com/sstadick/perbase

    """
    try:
        df = pd.read_csv(table_path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()
    df["REF_BASE_UPPER"] = df.REF_BASE.str.upper()
    idx, cols = pd.factorize(df["REF_BASE_UPPER"])
    df["REF_COUNT"] = df.reindex(cols, axis=1).to_numpy()[np.arange(len(df)), idx]
    df["ACCURACY"] = df.REF_COUNT / df.DEPTH
    df["Q"] = df.ACCURACY.apply(acc_to_Q)
    df["CHROM"] = df["REF"]

    return df


def compute_accuracy(
    row: pd.Series, insertions=True, deletions=True, by_strand=False
) -> float:
    print("compute_accuracy)")
    """Compute base call accuracy as reference-calls/coverage"""
    # designed to be used with df.apply(compute_accuracy, axis=1)
    # insertions/deletions determine wether or not ins and dels will be counted in the COV

    # expected row keys:
    #   ['A','T','C','G','N','a','t','c','g','n','fINS','fDEL','rINS','rDEL']

    ref = row.REF
    if ref == "N":
        return np.nan

    if not by_strand:
        ref_calls = row[ref.upper()] + row[ref.lower()]
        cov = row[
            [
                "A",
                "T",
                "C",
                "G",
                "N",
                "a",
                "t",
                "c",
                "g",
                "n",
                "fINS",
                "fDEL",
                "rINS",
                "rDEL",
            ]
        ].sum()
        if not insertions:
            cov -= row["fINS"] + row["rINS"]
        if not deletions:
            cov -= row["fDEL"] + row["rDEL"]

        if cov > 0:
            return ref_calls / cov
        else:
            return 0

    else:
        # Forward
        fref_calls = row[ref.upper()]
        fcov = row[["A", "T", "C", "G", "N", "fINS", "fDEL"]].sum()
        if not insertions:
            fcov -= row["fINS"]
        if not deletions:
            fcov -= row["fDEL"]
        # Reverse
        rref_calls = row[ref.lower()]
        rcov = row[["a", "t", "c", "g", "n", "rINS", "rDEL"]].sum()
        if not insertions:
            rcov -= row["rINS"]
        if not deletions:
            rcov -= row["rDEL"]

        return pd.Series([fref_calls / fcov, rref_calls / rcov])


def plot_compare_accuracy(
    roi: List[chromosomal_region],
    df_list: List[pd.DataFrame],
    my_title="untitled",
) -> list[go.Figure]:
    """
    Plot the Q score against position for each region of interest.
    Returns a list of Plotly figures.
    """
    plots = []
    for chrom, start, stop in roi:
        df1 = df_list[0]
        df2 = df_list[1]

        df1_region = df1.query("CHROM == @chrom and @start <= POS <= @stop")
        df2_region = df2.query("CHROM == @chrom and @start <= POS <= @stop")

        if df1_region.empty or df2_region.empty:
            continue

        fig = go.Figure()

        # Raw Q line
        fig.add_trace(
            go.Scatter(
                x=df1_region["POS"],
                y=df1_region["Q"],
                mode="lines",
                name="Raw",
                line=dict(color="steelblue", width=2),
                hovertemplate="CHROM=%{customdata[0]}<br>POS=%{x}<br>Q_raw=%{y}<extra></extra>",
                customdata=np.stack([df1_region["CHROM"]], axis=-1),
            )
        )

        # Consensus Q line
        fig.add_trace(
            go.Scatter(
                x=df2_region["POS"],
                y=df2_region["Q"],
                mode="lines",
                name="Consensus",
                line=dict(color="orange", width=2),
                hovertemplate="CHROM=%{customdata[0]}<br>POS=%{x}<br>Q_cons=%{y}<extra></extra>",
                customdata=np.stack([df2_region["CHROM"]], axis=-1),
            )
        )

        fig.update_layout(
            title=f"{my_title} {chrom}:{start}-{stop}",
            xaxis_title="Position",
            yaxis_title="Q Score",
            legend=dict(x=0.01, y=0.99, bgcolor="rgba(255,255,255,0.6)"),
            width=700,
            height=400,
            margin=dict(l=40, r=20, t=60, b=40),
        )
        plots.append(fig)
    return plots


def make_qscore_scatter(df1, df2, csv_merge=None):
    """
    Scatter of Q_raw vs Q_cons with marginal histograms.
    """
    pre = df1[["CHROM", "POS", "Q"]]
    post = df2[["CHROM", "POS", "Q"]]

    df_merge = pd.merge(pre, post, how="left", on=["CHROM", "POS"])
    df_merge.rename(columns={"Q_x": "Q_raw", "Q_y": "Q_cons"}, inplace=True)

    if csv_merge:
        df_merge.to_csv(csv_merge, index=False)

    fig = make_subplots(
        rows=2,
        cols=2,
        column_widths=[0.8, 0.2],
        row_heights=[0.2, 0.8],
        specs=[
            [{}, {"rowspan": 1, "colspan": 1}],
            [{"colspan": 1, "rowspan": 1}, None],
        ],
        horizontal_spacing=0.02,
        vertical_spacing=0.02,
    )

    # Scatter
    fig.add_trace(
        go.Scatter(
            x=df_merge["Q_raw"],
            y=df_merge["Q_cons"],
            mode="markers",
            marker=dict(color="orange", size=6, opacity=0.5),
            name="Q scores",
            hovertemplate="CHROM=%{customdata[0]}<br>POS=%{customdata[1]}<br>Raw=%{x}<br>Cons=%{y}<extra></extra>",
            customdata=np.stack([df_merge["CHROM"], df_merge["POS"]], axis=-1),
        ),
        row=2,
        col=1,
    )

    # Diagonal line y=x
    fig.add_trace(
        go.Scatter(
            x=[0, 50],
            y=[0, 50],
            mode="lines",
            line=dict(color="grey", width=1, dash="dash"),
            showlegend=False,
        ),
        row=2,
        col=1,
    )

    # Top histogram (Q_raw)
    fig.add_trace(
        go.Histogram(
            x=df_merge["Q_raw"],
            nbinsx=50,
            marker_color="steelblue",
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    # Right histogram (Q_cons)
    fig.add_trace(
        go.Histogram(
            y=df_merge["Q_cons"],
            nbinsy=50,
            marker_color="orange",
            showlegend=False,
        ),
        row=2,
        col=2,
    )

    fig.update_layout(
        title="Q score scatter, variant unaware",
        width=700,
        height=700,
        xaxis2={"showticklabels": False},
        yaxis1={"showticklabels": False},
        xaxis_title="ONT Q score",
        yaxis_title="Cyclomics Q score",
        bargap=0.05,
    )

    return fig


def main(perbase_path1, perbase_path2, output_plot_file, priority_limit: int):
    tab_name = "Consensus quality"
    json_obj = {}
    json_obj["additional_info"] = {}  # NOTE: Do we need this?
    json_obj[tab_name] = {}
    json_obj[tab_name]["name"] = tab_name
    json_obj[tab_name]["priority"] = TAB_PRIORITY

    if TAB_PRIORITY < priority_limit:
        df1 = perbase_table_to_df(perbase_path1)
        df2 = perbase_table_to_df(perbase_path2)

        if df1.empty or df2.empty:
            json_obj[tab_name]["script"], json_obj[tab_name]["div"] = (
                [],
                ["<h1>No loci were covered deeply enough to compare accuracy.</h1>"],
            )

        else:
            roi = get_roi_pileup_df(df1)

            q_score_plot = make_qscore_scatter(df1, df2)
            positional_accuracy_figs = plot_compare_accuracy(
                roi, [df1, df2], "Variant unaware Q"
            )

            scripts, divs = plotly_components([q_score_plot, *positional_accuracy_figs])
            json_obj[tab_name]["script"] = scripts
            json_obj[tab_name]["div"] = divs

    with open(Path(output_plot_file).with_suffix(".json"), "w", encoding="utf-8") as f:
        json.dump(json_obj, f, ensure_ascii=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Create hist plot from a regex for fastq and fastq.gz files."
    )

    parser.add_argument("pileup_split", type=Path)
    parser.add_argument("pileup_consensus", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("priority_limit", type=int, default=89)

    args = parser.parse_args()
    main(args.pileup_split, args.pileup_consensus, args.output, args.priority_limit)
