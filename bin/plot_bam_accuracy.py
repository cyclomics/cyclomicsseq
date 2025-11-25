#!/usr/bin/env python

import argparse
import json
import math
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


def acc_to_Q(acc: float, max_Q=50) -> int:
    """
    Convert an accuracy (eg 0.99) to a Q score
    """
    # perfect acc is not Q score convertable
    if acc == 1:
        return max_Q
    else:
        return -10 * math.log10(1 - acc)


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
                hovertemplate="CHROM: %{customdata[0]}<br>POS: %{x}<br>Raw Q: %{y:.0f}<extra></extra>",
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
                hovertemplate="CHROM: %{customdata[0]}<br>POS: %{x}<br>Consensus Q: %{y:.0f}<extra></extra>",
                customdata=np.stack([df2_region["CHROM"]], axis=-1),
            )
        )

        fig.update_layout(
            title=f"Variant unaware positional Q scores: {chrom}:{start}-{stop}",
            template="plotly_white",
            xaxis_title=f"Position in {chrom}",
            yaxis_title="Q score",
            legend=dict(x=0.01, y=0.99, bgcolor="rgba(255,255,255,0.6)"),
            width=700,
            height=400,
            margin=dict(l=40, r=20, t=60, b=40),
        )
        fig.update_xaxes(tickformat=",d", tickangle=-30)

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
        column_widths=[0.75, 0.15],
        row_heights=[0.15, 0.75],
        specs=[
            [{}, {"rowspan": 1, "colspan": 1}],
            [{"colspan": 1, "rowspan": 1}, {}],
        ],
        horizontal_spacing=0.05,
        vertical_spacing=0.05,
    )

    # Scatter
    fig.add_trace(
        go.Scatter(
            x=df_merge["Q_raw"],
            y=df_merge["Q_cons"],
            mode="markers",
            marker=dict(color="orange", size=6, opacity=0.5),
            showlegend=False,
            hoverinfo="skip",
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
            hoverinfo="skip",
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
            histnorm="percent",
            hovertemplate="Raw Q: %{x:.0f}<br>%{y:.1f}<extra></extra>",
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
            histnorm="percent",
            hovertemplate="Consensus Q: %{y:.0f}<br>%{x:.1f}<extra></extra>",
        ),
        row=2,
        col=2,
    )

    fig.update_layout(
        title="Variant unaware Q score distribution",
        template="plotly_white",
        width=700,
        height=700,
        xaxis2={"showticklabels": False},
        yaxis1={"showticklabels": True},
        bargap=0.05,
    )

    # Update axes on main plot
    fig.update_xaxes(title_text="Raw Q score", row=2, col=1)
    fig.update_yaxes(title_text="Consensus Q score", row=2, col=1)

    # Update axes on top histogram
    fig.update_xaxes(showticklabels=False, row=1, col=1)
    fig.update_yaxes(title_text="Percent", row=1, col=1)

    # Update axes right histogram
    fig.update_xaxes(title_text="Percent", row=2, col=2)
    fig.update_yaxes(showticklabels=False, row=2, col=2)

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
            positional_accuracy_figs = plot_compare_accuracy(roi, [df1, df2])

            scripts, divs = plotly_components([q_score_plot, *positional_accuracy_figs])
            json_obj[tab_name]["script"] = scripts
            json_obj[tab_name]["div"] = divs

    with open(Path(output_plot_file).with_suffix(".json"), "w", encoding="utf-8") as f:
        json.dump(json_obj, f, ensure_ascii=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create hist plot from a regex for fastq and fastq.gz files."
    )

    parser.add_argument("pileup_split", type=Path)
    parser.add_argument("pileup_consensus", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("priority_limit", type=int, default=89)

    args = parser.parse_args()
    main(args.pileup_split, args.pileup_consensus, args.output, args.priority_limit)
