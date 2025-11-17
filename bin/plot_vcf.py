#!/usr/bin/env python

import io
import json
from pathlib import Path
from typing import Any, List, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotting_defaults import cyclomics_defaults, plotly_components

chromosomal_region = Tuple[str, int, int]
TAB_PRIORITY = 3


def relaxed_float(x: Any) -> float:
    """Return a float, with value error catch."""
    try:
        return float(x)
    except ValueError:
        return x


def get_roi_pileup_df(
    df: pd.DataFrame, distance: int = 100
) -> List[chromosomal_region]:
    """Find chromosomal stretches where variants are < distance apart."""

    def get_continous_strech(array: np.ndarray, distance: int):
        m = np.concatenate(([True], array[1:] > array[:-1] + distance, [True]))
        idx = np.flatnonzero(m)
        l = array.tolist()
        return [l[i:j] for i, j in zip(idx[:-1], idx[1:])]

    assemblies = df.CHROM.unique()
    roi = []
    for chrom in assemblies:
        df_tmp = df[df.CHROM == chrom]
        stretches = get_continous_strech(df_tmp.POS.to_numpy(), distance)
        for stretch in stretches:
            roi.append((chrom, stretch[0], stretch[-1]))

    return roi


def read_vcf(path: Path) -> pd.DataFrame:
    """Read a VCF file into a DataFrame."""
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

    df.drop(columns=["FORMAT", "SAMPLE1"], inplace=True)
    return df


def update_vcf_with_truth(true_vcf: Path, df: pd.DataFrame) -> pd.DataFrame:
    """Add truth set labels to VCF data."""
    true_df = read_vcf(true_vcf)[["CHROM", "POS"]]
    true_df["Variant"] = 1
    df = pd.merge(df, true_df, how="left", on=["CHROM", "POS"])
    df["Variant"].fillna(0, inplace=True)
    return df


def make_scatter_plots(data: pd.DataFrame, roi: List[chromosomal_region]):
    """Make interactive scatter plots for VAF, direction ratios, and depth."""
    plots = []

    for chrom, start, stop in roi:
        data_relevant = data[
            (data.CHROM == chrom) & (data.POS >= start) & (data.POS <= stop)
        ].copy()

        data_snp = data_relevant[
            data_relevant.REF.str.len() == data_relevant.ALT.str.len()
        ]
        data_indel = data_relevant[
            data_relevant.REF.str.len() != data_relevant.ALT.str.len()
        ]

        # Create 3 stacked subplots for each region
        fig = make_subplots(
            rows=3,
            cols=1,
            subplot_titles=(
                f"Positional VAF {chrom}:{start}-{stop}",
                f"Support ratios per direction {chrom}:{start}-{stop}",
                f"Depth pre and post base-quality filtering {chrom}:{start}-{stop}",
            ),
            vertical_spacing=0.1,
        )

        # --- 1. Positional VAF ---
        fig.add_trace(
            go.Scatter(
                x=data_snp["POS"],
                y=data_snp["VAF"],
                mode="markers",
                name="SNP",
                marker=dict(color="blue", size=6),
                hovertext=[
                    f"CHROM: {r.CHROM}<br>POS: {r.POS}<br>REF→ALT: {r.REF}->{r.ALT}<br>VAF: {r.VAF}"
                    for _, r in data_snp.iterrows()
                ],
                hoverinfo="text",
            ),
            row=1,
            col=1,
        )
        fig.add_trace(
            go.Scatter(
                x=data_indel["POS"],
                y=data_indel["VAF"],
                mode="markers",
                name="InDel",
                marker=dict(color="red", size=6),
                hovertext=[
                    f"CHROM: {r.CHROM}<br>POS: {r.POS}<br>REF→ALT: {r.REF}->{r.ALT}<br>VAF: {r.VAF}"
                    for _, r in data_indel.iterrows()
                ],
                hoverinfo="text",
            ),
            row=1,
            col=1,
        )

        # --- 2. Forward vs reverse ratio ---
        fig.add_trace(
            go.Scatter(
                x=data_snp["FWDR"],
                y=data_snp["REVR"],
                mode="markers",
                name="SNP",
                marker=dict(color="blue", size=6),
                hovertext=[
                    f"FWDR: {r.FWDR}<br>REVR: {r.REVR}<br>VAF: {r.VAF}<br>DP: {r.DP}"
                    for _, r in data_snp.iterrows()
                ],
                hoverinfo="text",
            ),
            row=2,
            col=1,
        )
        fig.add_trace(
            go.Scatter(
                x=data_indel["FWDR"],
                y=data_indel["REVR"],
                mode="markers",
                name="InDel",
                marker=dict(color="red", size=6),
                hovertext=[
                    f"FWDR: {r.FWDR}<br>REVR: {r.REVR}<br>VAF: {r.VAF}<br>DP: {r.DP}"
                    for _, r in data_indel.iterrows()
                ],
                hoverinfo="text",
            ),
            row=2,
            col=1,
        )

        # --- 3. Depth pre/post filtering ---
        fig.add_trace(
            go.Scatter(
                x=data_relevant["POS"],
                y=data_relevant["DP"],
                mode="markers",
                name="Depth",
                marker=dict(color="orange", size=5),
                hovertext=[
                    f"POS: {r.POS}<br>DP: {r.DP}<br>DPQ: {r.DPQ}"
                    for _, r in data_relevant.iterrows()
                ],
                hoverinfo="text",
            ),
            row=3,
            col=1,
        )
        fig.add_trace(
            go.Scatter(
                x=data_relevant["POS"],
                y=data_relevant["DPQ"],
                mode="markers",
                name="Depth after Q filtering",
                marker=dict(color="green", size=5),
                hovertext=[
                    f"POS: {r.POS}<br>DP: {r.DP}<br>DPQ: {r.DPQ}"
                    for _, r in data_relevant.iterrows()
                ],
                hoverinfo="text",
            ),
            row=3,
            col=1,
        )

        # --- Layout ---
        fig.update_layout(
            width=cyclomics_defaults.width,
            height=900,
            showlegend=True,
            legend=dict(orientation="h", y=-0.1),
            hovermode="closest",
            template="plotly_white",
        )

        fig.update_xaxes(title_text="Position", row=1, col=1)
        fig.update_yaxes(title_text="VAF", row=1, col=1)
        fig.update_xaxes(title_text="Forward ratio", row=2, col=1)
        fig.update_yaxes(title_text="Reverse ratio", row=2, col=1)
        fig.update_xaxes(title_text="Position", row=3, col=1)
        fig.update_yaxes(title_text="Depth", row=3, col=1)

        plots.append(fig)

    return plots


def main(vcf_file: str, plot_file: str, priority_limit: int):
    """Generate interactive variant plots and export JSON for embedding."""
    data = read_vcf(vcf_file)
    tab_name = "Variants"
    json_obj = {tab_name: {"name": tab_name, "priority": TAB_PRIORITY}}

    if TAB_PRIORITY < priority_limit:
        if data.empty:
            json_obj[tab_name]["script"], json_obj[tab_name]["div"] = (
                [],
                ["<h1>No variants were found.</h1>"],
            )
        else:
            roi = get_roi_pileup_df(data)
            figs = make_scatter_plots(data, roi)
            html_blocks = [plotly_components(fig) for fig in figs]
            scripts = "".join(script for script, _ in html_blocks)
            divs = "".join(div for _, div in html_blocks)
            json_obj[tab_name]["script"], json_obj[tab_name]["div"] = scripts, divs

    with open(Path(plot_file).with_suffix(".json"), "w") as f:
        f.write(json.dumps(json_obj, indent=2))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Create variant scatter plots.")
    parser.add_argument("vcf_file")
    parser.add_argument("plot_file")
    parser.add_argument("priority_limit", type=int, default=89)
    args = parser.parse_args()

    main(args.vcf_file, args.plot_file, args.priority_limit)
