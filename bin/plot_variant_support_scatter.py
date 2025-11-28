#!/usr/bin/env python

import io
import json
from pathlib import Path
from typing import Any

import pandas as pd
import plotly.graph_objects as go
import pysam
from plotly.subplots import make_subplots
from plotting_defaults import plotly_components

TAB_PRIORITY = 95


def relaxed_float(x: Any) -> float:
    """Return a float, with value error catch."""
    try:
        return float(x)
    except ValueError:
        return x


def read_vcf(path: Path) -> pd.DataFrame:
    """Read VCF file into a DataFrame with formatted genotype columns."""
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


def determine_variant_is_snp(vcf_entry: pd.Series) -> bool:
    """Return True if variant is a SNP."""
    return vcf_entry.REF.upper() in ["A", "C", "G", "T"] and vcf_entry.ALT.upper() in [
        "A",
        "C",
        "G",
        "T",
    ]


def get_relevant_YN_support(
    YN: str, query_position: int, reverse: bool, entry_split="|", subentry_split=","
) -> str:
    """Extract relevant YN subentry for a given position."""
    YN_struct = YN.split(entry_split)
    YN_support = YN_struct[2:]
    if reverse:
        YN_support = YN_support[::-1]

    YN_pos, YN_skip = 0, 0
    for base in YN_support:
        sub_struct = base.split(subentry_split)
        if sub_struct[0].startswith("D"):
            YN_skip += 1
            continue
        if YN_pos == query_position:
            return sub_struct
        YN_pos += 1
    return ""


def count_base_support(
    df: pd.DataFrame, letters=["A", "C", "D", "G", "T"]
) -> pd.DataFrame:
    """Count bases supporting each letter from YN_support."""

    def count_initials(lst, letter):
        return sum(1 for x in lst if x and x[0].upper() == letter)

    for letter in letters:
        df[f"{letter}_count"] = df["YN_support"].apply(
            lambda lst: count_initials(lst, letter)
        )
    return df


def _add_reference_information_snp(
    df: pd.DataFrame, vcf_entry: pd.Series
) -> pd.DataFrame:
    """Compute reference/alt allele fractions."""
    ref_col = f"{vcf_entry.REF.upper()}_count"
    alt_col = f"{vcf_entry.ALT.upper()}_count"

    df["total_count"] = df[["A_count", "C_count", "G_count", "T_count"]].sum(axis=1)
    df["ref_frac_all"] = df[ref_col] / df["total_count"]
    df["alt_frac_all"] = df[alt_col] / df["total_count"]

    df["ref_alt_total"] = df[ref_col] + df[alt_col]
    df["ref_frac_ref_alt"] = df[ref_col] / df["ref_alt_total"]
    df["alt_frac_ref_alt"] = df[alt_col] / df["ref_alt_total"]

    df.fillna(0, inplace=True)
    return df


def get_read_support_snp(vcf_entry: pd.Series, bam: pysam.AlignmentFile):
    """Collect read support for SNP."""
    contig, start, stop = vcf_entry.CHROM, vcf_entry.POS - 1, vcf_entry.POS
    read_support = []
    for pileupcolumn in bam.pileup(contig, start, stop, stepper="all", truncate=True):
        for pileupread in pileupcolumn.pileups:
            if pileupread.query_position is None:
                continue
            YN_support = get_relevant_YN_support(
                pileupread.alignment.get_tag("YN"),
                pileupread.query_position,
                reverse=not pileupread.alignment.is_forward,
            )
            if not YN_support:
                continue
            read_support.append(
                [
                    pileupread.alignment.query_name,
                    pileupread.alignment.is_forward,
                    pileupread.alignment.get_tag("YM"),
                    YN_support[0][0],
                    YN_support[1:],
                ]
            )
    return read_support


def plot_variant(df: pd.DataFrame, vcf_entry: pd.Series):
    """Plot SNP variant with scatter + histograms using Plotly."""
    igv_colors = {
        "A": "#00C000",
        "C": "#0000C0",
        "G": "#FFA500",
        "T": "#FF0000",
        "D": "#808080",
    }

    fig = make_subplots(
        rows=2,
        cols=2,
        column_widths=[0.8, 0.2],
        row_heights=[0.2, 0.8],
        specs=[[{"type": "xy"}, None], [{"type": "xy"}, {"type": "xy"}]],
        horizontal_spacing=0.02,
        vertical_spacing=0.02,
    )

    # Scatter
    colors = [igv_colors.get(base, "#888888") for base in df["YN_base"]]
    fig.add_trace(
        go.Scatter(
            x=df["total_count"],
            y=df["alt_frac_all"],
            mode="markers",
            marker=dict(size=6, color=colors, opacity=0.6),
            text=df["YN_base"],
            name="Scatter",
        ),
        row=2,
        col=1,
    )

    # X histogram
    fig.add_trace(
        go.Histogram(
            x=df["total_count"],
            nbinsx=25,
            marker_color="gray",
            name="X hist",
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    # Y histogram
    fig.add_trace(
        go.Histogram(
            y=df["alt_frac_all"],
            nbinsy=53,
            marker_color="gray",
            name="Y hist",
            showlegend=False,
        ),
        row=2,
        col=2,
    )

    # Layout
    fig.update_layout(
        title=f"Variant Support Scatter Plot for {vcf_entry.CHROM}:{vcf_entry.POS} "
        f"({vcf_entry.REF.upper()} → {vcf_entry.ALT.upper()})",
        width=900,
        height=700,
        xaxis_title="Total count (repeats)",
        yaxis_title="Alt allele fraction",
        bargap=0.05,
        hovermode="closest",
        showlegend=False,
    )
    fig.update_yaxes(range=[-0.02, 1.02], row=2, col=1)
    fig.update_xaxes(range=[0, 25], row=2, col=1)

    return fig


def main(vcf_path: Path, bam_path: Path, output_path: Path, priority_limit: int = 0):
    """Main: generate variant scatter plot JSON output."""
    tab_name = "Variant support scatter plots"
    json_obj = {tab_name: {"name": tab_name, "priority": TAB_PRIORITY}}
    add_info = {}

    if TAB_PRIORITY > priority_limit:
        with open(Path(output_path).with_suffix(".json"), "w") as f:
            f.write(json.dumps(json_obj))
        return

    vcf_path, bam_path, output_path = Path(vcf_path), Path(bam_path), Path(output_path)
    bam = pysam.AlignmentFile(bam_path, "rb")
    vcf = read_vcf(vcf_path)

    figs = []
    for _, vcf_entry in vcf.iterrows():
        if determine_variant_is_snp(vcf_entry):
            read_support = get_read_support_snp(vcf_entry, bam)
            df = pd.DataFrame(
                read_support,
                columns=["read_name", "fwd", "YM", "YN_base", "YN_support"],
            )
            df = count_base_support(df)
            df = _add_reference_information_snp(df, vcf_entry)
            figs.append(plot_variant(df, vcf_entry))

    if not figs:
        json_obj[tab_name]["script"], json_obj[tab_name]["div"] = (
            "",
            "<h1>No SNP variants were plottable (indels skipped).</h1>",
        )
    else:
        # Combine plots vertically
        html_blocks = [plotly_components(fig) for fig in figs]
        scripts = "".join(script for script, _ in html_blocks)
        divs = "".join(div for _, div in html_blocks)
        json_obj[tab_name]["script"], json_obj[tab_name]["div"] = scripts, divs

    json_obj["additional_info"] = add_info
    with open(output_path, "w") as f:
        f.write(json.dumps(json_obj, indent=2))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plot variant support scatter plot.")
    parser.add_argument("vcf_path", type=str, help="Path to the VCF file.")
    parser.add_argument("bam_path", type=str, help="Path to the BAM file.")
    parser.add_argument("output_path", type=str, help="Path to save the output JSON.")
    parser.add_argument("priority_limit", type=int, default=0)
    args = parser.parse_args()

    main(args.vcf_path, args.bam_path, args.output_path, args.priority_limit)
