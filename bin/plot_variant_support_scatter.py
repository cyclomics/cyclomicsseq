import io
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pysam
from bokeh.embed import components
from bokeh.layouts import column, gridplot
from bokeh.models import CategoricalColorMapper, ColumnDataSource
from bokeh.plotting import figure
from bokeh.transform import jitter

TAB_PRIORITY = 95


def relaxed_float(x: Any) -> float:
    """Return a float, with value error catch"""
    try:
        my_float = float(x)
    except ValueError:
        my_float = x
    return my_float


def read_vcf(path: Path) -> pd.DataFrame:
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


def get_relevant_YN_support(
    YN: str, query_position: int, reverse: bool, entry_split="|", subentry_split=","
) -> str:
    """
    Get the relevant YN support for a given query position.
    """
    YN_struct = YN.split(entry_split)
    YN_len = YN_struct[0]
    YN_bases = YN_struct[1]

    YN_support = YN_struct[2:]
    YN_pos = 0
    YN_skip = 0
    relevant_YN_struct = None

    if reverse:
        YN_support = YN_support[::-1]

    for i, base in enumerate(YN_support):
        sub_struct = base.split(subentry_split)

        if sub_struct[0].startswith("D"):
            YN_skip += 1
            continue

        elif YN_pos == query_position:

            relevant_YN_struct = sub_struct
            break

        YN_pos += 1

    if not relevant_YN_struct:
        relevant_YN_struct = ""
    return relevant_YN_struct


def count_base_support(
    df: pd.DataFrame, letters=["A", "C", "D", "G", "T"]
) -> pd.DataFrame:
    # inline func to apply to the YN_support column
    def count_initials(lst, letter):
        return sum(1 for x in lst if x and x[0].upper() == letter)

    # Apply the function per row
    for letter in letters:
        df[f"{letter}_count"] = df["YN_support"].apply(
            lambda lst: count_initials(lst, letter)
        )

    return df


def _add_reference_information_snp(df: pd.DataFrame, vcf_entry: str) -> pd.DataFrame:
    ref_base = vcf_entry.REF.upper()
    alt_base = vcf_entry.ALT.upper()

    # Make sure the base names match column format
    ref_col = f"{ref_base}_count"
    alt_col = f"{alt_base}_count"

    # Add the total fractions
    df["total_count"] = df[["A_count", "C_count", "G_count", "T_count"]].sum(axis=1)
    # Fractional reference support over all bases
    df["ref_frac_all"] = df[ref_col] / df["total_count"]
    # Fractional alt support over all bases
    df["alt_frac_all"] = df[alt_col] / df["total_count"]

    # Total over just ref + alt
    df["ref_alt_total"] = df[ref_col] + df[alt_col]
    # Fractional reference support over ref+alt
    df["ref_frac_ref_alt"] = df[ref_col] / df["ref_alt_total"]
    # Fractional alt support over ref+alt
    df["alt_frac_ref_alt"] = df[alt_col] / df["ref_alt_total"]

    # fix the NaN values due to division by zero
    df.fillna(0, inplace=True)
    return df


def _add_reference_information_insert(df: pd.DataFrame, vcf_entry: str) -> pd.DataFrame:
    return df


def _add_reference_information_deletion(
    df: pd.DataFrame, vcf_entry: str
) -> pd.DataFrame:
    ref_base = vcf_entry.REF.upper()
    alt_base = vcf_entry.ALT.upper()

    deleted_bases = ref_base[1:]

    return df


def determine_variant_is_snp(vcf_entry: pd.Series) -> bool:
    """
    Simply check if the variant is a Single nucleotide polymorphism.
    e.g. N-> N and not NN -> N or N-> NN
    """
    return vcf_entry.REF.upper() in ["A", "C", "G", "T"] and vcf_entry.ALT.upper() in [
        "A",
        "C",
        "G",
        "T",
    ]


def print_read_stuff(
    pileupread: pysam.PileupRead, pileupcolumn, YN_support: str = None
):
    """
    Debug support func to print read and positional info information.
    """
    read = pileupread.alignment
    print(f"Read name: {read.query_name}")
    print(f"read fwd: {read.is_forward}")
    print(f"Read start: {read.reference_start}")
    print(f"Read CIGAR: {read.cigartuples}")
    print(f"Reference pos: {pileupcolumn.reference_pos}")
    print(f"Query pos: {pileupread.query_position}")
    print(f"Query seq (start): {read.query_sequence[:10]}")
    print(f"YN support: {YN_support if YN_support else 'None'}")
    print(f"Base: {pileupread.alignment.query_sequence[pileupread.query_position]}")
    print("---")


def plot_variant(df: pd.DataFrame, vcf_entry: pd.Series):
    # output_file("my_plot.html")

    # Color setup
    igv_colors = {
        "A": "#00C000",  # Green
        "C": "#0000C0",  # Blue
        "G": "#FFA500",  # Orange
        "T": "#FF0000",  # Red
        "D": "#808080",  # Gray for deletion
    }
    categories = list(igv_colors.keys())
    palette = [igv_colors[base] for base in categories]
    color_mapper = CategoricalColorMapper(factors=categories, palette=palette)

    source = ColumnDataSource(df)

    # Main scatter plot
    p_main = figure(
        title=f"Variant Support Scatter Plot for {vcf_entry.CHROM}:{vcf_entry.POS} ({vcf_entry.REF.upper()} -> {vcf_entry.ALT.upper()})",
        x_axis_label="Number of repeats in read position",
        y_axis_label="Support for the variant",
        width=800,
        height=600,
        x_range=(0, 25),
        y_range=(-0.02, 1.02),
        tools="pan,wheel_zoom,box_zoom,reset",
    )

    p_main.circle(
        x=jitter("total_count", width=0.45),
        y=jitter("alt_frac_all", width=0.02),
        source=source,
        size=6,
        alpha=0.6,
        color={"field": "YN_base", "transform": color_mapper},
        legend_field="YN_base",
    )
    p_main.legend.title = "YN_base"
    p_main.legend.location = "top_right"

    # X histogram
    hist_x, edges_x = np.histogram(df["total_count"], bins=25, range=(0, 25))
    p_hist_x = figure(
        width=800, height=150, x_range=p_main.x_range, tools="", toolbar_location=None
    )
    p_hist_x.quad(
        top=hist_x,
        bottom=0,
        left=edges_x[:-1],
        right=edges_x[1:],
        fill_color="gray",
        line_color="white",
    )
    p_hist_x.xaxis.visible = False
    p_hist_x.yaxis.axis_label = "Frequency"

    # Y histogram
    hist_y, edges_y = np.histogram(
        df["alt_frac_all"], bins=53, range=(-0.02, 1.02), density=True
    )
    hist_y = hist_y / np.sum(hist_y)
    hist_y_display = np.where(hist_y == 0, 1e-6, hist_y)
    # Create Bokeh figure
    p_hist_y = figure(
        width=150,
        height=600,
        y_range=(-0.02, 1.02),  # same as main plot Y
        x_range=(1e-4, 1.4),  # same as main plot X
        toolbar_location=None,
        x_axis_type="log",
        title="Y-axis Histogram (normalized + log1p)",
    )

    p_hist_y.quad(
        left=1e-6,
        right=hist_y_display,
        bottom=edges_y[:-1],
        top=edges_y[1:],
        fill_color="gray",
        line_color="white",
    )
    p_hist_y.yaxis.visible = False
    p_hist_y.yaxis.axis_label = "alt_frac_all"
    p_hist_y.xaxis.axis_label = "log1p(freq density)"

    # Combine plots using gridplot
    layout = gridplot([[p_hist_x, None], [p_main, p_hist_y]], merge_tools=False)
    return layout


def get_read_support_snp(vcf_entry: pd.Series, bam: pysam.AlignmentFile):
    """
    Get the read support for a given variant.
    """

    contig = vcf_entry.CHROM
    start = vcf_entry.POS - 1
    stop = vcf_entry.POS
    print(f"Plotting variant at {contig}:{start}-{stop}")
    read_support = []
    for pileupcolumn in bam.pileup(contig, start, stop, stepper="all", truncate=True):
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del:
                print(f"Deletion at position {pileupcolumn.reference_pos}")
            if pileupread.query_position is None:
                print(
                    f"Query position is None for base at position {pileupcolumn.reference_pos}"
                )
                continue

            YN_support = get_relevant_YN_support(
                pileupread.alignment.get_tag("YN"),
                pileupread.query_position,
                reverse=not pileupread.alignment.is_forward,
            )
            read_support.append(
                [
                    pileupread.alignment.query_name,
                    pileupread.alignment.is_forward,
                    pileupread.alignment.get_tag("YM"),
                    YN_support[0][0],
                    YN_support[1:],
                ]
            )

    print(f"found {len(read_support)} reads")
    return read_support


def get_relevant_YN_support_deletion(
    YN: str,
    query_position_next: int,
    reverse: bool,
    entry_split="|",
    subentry_split=",",
) -> str:
    """ALignment based consensus methods will have a different YN structure for deletions than sequence based methods."""
    deletion_support = []
    progress = 0
    for tag in YN.split(entry_split):
        if tag.startswith("D"):
            deletion_support.append(tag)
        else:
            progress += 1
            print(
                f"Progress: {progress} query_position_next: {query_position_next} D support was {deletion_support}"
            )
            if progress == query_position_next:
                return deletion_support
            else:
                deletion_support = []

    pass


def _get_read_support_deletion(vcf_entry: pd.Series, bam: pysam.AlignmentFile):
    contig = vcf_entry.CHROM
    # deletions are shown as ATC -> A if TC is deleted, so we need to start at the next base
    start = vcf_entry.POS
    stop = vcf_entry.POS + 1
    print(f"Plotting variant at {contig}:{start}-{stop}")
    read_support = []
    for pileupcolumn in bam.pileup(contig, start, stop, stepper="all", truncate=True):
        for read in pileupcolumn.pileups:
            if (
                read.alignment.query_name
                == "7f07dee5-7ec9-4260-a32a-535d0ea55928_I_0_chr7:55174724:55174833"
            ):
                print()
            if read.is_del:
                print(f"Deletion at position {pileupcolumn.reference_pos}")
                get_del_YN_support = get_relevant_YN_support_deletion(
                    read.alignment.get_tag("YN"),
                    read.query_position_or_next,
                    reverse=not read.alignment.is_forward,
                )
            if read.query_position is None:
                print(
                    f"Query position is None for base at position {pileupcolumn.reference_pos}"
                )
            else:
                # normal base at pos
                print(
                    f"Base at position {pileupcolumn.reference_pos} is {read.alignment.query_sequence[read.query_position]}"
                )
                read_support.append(
                    get_relevant_YN_support(
                        read.alignment.get_tag("YN"),
                        read.query_position,
                        reverse=not read.alignment.is_forward,
                    )
                )
            entry = []

    return read_support


def main(
    vcf_path: Path,
    bam_path: Path,
    output_path: Path,
    output_plot_file: Path = "output.html",
):
    """
    Main function to plot the variant support scatter plot.
    """
    # Tab info:
    tab_name = "Variant support scatter plots"
    add_info = {}
    json_obj = {}
    json_obj[tab_name] = {}
    json_obj[tab_name]["name"] = tab_name
    json_obj[tab_name]["priority"] = TAB_PRIORITY

    # Step 1: Collect all files
    vcf_path = Path(vcf_path)
    bam_file = Path(bam_path)
    output_path = Path(output_path)

    # Step 2:convert paths into related objects
    bam = pysam.AlignmentFile(bam_file, "rb")
    vcf = read_vcf(vcf_path)  # pandas df

    plots = []

    for index, vcf_entry in vcf.iterrows():
        variant_is_snp = determine_variant_is_snp(vcf_entry)
        if variant_is_snp:
            print(
                f"Plotting SNP variant at {vcf_entry.CHROM}:{vcf_entry.POS} ({vcf_entry.REF.upper()} -> {vcf_entry.ALT.upper()})"
            )
            read_support = get_read_support_snp(vcf_entry, bam)
            df = pd.DataFrame(
                read_support,
                columns=["read_name", "fwd", "YM", "YN_base", "YN_support"],
            )
            df = count_base_support(df)
            df = _add_reference_information_snp(df, vcf_entry)
            plots.append(plot_variant(df, vcf_entry))
        elif len(vcf_entry.REF) > len(vcf_entry.ALT):
            print(
                f"Skipping plotting deletion variant at {vcf_entry.CHROM}:{vcf_entry.POS} ({vcf_entry.REF.upper()} -> {vcf_entry.ALT.upper()})"
            )
            # read_support = _get_read_support_deletion(vcf_entry, bam)
            # df = _add_reference_information_deletion(df,vcf_entry)
        elif len(vcf_entry.REF) < len(vcf_entry.ALT):
            print(
                f"Skipping  plotting insertion variant at {vcf_entry.CHROM}:{vcf_entry.POS} ({vcf_entry.REF.upper()} -> {vcf_entry.ALT.upper()})"
            )
            # df = _add_`reference_information_insert(df,vcf_entry)
        else:
            exit(1, f"Unknown variant type for {vcf_entry.REF} -> {vcf_entry.ALT}")

    if len(plots) == 0:
        json_obj[tab_name]["script"], json_obj[tab_name]["div"] = (
            "",
            "<h1>No variants where plotable, indels are not plotted.</h1>",
        )
        return json_obj

    final_plot = column(*plots)
    final_plot = plots[0]
    json_obj[tab_name]["script"], json_obj[tab_name]["div"] = components(final_plot)
    json_obj["additional_info"] = add_info

    print("writing json")
    print(Path(output_plot_file))
    with open(Path(output_plot_file), "w") as f:
        f.write(json.dumps(json_obj))


if __name__ == "__main__":
    # test_vcf = "/home/dami/Software/cyclomicsseq/vis_test_data_cyc000492_bc2/variants/annotated/gDNA_Healthy_Controls_RCA_1_BF40c5B_barcode02_annotated.vcf"
    # test_bam = "/home/dami/Software/cyclomicsseq/vis_test_data_cyc000492_bc2/consensus_aligned/gDNA_Healthy_Controls_RCA_1_BF40c5B_barcode02.YM_gt_3.bam"
    # test_output = "output.html"
    # main(test_vcf, test_bam, test_output)

    # test_output = "output_az4.html"
    # main(test_vcf, test_bam, test_output)

    # Uncomment the following lines to run the script from the command line
    import argparse

    parser = argparse.ArgumentParser(description="Plot variant support scatter plot.")
    parser.add_argument("vcf_path", type=str, help="Path to the VCF file.")
    parser.add_argument("bam_path", type=str, help="Path to the BAM file.")
    parser.add_argument("output_path", type=str, help="Path to save the output plot.")

    args = parser.parse_args()

    main(args.vcf_path, args.bam_path, args.output_path)
