

import io
import os
import json
import datetime
from collections import defaultdict
import math
import glob

from pathlib import Path
from subprocess import check_output
from collections import Counter, defaultdict
from typing import Any


import pandas as pd
import numpy as np
import pysam

from bokeh.io import output_notebook
from bokeh.plotting import figure, show, output_file
from bokeh.models import Range1d, ColumnDataSource,CategoricalColorMapper
from bokeh.palettes import Category10
from bokeh.transform import jitter
from bokeh.io import export_png

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

def get_relevant_YN_support(YN:str, query_position:int, reverse:bool, entry_split = "|", subentry_split = ",") -> str:
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
            # print(f"YN_pos: {YN_pos}/{i} at {query_position} YN_skip: {YN_skip} sub_struct: {sub_struct}")
            # print(f"sub_struct: {YN_support[query_position -1]}")
            # return sub_struct
            relevant_YN_struct = sub_struct
            break
        
        YN_pos +=1

    # base_translator = str.maketrans("ACGTacgt", "tgcaTGCA")
    # if read is reversed, we need to flip all the bases to their complemntary bases, but keep their capitalization for reverse
    # if reverse and relevant_YN_struct:
        # relevant_YN_struct = [x.translate(base_translator) for x in relevant_YN_struct]
        
    if not relevant_YN_struct:
        # print(f"relevant_YN_struct: {relevant_YN_struct} YN_pos: {YN_pos} query_position: {query_position} YN_skip: {YN_skip}")
        # print(f"sub_struct: {YN_support[query_position -1]}")
        relevant_YN_struct = ""
    return relevant_YN_struct 


def count_base_support(df:pd.DataFrame, letters = ['A', 'C', 'D', 'G', 'T']) -> pd.DataFrame:
    # inline func to apply to the YN_support column
    def count_initials(lst, letter):
        return sum(1 for x in lst if x and x[0].upper() == letter)
    
    # Apply the function per row
    for letter in letters:
        df[f"{letter}_count"] = df["YN_support"].apply(lambda lst: count_initials(lst, letter))

    return df

def add_reference_information(df:pd.DataFrame, vcf_entry:str) -> pd.DataFrame:
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
    
    #fix the NaN values due to division by zero
    df.fillna(0, inplace=True)
    return df

def convert_to_df(read_support:list, vcf_entry: pd.Series) -> pd.DataFrame:
    """
    Convert the read support list to a pandas DataFrame.
    """
    df = pd.DataFrame(read_support, columns=["read_name","fwd", "YM", "YN_base", "YN_support"])
    
    df = count_base_support(df)
    df = add_reference_information(df,vcf_entry)
    return df

def print_read_stuff(pileupread:pysam.PileupRead,pileupcolumn, YN_support:str = None):
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

def _plot_variant(df:pd.DataFrame, vcf_entry: pd.Series):
    ref_base = vcf_entry.REF.upper()
    alt_base = vcf_entry.ALT.upper()

    contig = vcf_entry.CHROM
    position = vcf_entry.POS
    output_file("my_plot.html")
    # p = figure(title="Plot", width=400, height=400)
    # p.line([1, 2, 3], [4, 6, 2])
    p = figure(title=f"Variant Support Scatter Plot for {contig}:{position} ({ref_base} -> {alt_base})", 
               x_axis_label='Position', 
               y_axis_label='Support for the variant',
               width=1200, 
               height=600,
                x_range=(0, 25),
                y_range=(-0.02, 1.02),
               )
    # Prepare color mapper

    igv_colors = {
    'A': '#00C000',  # Green
    'C': '#0000C0',  # Blue
    'G': '#FFA500',  # Orange
    'T': '#FF0000',  # Red
    'D': '#808080',  # Gray for deletion
    }
    categories = list(igv_colors.keys())
    palette = [igv_colors[base] for base in categories]

    color_mapper = CategoricalColorMapper(factors=categories, palette=palette)
    source = ColumnDataSource(df)
    # p.line([1, 2, 3], [4, 6, 2], line_width=2)  # Adds a line glyph
    p.circle(x=jitter('total_count', width=0.45), 
             y=jitter('alt_frac_all', width=0.02), 
             source=source, 
             size=6, 
             alpha=0.6,
             color={'field': 'YN_base', 'transform': color_mapper},
             legend_field='YN_base',

             )
    p.legend.title = "YN_base"
    p.legend.location = "top_right"
    
    show(p)
    print('done')


def plot_variant(vcf_entry:pd.Series, bam: pysam.AlignmentFile):
    """
    Plot the variant support scatter plot for a given chromosome and position.
    """
    pass
    contig = vcf_entry.CHROM
    start = vcf_entry.POS - 1
    stop = vcf_entry.POS
    print(f"Plotting variant at {contig}:{start}-{stop}")
    read_support = []
    for pileupcolumn in bam.pileup(contig, start, stop, stepper="all", truncate=True):
        for pileupread in pileupcolumn.pileups:
            entry = []

            if pileupread.is_del:
                print(f"Deletion at position {pileupcolumn.reference_pos}")
                

            if pileupread.query_position is None:
                print(f"Query position is None for base at position {pileupcolumn.reference_pos}")
                # continue

            # print_read_stuff(pileupread=pileupread,pileupcolumn=pileupcolumn)
            YN_support = get_relevant_YN_support(pileupread.alignment.get_tag("YN"), pileupread.query_position,reverse= not pileupread.alignment.is_forward)
            # YN_support = get_relevant_YN_support(pileupread.alignment.get_tag("YN"), pileupread.query_position,reverse= False)
            # '1fc7aca5-22b3-4428-9e86-475312dfac1a_I_3_chr5:1295264:1295782'
            if pileupread.alignment.query_name == '1a07669d-fd34-4634-926c-8c07ed9f669c_I_0_chr5:1294641:1295766':
                print_read_stuff(pileupread=pileupread,pileupcolumn=pileupcolumn, YN_support=YN_support)
            read_support.append([pileupread.alignment.query_name,pileupread.alignment.is_forward, pileupread.alignment.get_tag("YM"), YN_support[0][0], YN_support[1:] ])
            pass
   
    print(f"found {len(read_support)} reads")
    df = convert_to_df(read_support, vcf_entry)

    _plot_variant(df, vcf_entry)

    print()

def main(vcf_path:Path,bam_path:Path,output_path:Path):
    """
    Main function to plot the variant support scatter plot.
    """
    # Step 1: Collect all files
    vcf_path = Path(vcf_path)
    bam_file = Path(bam_path)
    output_path = Path(output_path)

    # Step 2:convert paths into related objects
    bam = pysam.AlignmentFile(bam_file, "rb")
    vcf = read_vcf(vcf_path) # pandas df

    chrom = 'chr5'
    pos = 1295601
    for index, vcf_entry in vcf.iterrows():
        plot_variant(vcf_entry,bam)
        
    
    # # Step 4: Create a scatter plot
    # p = figure(title="Variant Support Scatter Plot", x_axis_label='Position', y_axis_label='Support')
    
    # # Step 5: Add data to the plot
    # for index, row in vcf_data.iterrows():
    #     pos = row['POS']
    #     support = bam_file.count(row['CHROM'], pos-1, pos)
    #     p.circle(pos, support, size=10, color="navy", alpha=0.5)
    
    # # Step 6: Save the plot
    # export_png(p, filename=str(output_path))



if __name__ == "__main__":
    test_vcf = "/home/dami/Software/cyclomicsseq/vis_test_data_cyc000492_bc2/variants/annotated/gDNA_Healthy_Controls_RCA_1_BF40c5B_barcode02_annotated.vcf"
    test_bam = "/home/dami/Software/cyclomicsseq/vis_test_data_cyc000492_bc2/consensus_aligned/gDNA_Healthy_Controls_RCA_1_BF40c5B_barcode02.YM_gt_3.bam"
    test_output = "output.html"
    main(test_vcf, test_bam, test_output)

    # Uncomment the following lines to run the script from the command line
    # import argparse

    # parser = argparse.ArgumentParser(description="Plot variant support scatter plot.")
    # parser.add_argument("vcf_path", type=str, help="Path to the VCF file.")
    # parser.add_argument("bam_path", type=str, help="Path to the BAM file.")
    # parser.add_argument("output_path", type=str, help="Path to save the output plot.")

    # args = parser.parse_args()

    # main(args.vcf_path, args.bam_path, args.output_path)