#!/usr/bin/env python

import argparse
from pathlib import Path
import yaml

import pandas as pd

from vcf.vcf_file import Vcf


def parse_arguments():
    """
    Argument parser.
    """
    parser = argparse.ArgumentParser(description="Filter a VCF")

    parser.add_argument(
        "-i",
        "--input_vcf",
        type=Path,
        required=True,
        help="Input VCF file with variant evidence to filter.",
    )
    parser.add_argument(
        "-o",
        "--output_vcf",
        type=Path,
        required=True,
        help="Output VCF file to which passed variants will be written.",
    )
    parser.add_argument(
        "-p",
        "--perbase_table",
        type=Path,
        required=True,
        help="Input Perbase TSV table with positional depths.",
    )
    parser.add_argument(
        "-p",
        "--dynamic_vaf_params",
        type=Path,
        required=True,
        help="Input file with params for Dynamic VAF filtering.",
    )
    parser.add_argument(
        "--min_dir_ratio",
        type=float,
        required=False,
        default=0.001,
        help="Minimum ratio of variant-supporting reads in each direction (default: 0.001).",
    )
    parser.add_argument(
        "--min_dir_count",
        type=float,
        required=False,
        default=5,
        help="Minimum number of variant-supporting reads in each direction (default: 50).",
    )
    parser.add_argument(
        "--min_dpq",
        type=float,
        required=False,
        default=5_000,
        help="Minimum positional depth after Q filtering (default: 1_000).",
    )
    parser.add_argument(
        "--min_dpq_n",
        type=int,
        required=False,
        default=25,
        help="Number of flanking nucleotides to the each position that will determine the window size for local maxima calculation (default = 4).",
    )
    parser.add_argument(
        "--min_dpq_ratio",
        type=float,
        required=False,
        default=0.3,
        help="Ratio of local depth maxima that will determine the minimum depth at each position (default = 0.3).",
    )
    parser.add_argument(
        "--min_vaf",
        type=float,
        required=False,
        default=0.003,
        help="Minimum variant allele frequency (default: 0.002).",
    )
    parser.add_argument(
        "--min_rel_ratio",
        type=float,
        required=False,
        default=0.3,
        help="Minimum relative ratio between forward and reverse variant-supporting reads (default: 0.3).",
    )
    parser.add_argument(
        "--min_abq",
        type=float,
        required=False,
        default=70,
        help="Minimum average base quality (default: 70).",
    )

    return parser.parse_args()


def get_depth_table(perbase_tsv: Path) -> pd.DataFrame:
    """
    Read in a perbase table and returns it as a pandas DataFrame with reference, position and depth.
    """
    try:
        df = pd.read_table(perbase_tsv, sep="\t")
    except pd.errors.EmptyDataError:
        return None
    return df[["REF", "POS", "DEPTH"]]


def get_dynamic_vaf_params(params_yml: Path) -> dict:
    with open(params_yml, "r") as file:
        params = yaml.safe_load(file)

    return params


if __name__ == "__main__":
    dev = False
    if not dev:
        args = parse_arguments()

        depth_table = get_depth_table(args.perbase_table)
        dynamic_vaf_params = get_dynamic_vaf_params(args.dynamic_vaf_params)

        vcf = Vcf(args.input_vcf)
        vcf.filter(
            depth_table,
            args.min_dir_ratio,
            args.min_dir_count,
            args.min_dpq,
            args.min_dpq_n,
            args.min_dpq_ratio,
            args.min_vaf,
            args.min_rel_ratio,
            args.min_abq,
        )
        vcf.write(args.output_vcf)

    if dev:
        depth_table = get_depth_table("./L7154-000200.tsv")
        dynamic_vaf_params = get_dynamic_vaf_params("./bin/vcf/dynamic_vaf_params.yml")
        vcf = Vcf("./200.norm.vcf")
        vcf.filter(depth_table, dynamic_vaf_params)
        vcf.write("./filtered_200.norm.vcf")
