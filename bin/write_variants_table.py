#!/usr/bin/env python
import json
from pathlib import Path
import pandas as pd
import io
import re


def load_vcf(vcf_file: Path) -> pd.DataFrame:
    """Loads a VCF file as a Pandas DataFrame."""
    with open(vcf_file, "r") as f:
        header = []
        lines = []
        for l in f:
            if l.startswith("##"):
                header.append(l)
            else:
                lines.append(l)

    vcf_header = header
    variants_df = pd.read_csv(
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
            "Sample1": str,
        },
        sep="\t",
    ).rename(columns={"#CHROM": "CHROM"})

    return variants_df


def restructure_annotations(variants_df: pd.DataFrame) -> pd.DataFrame:
    """Restructures variants dataframe to have readable annotations."""
    chrom = variants_df["CHROM"]
    pos = variants_df["POS"]
    location = pd.Series([f"{c}:{p}" for c, p in zip(chrom, pos)])

    ref = variants_df["REF"]
    alt = variants_df["ALT"]

    sample1 = variants_df["Sample1"].str.split(":")
    vaf = sample1.str[3]

    info = variants_df["INFO"].str.split(";")

    #  if we only have vcf files with empty INFO columns (eg only backbone variants, or non cosmic mutations)
    if len(info.shape) == 1:
        type = pd.Series(["N/A"] * len(location))
        consequence = pd.Series(["N/A"] * len(location))
        symbol = pd.Series(["N/A"] * len(location))
        impact = pd.Series(["N/A"] * len(location))
        biotype = pd.Series(["N/A"] * len(location))
        sift = pd.Series(["N/A"] * len(location))
        polyphen = pd.Series(["N/A"] * len(location))
        cosmic_ids = pd.Series(["N/A"] * len(location))
        
    else:
        type = info.str[0].str.split("=").str[1]
        consequence = info.str[1].str.split("=").str[1]
        symbol = info.str[4].str.split("=").str[1]
        impact = info.str[5].str.split("=").str[1]
        biotype = info.str[6].str.split("=").str[1]
        sift = info.str[9].str.split("=").str[1]
        polyphen = info.str[10].str.split("=").str[1]
        cosmic_ids = info.str[2].str.split("=").str[1]

    annot_columns = [
        "Location",
        "Ref",
        "Alt",
        "VAF",
        "Type",
        "Symbol",
        "Biotype",
        "Consequence",
        "Impact",
        "SIFT",
        "PolyPhen",
        "COSMIC",
    ]

    annot_data = [
        location,
        ref,
        alt,
        vaf,
        type,
        symbol,
        biotype,
        consequence,
        impact,
        sift,
        polyphen,
        cosmic_ids,
    ]

    annotation_df = pd.concat(annot_data, axis=1)
    annotation_df.columns = annot_columns
    annotation_df["COSMIC"] = annotation_df["COSMIC"].replace(
        to_replace=",", value=", ", regex=True
    )
    annotation_df = annotation_df.replace(to_replace=r"\(\)", value="")
    annotation_df = annotation_df.replace(to_replace=r"None \(None\)", value="N/A")
    annotation_df = annotation_df.replace(to_replace="None", value="N/A")
    annotation_df = annotation_df.replace(to_replace=r"nan \(nan\)", value="N/A")
    annotation_df = annotation_df.replace(to_replace="nan", value="N/A")
    annotation_df = annotation_df.replace(to_replace="NaN", value="N/A")

    return annotation_df


def main(vcf_file: Path, variant_table_file: Path, tab_name: str):
    """ """
    variants_df = load_vcf(vcf_file)
    annotation_df = restructure_annotations(variants_df)
    vcf_table = annotation_df.to_html(na_rep="N/A", escape=False)

    with open(Path(variant_table_file).with_suffix(".json"), "w") as f:
        json_obj = {}
        json_obj[tab_name] = {}
        json_obj[tab_name]["name"] = tab_name
        json_obj[tab_name]["script"] = ""
        json_obj[tab_name]["div"] = vcf_table

        f.write(json.dumps(json_obj))


if __name__ == "__main__":
    dev = False
    if not dev:
        import argparse

        parser = argparse.ArgumentParser(
            description="Create hist plot from a regex for fastq and fastq.gz files."
        )

        parser.add_argument("vcf_file")
        parser.add_argument("variant_table_file")
        parser.add_argument("tab_name", default="Variant Table")
        args = parser.parse_args()

        main(args.vcf_file, args.variant_table_file, args.tab_name)

    else:
        # vcf_file = "/scratch/projects/ROD_0908_63_variantcalling/results/PR_test/variants/FAS12641_annotated.vcf"
        vcf_file = "ABZ922_annotated.vcf"
        variant_table_file = "variant_table.json"
        tab_name = "variant_table"

        main(vcf_file, variant_table_file, tab_name)
