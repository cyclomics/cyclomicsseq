#!/usr/bin/env python
import io
import json
from pathlib import Path

import pandas as pd

TAB_PRIORITY = 70


def load_vcf(vcf_file: Path) -> pd.DataFrame:
    """Loads a VCF file as a Pandas DataFrame."""
    with open(vcf_file, "r") as fh:
        lines = []
        for line in fh:
            if line.startswith("##"):
                continue
            else:
                lines.append(line)

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
            "SAMPLE1": str,
        },
        sep="\t",
    ).rename(columns={"#CHROM": "CHROM"})

    return df


def restructure_annotations(
    df: pd.DataFrame, variant_decimal_points: int = 3
) -> pd.DataFrame:
    """Restructures variants dataframe to have readable annotations."""
    chrom = df["CHROM"]
    pos = df["POS"]
    location = pd.Series([f"{c}:{p}" for c, p in zip(chrom, pos)])

    ref = df["REF"]
    alt = df["ALT"]

    vaf = df["SAMPLE1"].str.split(":").str[3]
    vaf = (vaf.astype(float) * 100).round(variant_decimal_points).astype(str) + "%"

    columns = [
        "Location",
        "Ref",
        "Alt",
        "VAF (%)",
    ]

    data = [
        location,
        ref,
        alt,
        vaf,
    ]

    contaminant_table = pd.concat(data, axis=1)
    contaminant_table.columns = columns

    return contaminant_table


def main(
    vcf_file: Path, contaminant_table_file: Path, tab_name: str, priority_limit: int
):
    """
    Write contaminants in VCF file to a JSON file to be loaded into the HTML report.

    Args:
        vcf_file: Path to the VCF file containing the synthetic mutant contaminant
            information.
        contaminant_table_file: Name of the contaminant_table JSON file to write to
            from a pandas DataFrame.
        tab_name: Name of the variant table tab to add to the report, str.
    """
    json_obj = {}
    json_obj[tab_name] = {}
    json_obj[tab_name]["name"] = tab_name
    json_obj[tab_name]["priority"] = TAB_PRIORITY

    if TAB_PRIORITY < priority_limit:
        version_notice = "<br><br><p>Variant annotation currently only supported with human genome version GRCh38.p14 and with variant calling option '--variant_calling validate'</p>"

        vcf_df = load_vcf(vcf_file)
        contaminants_df = restructure_annotations(vcf_df)
        contaminant_table = contaminants_df.to_html(na_rep="N/A")
        contaminant_table = contaminant_table.replace(
            'class="dataframe"', 'class="table table-sm table-hover table-striped"'
        )
        contaminant_table = contaminant_table.replace('border="1"', "")

        json_obj[tab_name]["script"] = ""
        json_obj[tab_name]["div"] = (
            "<div>" + contaminant_table + version_notice + "</div>"
        )

    with open(Path(contaminant_table_file).with_suffix(".json"), "w") as f:
        f.write(json.dumps(json_obj))


if __name__ == "__main__":
    dev = False
    if not dev:
        import argparse

        parser = argparse.ArgumentParser(
            description="Create hist plot from a regex for fastq and fastq.gz files."
        )

        parser.add_argument("vcf_file")
        parser.add_argument("contaminant_table_file")
        parser.add_argument("tab_name", default="Contaminants")
        parser.add_argument("priority_limit", type=int, default=89)
        args = parser.parse_args()

        main(
            args.vcf_file,
            args.contaminant_table_file,
            args.tab_name,
            args.priority_limit,
        )

    else:
        vcf_file = "testing_37_TP53_fourMutSet/variants/fastq_filtered_annotated_contaminants.vcf"
        contaminant_table_file = "contaminant_table.json"
        tab_name = "Contaminants"
        priority_limit = 9999

        main(vcf_file, contaminant_table_file, tab_name, priority_limit)
