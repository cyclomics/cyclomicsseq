#!/usr/bin/env python
import argparse
import io
import json
from pathlib import Path

import pandas as pd
from plotting_defaults import human_format

TAB_PRIORITY = 2


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
            "SAMPLE1": str,
        },
        sep="\t",
    ).rename(columns={"#CHROM": "CHROM"})

    return variants_df


def restructure_annotations(
    variants_df: pd.DataFrame, variant_decimal_points=2
) -> pd.DataFrame:
    """Restructures variants dataframe to have readable annotations."""
    chrom = variants_df["CHROM"]
    pos = variants_df["POS"]
    location = pd.Series([f"{c}:{p}" for c, p in zip(chrom, pos)])

    ref = variants_df["REF"]
    alt = variants_df["ALT"]

    sample1 = variants_df["SAMPLE1"].str.split(":")
    vaf = sample1.str[3]
    # convert fraction to percentage
    vaf = (vaf.astype(float) * 100).round(variant_decimal_points).astype(str) + "%"

    coverage = sample1.str[0]
    coverage = coverage.apply(human_format)

    info = variants_df["INFO"].str.split(";")

    # Mask: row has annotation
    annot_mask = (info.str[0] == "ANNOTATION") & (info.str[1] != ".")

    # Default values (N/A)
    var_type = pd.Series(["N/A"] * len(info))
    consequence = pd.Series(["N/A"] * len(info))
    symbol = pd.Series(["N/A"] * len(info))
    biotype = pd.Series(["N/A"] * len(info))
    clinsig = pd.Series(["N/A"] * len(info))
    hgvsc = pd.Series(["N/A"] * len(info))
    hgvsp = pd.Series(["N/A"] * len(info))
    cosmic_ids = pd.Series(["N/A"] * len(info))

    # Fill only for annotated rows
    var_type.loc[annot_mask] = info.loc[annot_mask].str[1].str.split("=").str[1]
    consequence.loc[annot_mask] = info.loc[annot_mask].str[2].str.split("=").str[1]
    cosmic_ids.loc[annot_mask] = info.loc[annot_mask].str[3].str.split("=").str[1]
    symbol.loc[annot_mask] = info.loc[annot_mask].str[5].str.split("=").str[1]
    biotype.loc[annot_mask] = info.loc[annot_mask].str[7].str.split("=").str[1]
    clinsig.loc[annot_mask] = info.loc[annot_mask].str[12].str.split("=").str[1]
    hgvsc.loc[annot_mask] = info.loc[annot_mask].str[13].str.split("=").str[1]
    hgvsp.loc[annot_mask] = info.loc[annot_mask].str[14].str.split("=").str[1]

    annot_columns = [
        "Location",
        "Ref",
        "Alt",
        "VAF (%)",
        "Coverage",
        "Symbol",
        "Type",
        "Codon change",
        "Aminoacid change",
        "Biotype",
        "Consequence",
        "Clinical significance",
        "COSMIC",
    ]

    annot_data = [
        location,
        ref,
        alt,
        vaf,
        coverage,
        symbol,
        var_type,
        hgvsc,
        hgvsp,
        biotype,
        consequence,
        clinsig,
        cosmic_ids,
    ]

    annotation_df = pd.concat(annot_data, axis=1)
    annotation_df.columns = annot_columns
    annotation_df["COSMIC"] = annotation_df["COSMIC"].replace(",", ", ", regex=True)
    annotation_df["Biotype"] = annotation_df["Biotype"].replace("_", " ", regex=True)
    annotation_df["Consequence"] = (
        annotation_df["Consequence"]
        .replace("_variant", "", regex=True)
        .replace("_", " ", regex=True)
    )
    annotation_df["Clinical significance"] = annotation_df[
        "Clinical significance"
    ].replace("_", " ", regex=True)
    annotation_df = annotation_df.replace(to_replace=r"\(\)", value="")
    annotation_df = annotation_df.replace(to_replace=r"None \(None\)", value="N/A")
    annotation_df = annotation_df.replace(to_replace="None", value="N/A")
    annotation_df = annotation_df.replace(to_replace=r"nan \(nan\)", value="N/A")
    annotation_df = annotation_df.replace(to_replace="nan", value="N/A")
    annotation_df = annotation_df.replace(to_replace="NaN", value="N/A")

    return annotation_df


def main(vcf_file: Path, variant_table_file: Path, tab_name: str, priority_limit: int):
    """
    Write variants in VCF file to a JSON file to be loaded into the HTML report.

    Args:
        vcf_file: Path to the VCF file containing the variant information.
        variant_table_file: Name of the variant_table JSON file to write to
            from a pandas DataFrame.
        tab_name: Name of the variant table tab to add to the report, str.
    """
    json_obj = {}
    json_obj[tab_name] = {}
    json_obj[tab_name]["name"] = tab_name
    json_obj[tab_name]["priority"] = TAB_PRIORITY

    if TAB_PRIORITY < priority_limit:
        version_notice = "<br><br><p>Variant annotation currently only supported with human genome version GRCh38.p14</p>"

        variants_df = load_vcf(vcf_file)
        annotation_df = restructure_annotations(variants_df)
        vcf_table = annotation_df.to_html(na_rep="N/A", table_id="variant_table")
        vcf_table = vcf_table.replace(
            'class="dataframe"',
            'class="table table-striped table-sm table-hover table-responsive"',
        )
        vcf_table = vcf_table.replace('border="1"', "")
        # f"width={cyclomics_defaults.width}")

        json_obj[tab_name]["script"] = (
            "<script>"
            "$(document).ready(function() {"
            "$('#variant_table').DataTable({"
            "paging: true, searching: true, ordering: true, scrollX: true"
            "});"
            "});"
            "</script>"
        )
        json_obj[tab_name]["div"] = (
            '<div style="font-size: 12px; max-width:100%; overflow-x:auto; white-space:nowrap;">'
            + vcf_table
            + version_notice
            + "</div>"
        )

    with open(Path(variant_table_file).with_suffix(".json"), "w") as f:
        f.write(json.dumps(json_obj))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create hist plot from a regex for fastq and fastq.gz files."
    )

    parser.add_argument("vcf_file")
    parser.add_argument("variant_table_file")
    parser.add_argument("tab_name", default="Variant Table")
    parser.add_argument("priority_limit", type=int, default=89)
    args = parser.parse_args()

    main(args.vcf_file, args.variant_table_file, args.tab_name, args.priority_limit)
