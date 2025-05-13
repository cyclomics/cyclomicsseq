#!/usr/bin/env python
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


def get_unique_values(lst):
    # Flatten one level and filter out unhashable items (e.g. inner lists)
    flat = (item for sub in lst for item in (sub if isinstance(sub, list) else [sub]))

    # Filter out None and unhashable items, keep only hashables
    unique_vals = {
        item
        for item in flat
        if item is not None and isinstance(item, (str, int, float, bool))
    }

    return ", ".join(unique_vals) if unique_vals else None


def restructure_annotations(
    variants_df: pd.DataFrame, variant_decimal_points=3
) -> pd.DataFrame:
    """Restructures variants dataframe to have readable annotations."""
    chrom = variants_df["CHROM"]
    pos = variants_df["POS"]
    location = pd.Series(
        [f"{c}:{p}" for c, p in zip(chrom, pos)], name="LOC", dtype=object
    )

    ref = variants_df["REF"]
    alt = variants_df["ALT"]

    sample1 = variants_df["SAMPLE1"].str.split(":")
    vaf = sample1.str[3]
    # convert fraction to percentage
    vaf = (vaf.astype(float) * 100).round(variant_decimal_points).astype(str) + "%"

    coverage = sample1.str[0]
    coverage = coverage.apply(human_format)

    info = variants_df["INFO"].str.split(";")

    if (info.str[0] == "ANNOTATION").all() and (info.str[1] != ".").all():
        var_type = info.str[1].str.split("=").str[1]
        consequence = info.str[2].str.split("=").str[1]
        signal = (
            info.str[3].str.split("=").str[1].str.split(",").apply(get_unique_values)
        )
        clinvar = info.str[4].str.split("=").str[1]
        cosmic_ids = info.str[5].str.split("=").str[1]
        legacy_ids = info.str[6].str.split("=").str[1]
        symbol = info.str[7].str.split("=").str[1]
        impact = info.str[8].str.split("=").str[1]
        biotype = info.str[9].str.split("=").str[1]
        aa_change = info.str[10].str.split("=").str[1]
        sift = info.str[12].str.split("=").str[1]
        polyphen = info.str[13].str.split("=").str[1]

    else:
        var_type = pd.Series(["N/A"] * len(location))
        consequence = pd.Series(["N/A"] * len(location))
        signal = pd.Series(["N/A"] * len(location))
        clinvar = pd.Series(["N/A"] * len(location))
        cosmic_ids = pd.Series(["N/A"] * len(location))
        legacy_ids = pd.Series(["N/A"] * len(location))
        symbol = pd.Series(["N/A"] * len(location))
        impact = pd.Series(["N/A"] * len(location))
        biotype = pd.Series(["N/A"] * len(location))
        aa_change = pd.Series(["N/A"] * len(location))
        sift = pd.Series(["N/A"] * len(location))
        polyphen = pd.Series(["N/A"] * len(location))

    annot_columns = [
        "Location",
        "Ref",
        "Alt",
        "Var (%)",
        "Coverage",
        "Type",
        "Symbol",
        "Biotype",
        "AA change",
        "Consequence",
        "Impact",
        "SIFT",
        "PolyPhen",
        "Signal",
        "ClinVar",
        "COSMIC",
        "COSMIC legacy",
    ]

    annot_data = [
        location,
        ref,
        alt,
        vaf,
        coverage,
        var_type,
        symbol,
        biotype,
        aa_change,
        consequence,
        impact,
        sift,
        polyphen,
        signal,
        clinvar,
        cosmic_ids,
        legacy_ids,
    ]

    annotation_df = pd.concat(annot_data, axis=1)
    annotation_df.columns = annot_columns
    annotation_df["ClinVar"] = annotation_df["ClinVar"].replace(
        to_replace=",", value=", ", regex=True
    )
    annotation_df["COSMIC"] = annotation_df["COSMIC"].replace(
        to_replace=",", value=", ", regex=True
    )
    annotation_df["COSMIC legacy"] = annotation_df["COSMIC legacy"].replace(
        to_replace=",", value=", ", regex=True
    )
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
        version_notice = "<br><br><p>Variant annotation currently only supported with human genome version GRCh38.p14 and with variant calling option '--variant_calling validate'</p>"

        variants_df = load_vcf(vcf_file)
        annotation_df = restructure_annotations(variants_df)
        vcf_table = annotation_df.to_html(na_rep="N/A")
        vcf_table = vcf_table.replace(
            'class="dataframe"', 'class="table table-sm table-hover table-striped"'
        )
        vcf_table = vcf_table.replace('border="1"', "")
        # f"width={cyclomics_defaults.width}")

        json_obj[tab_name]["script"] = ""
        json_obj[tab_name]["div"] = "<div>" + vcf_table + version_notice + "</div>"

    with open(Path(variant_table_file).with_suffix(".json"), "w") as f:
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
        parser.add_argument("priority_limit", type=int, default=89)
        args = parser.parse_args()

        main(args.vcf_file, args.variant_table_file, args.tab_name, args.priority_limit)

    else:
        vcf_file = "/data/projects/ROD_tmp/65/08147aed79107ea52881035e4be36a/FAW08675_filtered_annotated.vcf"
        variant_table_file = "variant_table.json"
        tab_name = "variant_table"
        priority_limit = 9999

        main(vcf_file, variant_table_file, tab_name, priority_limit)
