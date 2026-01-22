#!/usr/bin/env python

import argparse
import io
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

FORMAT = [
    "DP",
    "DPQ",
    "FREQ",
    "VAF",
    "FWDC",
    "FWDR",
    "REVC",
    "REVR",
    "TOTC",
    "TOTR",
    "SAME",
    "OBSR",
    "ABQ",
    "OBQ",
    "HCR",
]

FORMAT_HEADERS = [
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n',
    '##FORMAT=<ID=DPQ,Number=1,Type=Integer,Description="Depth quality">\n',
    '##FORMAT=<ID=FREQ,Number=1,Type=Float,Description="Reference allele fraction">\n',
    '##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele fraction">\n',
    '##FORMAT=<ID=FWDC,Number=1,Type=Integer,Description="Forward variant count">\n',
    '##FORMAT=<ID=FWDR,Number=1,Type=Float,Description="Forward variant ratio">\n',
    '##FORMAT=<ID=REVC,Number=1,Type=Integer,Description="Reverse variant count">\n',
    '##FORMAT=<ID=REVR,Number=1,Type=Float,Description="Reverse variant ratio">\n',
    '##FORMAT=<ID=TOTC,Number=1,Type=Integer,Description="Total variant count">\n',
    '##FORMAT=<ID=TOTR,Number=1,Type=Float,Description="Total variant ratio">\n',
    '##FORMAT=<ID=SAME,Number=1,Type=String,Description="Strand agreement">\n',
    '##FORMAT=<ID=OBSR,Number=1,Type=String,Description="Observed strands">\n',
    '##FORMAT=<ID=ABQ,Number=1,Type=Float,Description="Alt base quality">\n',
    '##FORMAT=<ID=OBQ,Number=1,Type=Float,Description="Overall base quality">\n',
    '##FORMAT=<ID=HCR,Number=1,Type=String,Description="High-confidence region">\n',
]


class VcfFile:
    def __init__(self, vcf_file):
        self.vcf_header = ""
        self.vcf = self.read_vcf(vcf_file)

    def read_vcf(self, path: Path) -> pd.DataFrame:
        """Read in a VCF file and return it as Pandas DataFrame"""
        with open(path, "r") as f:
            header = []
            lines = []
            for l in f:
                if l.startswith("##"):
                    header.append(l)
                else:
                    lines.append(l)

        self.vcf_header = header
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
                "FORMAT": str,
                "SAMPLE1": str,
            },
            sep="\t",
        ).rename(columns={"#CHROM": "CHROM"})

        return df

    def write_vcf(self, output_path: Path):
        """Write output VCF file."""

        with open(output_path, "w") as new_vcf:
            new_header = [i for i in self.vcf_header if not i.startswith("##INFO")]
            new_header += FORMAT_HEADERS
            self.vcf_header = new_header
            for line in self.vcf_header:
                new_vcf.write(line)

            writeable_vcf = self.vcf.rename(columns={"CHROM": "#CHROM"})
            writeable_vcf = writeable_vcf[
                [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "INFO",
                    "FORMAT",
                    "SAMPLE1",
                ]
            ]

            # Write body
            new_vcf.writelines(writeable_vcf.to_csv(sep="\t", index=False))

    @staticmethod
    def extract_info_lofreq(info_str: str) -> Dict[str, str]:
        """Parse LoFreq INFO field and compute derived annotations.

        Args:
            info_str: INFO column string from VCF.

        Returns:
            Dictionary mapping FORMAT keys to string values.
        """
        raw: Dict[str, str] = {}
        for item in info_str.split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                raw[k] = v

        DP = int(raw.get("DP", 0))
        DPQ = DP
        AF = float(raw.get("AF", 0))
        FREQ = f"{1 - AF:.5f}"
        VAF = f"{AF:.5f}"

        dp4 = list(map(int, raw.get("DP4", "0,0,0,0").split(",")))
        FWDREF, REVREF, FWDC, REVC = dp4
        FWDR = FWDC / FWDREF if FWDREF else 0.0
        REVR = REVC / REVREF if REVREF else 0.0
        TOTC = FWDC + REVC
        TOTR = TOTC / (FWDREF + REVREF) if (FWDREF + REVREF) else 0.0

        SAME = raw.get("SAME", "0")
        OBSR = raw.get("OBSR", "0.0")
        ABQ = raw.get("ABQ", "0.0")
        OBQ = raw.get("OBQ", "0.0")
        HCR = raw.get("HCR", "0")

        return {
            "DP": str(DP),
            "DPQ": str(DPQ),
            "FREQ": str(FREQ),
            "VAF": str(VAF),
            "FWDC": str(FWDC),
            "FWDR": f"{FWDR:.5f}",
            "REVC": str(REVC),
            "REVR": f"{REVR:.5f}",
            "TOTC": str(TOTC),
            "TOTR": f"{TOTR:.5f}",
            "SAME": str(SAME),
            "OBSR": str(OBSR),
            "ABQ": str(ABQ),
            "OBQ": str(OBQ),
            "HCR": str(HCR),
        }

    def parse_columns(self, caller: str) -> pd.DataFrame:
        df = self.vcf.copy()

        if df.empty:
            return df

        if caller != "lofreq":
            return df

        def transform_row(info_str: str) -> Tuple[str, str, str]:
            vals = self.extract_info_lofreq(info_str)

            # FORMAT is fixed
            format_str = ":".join(FORMAT)

            # SAMPLE1 follows FORMAT order
            sample_vals: List[str] = [vals.get(k, ".") for k in FORMAT]
            sample_str = ":".join(sample_vals)

            return ".", format_str, sample_str

        new_cols = df["INFO"].apply(
            lambda s: pd.Series(transform_row(s), index=["INFO", "FORMAT", "SAMPLE1"])
        )

        df[["INFO", "FORMAT", "SAMPLE1"]] = new_cols

        self.vcf = df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Reformat VCF file to standardize annotations."
    )
    parser.add_argument(
        "caller",
        type=str,
        choices=["lofreq"],
        help="Variant caller used to generate the VCF file.",
    )
    parser.add_argument(
        "input_vcf",
        type=Path,
        help="Path to the input VCF file.",
    )
    parser.add_argument(
        "output_vcf",
        type=Path,
        help="Path to the output reformatted VCF file.",
    )
    args = parser.parse_args()

    vcf = VcfFile(args.input_vcf)
    vcf.parse_columns(args.caller)
    vcf.write_vcf(args.output_vcf)
