#!/usr/bin/env python

import argparse
import io
from pathlib import Path
from typing import Any, List

import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter a VCF")

    parser.add_argument(
        "-i",
        "--input-vcf",
        type=Path,
        required=True,
        help="Input VCF file with variant evidence to filter.",
    )
    parser.add_argument(
        "-o",
        "--output-vcf",
        type=Path,
        required=True,
        help="Output VCF file to which passed variants will be written.",
    )
    parser.add_argument(
        "-p",
        "--perbase-table",
        type=Path,
        required=True,
        help="Input Perbase TSV table with positional depths.",
    )
    parser.add_argument(
        "--min-dir-ratio",
        type=float,
        required=False,
        default=0.001,
        help="Minimum ratio of variant-supporting reads in each direction (default: 0.001).",
    )
    parser.add_argument(
        "--min-dir-count",
        type=float,
        required=False,
        default=50,
        help="Minimum number of variant-supporting reads in each direction (default: 50).",
    )
    parser.add_argument(
        "--min-dpq",
        type=float,
        required=False,
        default=1_000,
        help="Minimum positional depth after Q filtering (default: 1_000).",
    )
    parser.add_argument(
        "--min-dpq-n",
        type=int,
        required=False,
        default=4,
        help="Number of flanking nucleotides to the each position that will determine the window size for local maxima calculation (default = 4).",
    )
    parser.add_argument(
        "--min-dpq-ratio",
        type=float,
        required=False,
        default=0.3,
        help="Ratio of local depth maxima that will determine the minimum depth at each position (default = 0.3).",
    )
    parser.add_argument(
        "--min-vaf",
        type=float,
        required=False,
        default=0.002,
        help="Minimum variant allele frequency (default: 0.002).",
    )
    parser.add_argument(
        "--min-rel-ratio",
        type=float,
        required=False,
        default=0.3,
        help="Minimum relative ratio between forward and reverse variant-supporting reads (default: 0.3).",
    )
    parser.add_argument(
        "--min-abq",
        type=float,
        required=False,
        default=70,
        help="Minimum average base quality (default: 70).",
    )

    return parser.parse_args()


def get_depth_table(perbase_tsv: Path) -> pd.DataFrame:
    df = pd.read_table(perbase_tsv, sep="\t")
    return df[["POS", "DEPTH"]]


class VCF_file:
    def __init__(self, vcf_file):
        self.vcf_file = vcf_file
        self.vcf_header = ""
        self.vcf = self.read_vcf(self.vcf_file)

    @staticmethod
    def relaxed_float(x: Any) -> float:
        """
        Try except wrapper with value error catch for non floatable objects
        """
        try:
            my_float = float(x)
        except ValueError:
            my_float = float(0)
        return my_float

    def read_vcf(self, path: Path) -> pd.DataFrame:
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
            },
            sep="\t",
        ).rename(columns={"#CHROM": "CHROM"})

        if not df.empty:
            formats = df.FORMAT[0].split(":")
            for i, fmt in enumerate(formats):
                df[fmt] = df.Sample1.apply(
                    lambda x: self.relaxed_float(x.split(":")[i])
                    if (x.split(":")[i])
                    else 0
                )

        return df

    def write(self, path):
        with open(path, "w") as new_vcf:
            new_vcf.writelines(self.vcf_header)

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
                    "Sample1",
                ]
            ]
            # new_vcf.writelines((x.lstrip() for x in writeable_vcf.to_csv(sep='\t').split('\n')))
            new_vcf.writelines(writeable_vcf.to_csv(sep="\t", index=False))

    def apply_min_depth(
        self, depth_table: pd.DataFrame, n: int = 4, ratio: float = 0.3
    ) -> List[bool]:
        """
        Apply a minimum depth filter based on local maxima.

        Args:
            n: Number of positions to the left and to the right, which will
                be used to calculate local depth maxima (default = 4).
            ratio: Ratio of the local depth maxima that will determine the
                minimum depth at each position (default = 0.3).

        Returns:
            List of booleans, with size equal to the number of VCF entries,
                flagging whether a VCF entry passes or fails this filter.
        """

        min_depth = []
        for row in self.vcf.iterrows():
            pos = row[1]["POS"]
            pos_depth = int(depth_table[depth_table["POS"] == pos]["DEPTH"])

            # Calculate local maximum
            depth_range = depth_table.loc[
                depth_table["POS"].between(pos - n, pos + n), "DEPTH"
            ]
            local_max = max(depth_range) * ratio

            min_depth.append(pos_depth >= local_max)

        return min_depth

    def filter(
        self,
        depth_table: pd.DataFrame,
        min_dir_ratio: float = 0.001,
        min_dir_count: int = 5,
        min_dpq: int = 1_000,
        min_dpq_n: int = 4,
        min_dpq_ratio: float = 0.3,
        min_vaf: float = 0.002,
        min_rel_ratio: float = 0.3,
        min_abq: int = 70,
    ):
        # nothing to filter
        if self.vcf.empty:
            return

        result = []
        for row in self.vcf.iterrows():
            r = row[1]["REVR"]
            f = row[1]["FWDR"]
            low = min(f, r)
            high = max(f, r)

            if high == 0:
                result.append(0)
            else:
                result.append(low / high)
        self.vcf["RELR"] = result

        print("pre filter")
        print(self.vcf.shape)
        self.vcf = self.vcf[self.vcf["SAME"] > 0.9]
        print("SAME filter")
        self.vcf = self.vcf[self.vcf["ABQ"] > min_abq]
        print("ABQ filter")
        print(self.vcf.shape)
        # self.vcf = self.vcf[self.vcf["DPQ"] > min_dpq]
        self.vcf = self.vcf[vcf.apply_min_depth(depth_table, min_dpq_n, min_dpq_ratio)]
        print("DPQ filter")
        print(self.vcf.shape)
        self.vcf = self.vcf[self.vcf["VAF"] > min_vaf]
        print("VAF filter")
        print(self.vcf.shape)
        self.vcf = self.vcf[self.vcf["REVC"] > min_dir_count]
        print("REVC filter")
        print(self.vcf.shape)
        self.vcf = self.vcf[self.vcf["FWDC"] > min_dir_count]
        print("FWDC filter")
        print(self.vcf.shape)
        self.vcf = self.vcf[self.vcf["FWDR"] > min_dir_ratio]
        print("FWDR filter")
        print(self.vcf.shape)
        self.vcf = self.vcf[self.vcf["REVR"] > min_dir_ratio]
        print("REVR filter")
        print(self.vcf.shape)
        self.vcf = self.vcf[self.vcf["RELR"] > min_rel_ratio]
        print("relative ratio filter")
        print(self.vcf.shape)

        print(self.vcf)


if __name__ == "__main__":
    args = parse_arguments()

    depth_table = get_depth_table(args.perbase_table)

    vcf = VCF_file(args.input_vcf)
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

    # vcf = VCF_file('/home/dami/Software/cycloseq/tmp.vcf')
    # vcf.filter()
    # print('loaded')
