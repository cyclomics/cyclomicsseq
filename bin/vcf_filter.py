#!/usr/bin/env python

from abc import abstractmethod
from pathlib import Path
import pandas as pd
import io


class VCF_file:
    def __init__(self, vcf_file):
        self.vcf_file = vcf_file
        self.vcf_header = ""
        self.vcf = self.read_vcf(self.vcf_file)

    @staticmethod
    def relaxed_float(x):
        try:
            float(x)
        except ValueError:
            return float(0)
            
    def read_vcf(self, path):
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
                    lambda x: relaxed_float((x.split(":")[i] if (x.split(":")[i]) else 0))
                )
        return df

    def write(self, path):
        with open(path, "w") as new_vcf:
            new_vcf.writelines(self.vcf_header)

            writeable_vcf = self.vcf.rename(columns={"CHROM": "#CHROM"})
            writeable_vcf = writeable_vcf[
                ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
            ]
            # new_vcf.writelines((x.lstrip() for x in writeable_vcf.to_csv(sep='\t').split('\n')))
            new_vcf.writelines(writeable_vcf.to_csv(sep="\t", index=False))

    def filter(
        self,
        min_dir_ratio=0.001,
        min_dir_count=5,
        min_dqp=500,
        min_vaf=0.003,
        min_dir_ratio_ratio=0.1,
    ):
        # nothing to filter
        if self.vcf.empty:
            return

        ratio_ratios = (min_dir_ratio_ratio, 1 / min_dir_ratio_ratio)

        print("pre filter")
        print(self.vcf.shape)
        self.vcf = self.vcf[self.vcf["SAME"] > 0.9]
        print("SAME filter")
        print(self.vcf.shape)
        self.vcf = self.vcf[self.vcf["DPQ"] > min_dqp]
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
        self.vcf = self.vcf[self.vcf["FWDR"] / self.vcf["REVR"] > min(ratio_ratios)]
        self.vcf = self.vcf[self.vcf["FWDR"] / self.vcf["REVR"] < max(ratio_ratios)]
        print("relative ratio filter")
        print(self.vcf.shape)

        print(self.vcf)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Filter a vcf")

    parser.add_argument("variant_vcf", type=Path)
    parser.add_argument("file_out", type=Path)
    args = parser.parse_args()

    vcf = VCF_file(args.variant_vcf)
    vcf.filter()
    vcf.write(args.file_out)

    # vcf = VCF_file('/home/dami/Downloads/vcf_filter_error/dcffb7f76acf0f1aae134b54a7b99c/FAT55666_validated.vcf')
    # print('loaded')