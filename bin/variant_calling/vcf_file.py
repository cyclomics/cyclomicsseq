#!/usr/bin/env python

from functools import partial
from pathlib import Path
from typing import Any

import pandas as pd


class Vcf:
    """
    VCF file object where filters will be applied to.
    """

    def __init__(self, input_file):
        self.input_file = input_file
        self.header = ""
        self.table = self.read(self.input_file)

    @staticmethod
    def _relaxed_float(x: Any) -> float:
        """
        Try except wrapper with value error catch for non floatable objects.
        """
        try:
            my_float = float(x)
        except ValueError:
            my_float = float(0)
        return my_float

    @staticmethod
    def _split_info_annot(cell: str) -> dict:
        return dict(x.split("=") for x in cell.split(";"))

    @staticmethod
    def _join_info_annot(cell: dict) -> str:
        annotations = []
        for key, value in cell.items():
            annotations.append("=".join((key, value)))

        annotations = ";".join(annotations)
        return annotations

    def read(self, filepath: Path) -> pd.DataFrame:
        """
        Read in a VCF file and a return it as a pandas DataFrame.
        """
        # Store the VCF header so we can write it out later
        with open(filepath, "r") as file:
            header = []
            for n, line in enumerate(file):
                if line.startswith("##"):
                    header.append(line)
                else:
                    break

        self.header = header

        # Read VCF into dataframe, skip n commen
        df = pd.read_csv(
            filepath,
            skiprows=n,
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
                "sample1": str,
            },
            sep="\t",
        ).rename(columns={"#CHROM": "CHROM"})

        df.INFO = df.INFO.apply(func=self._split_info_annot)

        return df

    def write(self, path):
        """
        Write filtered VCF_file object as a new VCF file to path.
        """
        self.table.INFO = self.table.INFO.apply(func=self._join_info_annot)
        with open(path, "w") as file:
            file.writelines(self.header)
            self.table = self.table.rename(columns={"CHROM": "#CHROM"})
            self.table.columns = self.table.columns.str.upper()
            file.writelines(self.table.to_csv(sep="\t", index=False))

    @staticmethod
    def _apply_min_ao(row, ao_filter: int) -> bool:
        result = int(row["INFO"]["AO"]) > ao_filter

        return result

    @staticmethod
    def _apply_min_abq(row, abq_filter: int) -> bool:
        result = float(row["INFO"]["QA"]) / int(row["INFO"]["AO"]) > abq_filter

        return result

    @staticmethod
    def _apply_min_depth(row, depth_filter: int) -> bool:
        result = int(row["INFO"]["DP"]) > depth_filter

        return result

    @staticmethod
    def _apply_local_depth(
        row,
        depth_df: pd.DataFrame,
        n: int,
        ratio: float,
    ) -> bool:
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
        if depth_df is None:
            # What if there is simply a problem with depth_df?
            # Should it instead return True?
            return False

        chrom = row["CHROM"]
        pos = row["POS"]
        chr_depth = depth_df[depth_df["REF"] == chrom]
        #  Pos could be missing in the depth table, have backup for snp
        if pos in chr_depth["POS"]:
            pos_depth = int(chr_depth[chr_depth["POS"] == pos]["DEPTH"])
        else:
            try:
                pos_depth = int(row["INFO"]["DP"])
            except KeyError:
                pos_depth = 0

        # Calculate local maximum
        depth_range = chr_depth.loc[chr_depth["POS"].between(pos - n, pos + n), "DEPTH"]
        local_max = max(depth_range) * ratio

        result = pos_depth >= local_max

        return result

    @staticmethod
    def _apply_dynamic_min_vaf(row, params) -> bool:
        type = row["INFO"]["TYPE"]
        candidate_depth = int(row["INFO"]["DP"])
        candidate_vaf = int(row["INFO"]["AO"]) / candidate_depth

        if type == "snp":
            type_params = params["snp"]
        else:
            type_params = params["other"]

        result = False
        for param_depth, param_vaf in sorted(type_params.items(), reverse=True):
            if candidate_depth > param_depth and candidate_vaf > param_vaf:
                result = True
                break

        return result

    @staticmethod
    def _apply_max_sap(row, sap_filter: int) -> bool:
        result = float(row["INFO"]["SAP"]) < sap_filter

        return result

    @staticmethod
    def _apply_min_ratio(row, ratio_filter: float) -> bool:
        saf = int(row["INFO"]["SAF"])
        sar = int(row["INFO"]["SAR"])
        low = min(saf, sar)
        high = max(saf, sar)
        ratio = low / high if high > 0 else 0
        result = ratio > ratio_filter

        return result

    def filter(
        self,
        depth_table: pd.DataFrame,
        dynamic_vaf_params: dict,
        min_ao: int = 10,
        min_dpq: int = 5_000,
        min_dpq_n: int = 25,
        min_dpq_ratio: float = 0.3,
        max_sap: int = 60,
        min_rel_ratio: float = 0.3,
        min_abq: int = 70,
    ):
        """
        Filters the VCF DataFrame based on parsed or default parameters.

        Args:
            depth_table: Perbase depth table as a pandas DataFrame, used with
                min_dpq_n and min_dpq_ratio to calculate local depth maxima and
                apply a minimum depth filter based on local maxima.
            dynamic_vaf_params: Input file with params for Dynamic VAF filtering.
            min_ao: Minimum number of variant-supporting reads.
            min_dpq: Minimum positional depth after Q filtering.
            min_dpq_n: Number of flanking nucleotides to the each position that will
                determine the window size for local maxima calculation.
            min_dpq_ratio: Ratio of local depth maxima that will determine the
                minimum depth at each position.
            max_sap: Maximum Phred-scaled strand balance probability before
                alternative allele is considered strand-imbalanced.
            min_rel_ratio: Minimum relative ratio between forward and reverse
                variant-supporting reads.
            min_abq: Minimum average base quality.
        """
        # nothing to filter
        if self.table.empty:
            return

        # Define filtering partial functions
        min_ao_filter = partial(self._apply_min_ao, ao_filter=min_ao)
        min_abq_filter = partial(self._apply_min_abq, abq_filter=min_abq)
        min_depth_filter = partial(self._apply_min_depth, depth_filter=min_dpq)
        local_depth_filter = partial(
            self._apply_local_depth,
            depth_df=depth_table,
            n=min_dpq_n,
            ratio=min_dpq_ratio,
        )
        dynamic_vaf_filter = partial(
            self._apply_dynamic_min_vaf,
            params=dynamic_vaf_params,
        )
        max_sap_filter = partial(self._apply_max_sap, sap_filter=max_sap)
        min_ratio_filter = partial(self._apply_min_ratio, ratio_filter=min_rel_ratio)

        # List of filters to be applied
        filters = [
            min_ao_filter,
            min_abq_filter,
            min_depth_filter,
            local_depth_filter,
            dynamic_vaf_filter,
            min_ratio_filter,
        ]
        if max_sap > 0:
            filters.append(max_sap_filter)

        # Loop over list of filters and apply them in succession
        for filter_function in filters:
            self.table = self.table[self.table.apply(func=filter_function, axis=1)]
            print(f"After {filter_function.func.__name__}: {self.table.shape}")

        self.table["FILTER"] = "PASS"
