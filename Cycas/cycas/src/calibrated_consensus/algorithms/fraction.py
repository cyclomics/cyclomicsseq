import os
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import h5py
import joblib
import mappy
import numpy as np
import pandas as pd
import pysam
from loguru import logger
from sklearn.isotonic import IsotonicRegression

from ..utils import elongate_cigar, make_align_array
from .abstract import AbstractCalibratedConsensus


class FractionCalibratedConsensus(AbstractCalibratedConsensus):
    """
    Fraction-based calibrated consensus model.

    Details on the consensus algorithm:



    Details on the calibration algorithm:


    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.model = None

    def prep_file(
        self,
        hdf5_file: Path,
        output_dir: Path,
        reference_file: Path,
        min_alignment_quality: int,
        vcf_file: Optional[Path] = None,
        skip_variants: bool = True,
    ):
        if vcf_file is not None:
            if not skip_variants:
                raise NotImplementedError("Not skipping variants is not implemented")

        logger.info("Loading reference genome")
        aligner = mappy.Aligner(str(reference_file), preset="short")
        if vcf_file is not None:
            logger.info("Loading VCF file")
            vcf = pysam.VariantFile(vcf_file, "r")
        else:
            vcf = None

        output_file_metadata = os.path.join(
            output_dir, f"{hdf5_file.stem}.calib_metadata.csv"
        )
        output_file_data = os.path.join(output_dir, f"{hdf5_file.stem}.calib_data.csv")

        calib_dict = dict()
        rows = list()

        logger.info("File processing for calibration begins")
        with h5py.File(hdf5_file, "r") as h5_file:
            for sample_id in h5_file.keys():
                # should be sample_id, data, metadata
                data, metadata = self._extract_data(h5_file[sample_id])
                data = {"sample_id": sample_id, "data": data, "metadata": metadata}
                data["consensus_result"] = self._make_consensus(data["data"])
                evaluation_result = self.evaluate_consensus(
                    data, aligner, min_alignment_quality, vcf, skip_variants
                )
                row = evaluation_result["row"]
                calib = evaluation_result["calibration"]

                rows.append(row)
                for k, v in calib.items():
                    if k not in calib_dict:
                        calib_dict[k] = v
                    else:
                        calib_dict[k] += v

        logger.info("Saving results")
        pd.DataFrame(rows).to_csv(output_file_metadata, header=True, index=False)
        pd.DataFrame(calib_dict).transpose().to_csv(
            output_file_data, index=True, header=False
        )

        if vcf is not None:
            vcf.close()

    def evaluate_consensus(
        self, data, aligner, min_alignment_quality, vcf, skip_variants
    ) -> Dict[str, str]:
        # main result data structure
        evaluation_result = {
            "row": {
                "read_id": data["sample_id"],
                "consensus_len": len(data["consensus_result"]["consensus_sequence"]),
                "max_depth": data["consensus_result"]["consensus_max_depth"],
                "is_aligned": False,
                "mapping": "",
                "clip": "",
                "mapped_len": 0,
                "mapq": 0,
                "contains_variants": False,
                "distance": 0,
                "structure_classification": data["metadata"][
                    "structure_classification"
                ],
                "calibration_group": self._grouper(
                    data["metadata"]["structure_classification"]
                ),
                "five_prime_4mer": "",
                "three_prime_4mer": "",
            },
            "calibration": dict(),
        }

        consensus_seq = data["consensus_result"]["consensus_sequence"]
        consensus_length = len(consensus_seq)
        consensus_qual = data["consensus_result"]["consensus_quality_uncalibrated"]

        if consensus_length >= 8:
            evaluation_result["row"]["five_prime_4mer"] = consensus_seq[:4]
            evaluation_result["row"]["three_prime_4mer"] = consensus_seq[-4:]

        ##### ALIGNMENT ################################################################
        # first check if the consensus sequence aligns. If it does, record the result
        # and continue. If it does not align, early return.

        # align the consensus sequence against the reference
        aln = None
        for aln in aligner.map(consensus_seq):
            break

        # if we do not get an alignment, not much to be done here, early return
        if aln is None:
            return evaluation_result
        # if we have an alignment, then get the corresponding reference sequence
        else:
            ref_seq = aligner.seq(aln.ctg, aln.r_st, aln.r_en)
            if aln.strand == -1:
                ref_seq = mappy.revcomp(ref_seq)

        # record some data on the alignment
        strand = "F" if aln.strand == 1 else "R"
        evaluation_result["row"]["is_aligned"] = True
        evaluation_result["row"][
            "mapping"
        ] = f"{aln.ctg}:{aln.r_st}-{aln.r_en}:{strand}"
        evaluation_result["row"]["clip"] = f"{aln.q_st}-{aln.q_en}"
        evaluation_result["row"]["mapped_len"] = aln.q_en - aln.q_st
        evaluation_result["row"]["mapq"] = aln.mapq

        if aln.mapq < min_alignment_quality:
            return evaluation_result

        ##### VARIANTS #################################################################
        # Check if, given the alignment, we have an overlap with a variant and
        # proceed accordingly depending if we want to skip them or not.

        # check if there is a variant present if we have a VCF file
        if vcf is not None:
            rec = None
            for rec in vcf.fetch(aln.ctg, aln.r_st, aln.r_en):
                break
            if rec is not None:
                evaluation_result["row"]["contains_variants"] = True
        else:
            evaluation_result["row"]["contains_variants"] = None

        # if we skip variants, then early return
        if evaluation_result["row"]["contains_variants"] and skip_variants:
            return evaluation_result

        # TODO: so far we skip all variants, but we probably want to also evaluation
        # variant regions, in principle it should be fairly simple to modify the
        # reference sequence and include the variants. At the same time, it can be a
        # bit complicated, because what if its not an homozygous variant, then you can
        # still get a read that is correct. So probably best would be to return as a
        # result the VAF also to later filter non 100% VAF variants from the analysis.

        ##### ALIGNMENT_ARRAY ##########################################################
        # make an alignment array between the consensus and the reference for
        # evaluation of each base
        consensus_seq = consensus_seq[aln.q_st : aln.q_en]
        consensus_qual = consensus_qual[aln.q_st : aln.q_en]

        aln_array = make_align_array(
            cigar=elongate_cigar(aln.cigar_str),
            reference_seq=ref_seq,
            query_seq=consensus_seq,
            query_qual=consensus_qual,
        )
        # add the edit distance between the two sequences
        evaluation_result["row"]["distance"] = (aln_array[1] != "M").sum()

        if evaluation_result["row"]["distance"] >= 20:
            return evaluation_result

        # evaluate each base
        group_index = evaluation_result["row"]["calibration_group"]
        ref_base_i = aln.q_st

        # for each column in the alignment array
        for col in range(aln_array.shape[1]):
            # get the event information
            cigar_op = aln_array[1, col]
            consensus_base = aln_array[2, col]
            qual = aln_array[3, col]
            # if there is no quality, then it was a deletion, so record how many repeats
            # were in this read to calculate the chance of a deletion given a certain
            # amount of repeats
            if cigar_op == "D":
                qual = f"{evaluation_result['row']['max_depth']}/{evaluation_result['row']['max_depth']}"

            # keep track of the relative position of the base in terms of quantile
            rel_pos = ref_base_i / consensus_length
            # put bin 19 and 20 together cause otherwise the last bin is just the last
            # base
            rel_pos_bin = int(rel_pos * 20)
            if rel_pos_bin == 20:
                rel_pos_bin = 19

            # use a key for fast lookup in the future table. We index based on:
            # - the base ACGTN-
            # - the uncalibrated quality
            # - the group index, so the type of read 1D, 2D, 3D
            # - the relative position of the base in the read
            key = f"{consensus_base}_{qual}_{group_index}_{rel_pos_bin}"

            # if the key is still not in the dict, initialize an array of size 4
            # to keep track of matches, deletions, insertions and mismatches
            if key not in evaluation_result["calibration"]:
                evaluation_result["calibration"][key] = np.zeros((4,), dtype=int)

            # store the operation in the corret position
            if cigar_op == "M":
                evaluation_result["calibration"][key][0] += 1
                ref_base_i += 1
            elif cigar_op == "D":
                evaluation_result["calibration"][key][1] += 1
            elif cigar_op == "I":
                evaluation_result["calibration"][key][2] += 1
                ref_base_i += 1
            elif cigar_op == "X":
                evaluation_result["calibration"][key][3] += 1
                ref_base_i += 1
            else:
                raise ValueError(f"Unknown cigar operation: {cigar_op}")

        return evaluation_result

    def _make_consensus(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Wrapper around make_consensus to take care of compliant input and output"""

        msa = data["msa_seq"].astype(str)

        consensus_seq, consensus_qual = self.make_consensus(msa)

        data["consensus_sequence"] = consensus_seq
        data["consensus_quality_uncalibrated"] = consensus_qual
        data["consensus_max_depth"] = msa.shape[0]

        return data

    @staticmethod
    def make_consensus(msa: np.ndarray) -> Tuple[str, List[str]]:
        """Consensus algorithm that only considers a consensus base if it is present in
        >50% of the subreads. It also reports the fraction of subreads that support the
        consensus base as the uncalibrated quality score.

        Args:
            msa: a 2D numpy array where rows are each repeat and columns each base
                position. It is expected to be an array of strings. Importantly, '-'
                indicates gaps, and '_' indicate padding. Padding does not count towards
                majority voting.

        Returns:
            - a string with the consensus sequence
            - a list of strings with the fraction of supported bases as
                majority_base/total_bases_at_column
        """

        m = {"A": 0, "C": 1, "G": 2, "T": 3, "-": 4}

        consensus_seq = list()
        consensus_qual = list()
        # iterate over each column of the MSA
        for pos in range(msa.shape[1]):
            c = np.zeros(5, dtype=int)  # counter of each valid base
            for b in msa[:, pos]:
                # skip padding
                if b == "_":
                    continue
                # add counter to each base
                c[m[b]] += 1

            mc = np.max(c)  # max count of any base
            t = np.sum(c)  # total valid bases

            # if one base is more than half, then its the best base
            if mc > t / 2:
                best_base = "ACGT-"[np.argmax(c)]
                final_qual = f"{mc}/{t}"
            # otherwise report N
            else:
                best_base = "N"
                final_qual = f"{mc}/{t}"

            # if the best base is not a gap then we add it to the consensus
            if best_base != "-":
                consensus_seq.append(best_base)
                consensus_qual.append(final_qual)

        consensus_seq = "".join(consensus_seq)

        return consensus_seq, consensus_qual

    def fit_model(self, input_dir: Path, max_depth: int, min_coverage: int):
        logger.info("Reading calibration data")
        dfs = list()
        for file in os.listdir(input_dir):
            if not file.endswith(".calib_data.csv"):
                continue

            calib_data_file = os.path.join(input_dir, file)
            logger.debug(f"Reading: {calib_data_file}")
            df = pd.read_csv(calib_data_file, index_col=0, header=None)
            dfs.append(df)

        combined = pd.concat(dfs, ignore_index=False)

        # Group by the key and sum the rest of the columns
        result = combined.groupby(combined.index, as_index=True).sum()
        df = result.reset_index()  # 'index' becomes a new column
        df[["base", "qual", "type", "pos"]] = df[0].str.split("_", expand=True)
        df["pos"] = df["pos"].astype(int)

        df = df.drop([0], axis=1)
        df = df.rename({1: "M", 2: "D", 3: "I", 4: "X"}, axis=1)

        self._model = {}

        logger.info("Begin model fitting")
        read_types = np.unique(df["type"])
        for read_type in read_types:
            logger.info(f"Fitting model for read type: {read_type}")
            df_type = df[df["type"] == read_type].copy()
            read_type_models = self._fit_model_type(
                df=df_type,
                max_depth=max_depth,
                min_coverage=min_coverage,
            )
            self._model[read_type] = read_type_models

        self.model = deepcopy(self._model)
        del self._model
        logger.info("Model fitting finished")

    def _fit_model_type(
        self, df: pd.DataFrame, max_depth: int, min_coverage: int
    ) -> Dict[str, Any]:
        models = {
            "zero_discordance": dict(),
            "nonzero_discordance": dict(),
        }
        relative_positions = np.unique(df["pos"])

        # fit 0 discordance model
        data = df[df["base"].isin(["A", "C", "G", "T"])].copy()
        data["total"] = data["M"] + data["I"] + data["X"]
        data["error"] = data["I"] + data["X"]
        data["correct"] = data["M"]
        data[["concordance", "depth"]] = data["qual"].str.split("/", expand=True)
        data["concordance"] = data["concordance"].astype(int)
        data["depth"] = data["depth"].astype(int)

        data = (
            data.groupby(["depth", "concordance", "pos"])[["total", "error"]]
            .sum()
            .reset_index()
        )
        data["discordance"] = np.round(data["concordance"] / data["depth"], 2)
        data = data[data["discordance"] == 1.0]
        data = data[data["total"] >= min_coverage]

        data["Q"] = np.round(
            (10 * np.log10(1 / ((data["error"] + 1) / (data["total"] + 1)))), 0
        ).astype(int)
        data = data[data["depth"] <= max_depth]
        data = data.sort_values(["depth"], ascending=True)

        for p in relative_positions:
            x = np.array(data.loc[data["pos"] == p, "depth"])
            y = np.array(data.loc[data["pos"] == p, "Q"])

            iso = IsotonicRegression(increasing=True)
            iso.fit(x, y)
            iso.set_params(**{"y_min": np.min(x).item(), "y_max": np.max(x).item()})
            models["zero_discordance"][int(p)] = iso

        # fit >= 1 discordance
        data = df[df["base"].isin(["A", "C", "G", "T"])].copy()
        data["total"] = data["M"] + data["I"] + data["X"]
        data["error"] = data["I"] + data["X"]
        data["correct"] = data["M"]
        data[["concordance", "depth"]] = data["qual"].str.split("/", expand=True)
        data["concordance"] = data["concordance"].astype(int)
        data["depth"] = data["depth"].astype(int)

        data = (
            data.groupby(["depth", "concordance", "pos"])[["total", "error"]]
            .sum()
            .reset_index()
        )
        data["discordance"] = np.round(data["concordance"] / data["depth"], 2)
        data = data[data["discordance"] < 1.0]

        data = data[data["total"] >= 100]
        data = data[data["depth"] <= 50]

        data = (
            data.groupby(["discordance", "pos"])[["total", "error"]].sum().reset_index()
        )

        data["Q"] = np.round(
            (10 * np.log10(1 / ((data["error"] + 1) / (data["total"] + 1)))), 0
        ).astype(int)

        for p in relative_positions:
            x = np.array(data.loc[data["pos"] == p, "discordance"])
            y = np.array(data.loc[data["pos"] == p, "Q"])

            iso = IsotonicRegression(increasing=True)
            iso.fit(x, y)
            iso.set_params(**{"y_min": np.min(x).item(), "y_max": np.max(x).item()})
            models["nonzero_discordance"][int(p)] = iso

        return models

    def calibrate(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """Calibrate the consensus sequence quality scores based on the type of read
        and fitted model.

        Args:
            data: a dictionary with at least the following
                -   'consensus_quality_uncalibrated': a list of strings with the
                    uncalibrated quality scores in the format majority/total
                -   'metadata': a dictionary with at least the key 'structure_classification'
                    which indicates the type of read, e.g. 1D, 2D, 3D

        Returns:
            the input dictionary with one additional key:
                -   'consensus_quality_calibrated': a list of integers with the calibrated
                    quality scores

        Raises:
            ValueError: if the model is not loaded
        """

        if self.model is None:
            error_msg = "Model is None, did you forget to run `fit_model` or `load`?"
            logger.error(error_msg)
            raise ValueError(error_msg)

        consensus_qual = data["consensus_result"]["consensus_quality_uncalibrated"]
        read_type = self._grouper(data["metadata"]["structure_classification"])

        calibrated_qual = list()

        rel_pos = np.arange(len(consensus_qual)) / len(consensus_qual)
        rel_pos = (rel_pos * 20).astype(int)
        rel_pos[rel_pos == 20] = 19

        for i in range(len(consensus_qual)):
            cd = consensus_qual[i].split("/")
            concordance = int(cd[0])
            depth = int(cd[1])

            rp = rel_pos[i]

            if concordance == depth:
                val = depth
                model_str = "zero_discordance"
            else:
                val = concordance / depth
                model_str = "nonzero_discordance"

            iso = self.model[read_type][model_str][rp]
            params = iso.get_params()
            val = np.clip(val, params["y_min"], params["y_max"])
            q = np.round(iso.predict([val])[0]).astype(int).item()
            calibrated_qual.append(q)

        data["consensus_result"]["consensus_quality_calibrated"] = calibrated_qual
        return data

    def save(self, path: Path):
        logger.info(f"Saving models to: {path}")
        joblib.dump(self.model, os.path.join(path, "calibration_models.pkl"))
        logger.info("Saving succesful")

    def load(self, path: Path):
        logger.info(f"Loading models from: {path}")
        self.model = joblib.load(os.path.join(path, "calibration_models.pkl"))
        logger.info("Loading succesful")

    def _extract_data(self, sample_h5):
        """Processes a single sample from the hdf5 file

        Args:
            sample_h5: an h5py file object

        Returns:
            A tuple with the processed data.
            The first element is the sample id.
            The second element being a dictionary with the data.
            The third element being a dictionary with the attributes
        """

        direction_mapper = {b"0": "R", b"1": "F"}

        subreads = list(sample_h5.keys())
        x_len = len(sample_h5[subreads[0]]["nuc_array"][:])
        x = np.zeros((len(subreads), x_len), dtype="S1")
        direction = np.zeros((len(subreads),), dtype="S1")

        for i, subread in enumerate(sample_h5.keys()):
            x[i] = sample_h5[subread]["nuc_array"][:]
            direction[i] = direction_mapper[sample_h5[subread]["aln_data"][2]]

        data = {
            "msa_seq": x,
            "d": direction,
        }

        metadata = dict()
        for key in sample_h5.attrs.keys():
            metadata[key] = sample_h5.attrs[key]

        return data, metadata

    def _grouper(self, structure_classification: str) -> str:
        if "1D" in structure_classification:
            return "1D"
        elif "2D" in structure_classification:
            return "2D"
        elif "3D" in structure_classification:
            return "3D"
        else:
            return "unknown"
