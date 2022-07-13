#!/usr/bin/env python


import math
from collections import Counter
from typing import List, Tuple

import pandas as pd
import numpy as np
import re

from bokeh.io import save, output_file
from bokeh.plotting import figure, show
from bokeh.layouts import row, column

chromosomal_region = Tuple[str, int, int]


def get_roi_pileup_df(
    df: pd.DataFrame, distance: int = 100
) -> List[chromosomal_region]:
    print("get_roi_pileup_df)")
    """
    Given a dataframe with CHROM and POS columns, find all streches in chromosomes that are less than `distance` appart.
    reports 

    """

    def get_continous_strech(array: np.array, distance: int):
        print("get_continous_strech)")
        """
        https://stackoverflow.com/questions/47183828/pandas-how-to-find-continuous-values-in-a-series-whose-differences-are-within-a 
        """
        m = np.concatenate(([True], array[1:] > array[:-1] + distance, [True]))
        idx = np.flatnonzero(m)
        l = array.tolist()
        return [l[i:j] for i, j in zip(idx[:-1], idx[1:])]

    assemblies = df.CHROM.unique()
    roi = []
    for chrom in assemblies:
        df_tmp = df[df.CHROM == chrom]
        streches = get_continous_strech(df_tmp.POS.to_numpy(), distance)
        for strech in streches:
            roi.append((chrom, strech[0], strech[-1]))

    return roi


def acc_to_Q(acc: float, max_Q=50) -> int:
    """
    Convert an accuracy (eg 0.99) to a Q score
    """
    # perfect acc is not Q score convertable
    if acc == 1:
        return max_Q
    else:
        return -10 * math.log10(1 - acc)


def plot_roi_table(roi):
    print("plot_roi_table)")
    pass


def _parse_pileup(row):
    print("_parse_pileup)")
    """
    . (dot) means a base that matched the reference on the forward strand
    , (comma) means a base that matched the reference on the reverse strand
    </> (less-/greater-than sign) denotes a reference skip. This occurs, for example, if a base in the reference genome is intronic and a read maps to two flanking exons. If quality scores are given in a sixth column, they refer to the quality of the read and not the specific base.
    AGTCN (upper case) denotes a base that did not match the reference on the forward strand
    agtcn (lower case) denotes a base that did not match the reference on the reverse strand
    A sequence matching the regular expression \+[0-9]+[ACGTNacgtn]+ denotes an insertion of one or more bases starting from the next position. For example, +2AG means insertion of AG in the forward strand
    A sequence matching the regular expression \-[0-9]+[ACGTNacgtn]+ denotes a deletion of one or more bases starting from the next position. For example, -2ct means deletion of CT in the reverse strand
    ^ (caret) marks the start of a read segment and the ASCII of the character following `^' minus 33 gives the mapping quality
    $ (dollar) marks the end of a read segment
    * (asterisk) is a placeholder for a deleted base in a multiple basepair deletion that was mentioned in a previous line by the -[0-9]+[ACGTNacgtn]+ notation
    """
    pileup = row.PILEUP[:]
    starts = []
    while True:
        try:
            i = pileup.index("^")
            starts.append(pileup[i : i + 3])
            pileup = pileup[:i] + pileup[i + 3 :]
        except ValueError:
            break

    ends = []
    while True:
        try:
            i = pileup.index("$")
            ends.append(pileup[i - 1 : i + 1])
            pileup = pileup[: i - 1] + pileup[i + 1 :]
        except ValueError:
            break

    skips = []
    while True:
        try:
            i = pileup.index(">")
            skips.append(pileup[i])
            pileup = pileup[:i] + pileup[i + 1 :]
        except ValueError:
            break
    while True:
        try:
            i = pileup.index("<")
            skips.append(pileup[i])
            pileup = pileup[:i] + pileup[i + 1 :]
        except ValueError:
            break

    forw_ins = []
    p = re.compile("\+[0-9]+[ACGTN]+")
    matches = list(p.finditer(pileup))
    spans = reversed([m.span() for m in matches])
    forw_ins += [m.group() for m in matches]
    for i, j in spans:
        pileup = pileup[:i] + pileup[j:]

    rev_ins = []
    p = re.compile("\+[0-9]+[acgtn]+")
    matches = list(p.finditer(pileup))
    spans = reversed([m.span() for m in matches])
    rev_ins += [m.group() for m in matches]
    for i, j in spans:
        pileup = pileup[:i] + pileup[j:]

    forw_del = []
    p = re.compile("\-[0-9]+[ACGTN]+")
    matches = list(p.finditer(pileup))
    spans = reversed([m.span() for m in matches])
    forw_del += [m.group() for m in matches]
    for i, j in spans:
        pileup = pileup[:i] + pileup[j:]

    while True:
        try:
            i = pileup.index("*")
            forw_del.append(pileup[i])
            pileup = pileup[:i] + pileup[i + 1 :]
        except ValueError:
            break

    rev_del = []
    p = re.compile("\-[0-9]+[acgtn]+")
    matches = list(p.finditer(pileup))
    spans = reversed([m.span() for m in matches])
    rev_del += [m.group() for m in matches]
    for i, j in spans:
        pileup = pileup[:i] + pileup[j:]

    while True:
        try:
            i = pileup.index("#")  # del in reverse strand
            rev_del.append(pileup[i])
            pileup = pileup[:i] + pileup[i + 1 :]
        except ValueError:
            break

    # Put back start and end bases
    for seq in starts:
        pileup += seq[-1]
    for seq in ends:
        pileup += seq[0]

    # Count matches/mismatches
    forw_matches = pileup.count(".")
    rev_matches = pileup.count(",")
    pileup = pileup.replace(".", "").replace(",", "")
    snps = Counter(pileup)

    calls = Counter(
        {
            "A": 0,
            "T": 0,
            "C": 0,
            "G": 0,
            "N": 0,
            "a": 0,
            "t": 0,
            "c": 0,
            "g": 0,
            "n": 0,
            "fINS": 0,
            "fDEL": 0,
            "rINS": 0,
            "rDEL": 0,
            "SKIP": 0,
            "FMB": 0,
            "LMB": 0,
        }
    )

    calls.update(snps)
    calls.update({"FMB": len(starts)})  # first mapping base
    calls.update({"LMB": len(ends)})  # last mapping base
    calls.update({"fDEL": len(forw_del)})
    calls.update({"fINS": len(forw_ins)})
    calls.update({"rDEL": len(forw_del)})
    calls.update({"rINS": len(forw_ins)})
    calls.update({row.REF: forw_matches})
    calls.update({row.REF.lower(): rev_matches})
    return pd.Series(calls)


def pileup_to_df(pileup_path: str) -> pd.DataFrame:
    print("pileup_to_df)")

    with open(pileup_path, "r") as pu:

        lines = pu.readlines()
        lines = [x.split("\t") for x in lines]
        df1 = pd.DataFrame(
            data=lines,
            columns=["CHROM", "POS", "REF", "COV", "PILEUP", "BASEQ", "MAPQ", "READ"],
        )
        df1["POS"] = df1.POS.astype(int)
        df1["COV"] = df1.COV.astype(int)
        df1["CONCAT"] = [len(set(v.split(","))) for v in df1.READ.values]  # concatemers
        df1 = df1.join(df1.apply(_parse_pileup, axis=1), lsuffix="", rsuffix="_right")
        df1["ACCURACY"] = df1.apply(
            compute_accuracy, insertions=True, deletions=True, axis=1
        )
        df1["Q"] = df1.ACCURACY.apply(acc_to_Q)

        return df1


def perbase_table_to_df(table_path):
    """
    https://github.com/sstadick/perbase

    """
    try:
        df = pd.read_csv(table_path, sep="\t")
    except pd.errors.EmptyDataError:
        return None
    df["REF_BASE_UPPER"] = df.REF_BASE.str.upper()
    df["REF_COUNT"] = df.lookup(df.index, df["REF_BASE_UPPER"])
    df["ACCURACY"] = df.REF_COUNT / df.DEPTH
    df["Q"] = df.ACCURACY.apply(acc_to_Q)
    df["CHROM"] = df["REF"]

    print("done")
    return df


def compute_accuracy(
    row: pd.Series, insertions=True, deletions=True, by_strand=False
) -> float:
    print("compute_accuracy)")
    """Compute base call accuracy as reference-calls/coverage"""
    # designed to be used with df.apply(compute_accuracy, axis=1)
    # insertions/deletions determine wether or not ins and dels will be counted in the COV

    # expected row keys:
    #   ['A','T','C','G','N','a','t','c','g','n','fINS','fDEL','rINS','rDEL']

    ref = row.REF
    if ref == "N":
        return np.nan

    if not by_strand:
        ref_calls = row[ref.upper()] + row[ref.lower()]
        cov = row[
            [
                "A",
                "T",
                "C",
                "G",
                "N",
                "a",
                "t",
                "c",
                "g",
                "n",
                "fINS",
                "fDEL",
                "rINS",
                "rDEL",
            ]
        ].sum()
        if not insertions:
            cov -= row["fINS"] + row["rINS"]
        if not deletions:
            cov -= row["fDEL"] + row["rDEL"]

        if cov > 0:
            return ref_calls / cov
        else:
            return 0

    else:
        # Forward
        fref_calls = row[ref.upper()]
        fcov = row[["A", "T", "C", "G", "N", "fINS", "fDEL"]].sum()
        if not insertions:
            fcov -= row["fINS"]
        if not deletions:
            fcov -= row["fDEL"]
        # Reverse
        rref_calls = row[ref.lower()]
        rcov = row[["a", "t", "c", "g", "n", "rINS", "rDEL"]].sum()
        if not insertions:
            rcov -= row["rINS"]
        if not deletions:
            rcov -= row["rDEL"]

        return pd.Series([fref_calls / fcov, rref_calls / rcov])


def plot_compare_accuracy(
    roi: List[chromosomal_region], df_list: List[pd.DataFrame], my_title="untitled"
):
    print("plot_compare_accuracy)")
    """
    Plot the q score against the position for each region of interest.

    """
    plots = []
    for region in roi:
        chrom = region[0]
        start = region[1]
        stop = region[2]

        df1 = df_list[0]
        df1 = df1[df1.CHROM == chrom]
        df1 = df1[df1.POS >= start]
        df1 = df1[df1.POS <= stop]

        df2 = df_list[1]
        df2 = df2[df2.CHROM == chrom]
        df2 = df2[df2.POS >= start]
        df2 = df2[df2.POS <= stop]

        p_freq = figure(title=f"{my_title} {chrom}:{start}-{stop}", width=1500)
        p_freq.xaxis.axis_label = "Position"
        p_freq.yaxis.axis_label = "FREQ"
        p_freq.line("POS", "Q", source=df1, line_width=2, legend_label="Raw")
        p_freq.line(
            "POS",
            "Q",
            source=df2,
            line_width=2,
            color="orange",
            legend_label="Consensus",
        )
        plots.append(p_freq)

    return column(*plots)


def make_qscore_scatter(df1, df2):
    print("make_qscore_scatter)")
    pre = df1[["CHROM", "POS", "Q"]]
    post = df2[["CHROM", "POS", "Q"]]

    df_merge = pd.merge(
        pre, post, how="left", left_on=["CHROM", "POS"], right_on=["CHROM", "POS"]
    )

    q_scatter = figure(title=f"Q score scatter, variant unaware.", width=1500)
    q_scatter.xaxis.axis_label = "ONT Q score"
    q_scatter.yaxis.axis_label = "Cyclomics Q score"

    q_scatter.scatter("Q_x", "Q_y", source=df_merge, color="orange")
    q_scatter.line([0, 50], [0, 50], color="grey", alpha=0.5)

    return q_scatter


def main(pileup_split, pileup_cons, output_plot_file):
    print("main)")
    df1 = pileup_to_df(pileup_split)
    df2 = pileup_to_df(pileup_cons)
    roi = get_roi_pileup_df(df1)
    positional_accuracy = plot_compare_accuracy(roi, [df1, df2])
    print("done")
    q_score_plot = make_qscore_scatter(df1, df2)

    output_file(output_plot_file, title="Cyclomics accuracy")

    save(column([q_score_plot, positional_accuracy]))


def new_main(perbase_path1, perbase_path2, output_plot_file):

    df1 = perbase_table_to_df(perbase_path1)
    df2 = perbase_table_to_df(perbase_path2)

    if df1 == None or df2 == None:
        f = open(output_plot_file, "w")
        f.write("<h1>One of the pileups was not deep enouhh</h1>")
        f.close()
        return
    roi = get_roi_pileup_df(df1)

    positional_accuracy = plot_compare_accuracy(roi, [df1, df2])
    q_score_plot = make_qscore_scatter(df1, df2)

    output_file(output_plot_file, title="Cyclomics accuracy")
    save(column([q_score_plot, positional_accuracy]))


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Create hist plot from a regex for fastq and fastq.gz files."
    )

    parser.add_argument("pileup_split")
    parser.add_argument("pileup_consensus")
    parser.add_argument("output")

    args = parser.parse_args()
    new_main(args.pileup_split, args.pileup_consensus, args.output)

    # pileup_split = '/media/dami/a2bc89fb-be6b-4e23-912a-0c7137cd69ad/work/9c/1793ed932d613a9e65320558e8fdf9/splibams_merged.pileup'
    # pileup_cons = '/media/dami/a2bc89fb-be6b-4e23-912a-0c7137cd69ad/work/9c/1793ed932d613a9e65320558e8fdf9/FAS37222.pileup'
    # main(pileup_split, pileup_cons, 'testQ.html')

    # df1 = perbase_table_to_df('/media/dami/a2bc89fb-be6b-4e23-912a-0c7137cd69ad/results/Cyclomics_rc/000010/perbase_test_split.tsv')
    # df2 = perbase_table_to_df('/media/dami/a2bc89fb-be6b-4e23-912a-0c7137cd69ad/results/Cyclomics_rc/000010/perbase_test.tsv')
    # plots = plot_compare_accuracy([('chr17',7579983,7580193)], [df1,df2])
    # show(plots)
