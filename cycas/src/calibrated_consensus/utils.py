import re
from typing import List, Optional

import numpy as np


def elongate_cigar(short_cigar: str) -> str:
    """Converts a short cigar with format 9M1I10X16M to a long string of
    repeated MMMMMMMMMIXXXXXXXXXXMMMMMMMMMMMMMMMM

    Args:
        short_cigar: cigar string to be elongated
    """
    cigar_counts = re.split("M|I|D|N|S|H|P|=|X", short_cigar)
    cigar_strs = re.split("[0-9]", short_cigar)

    cigar_counts = [c for c in cigar_counts if c != ""]
    cigar_strs = [s for s in cigar_strs if s != ""]

    assert len(cigar_strs) == len(cigar_counts)

    longcigar = ""
    for c, s in zip(cigar_counts, cigar_strs):
        longcigar += s * int(c)
    return longcigar


def make_align_array(
    cigar: str,
    reference_seq: str,
    query_seq: str,
    query_qual: Optional[List[str]] = None,
) -> np.ndarray:
    """Create an alignment array between a reference and a query sequence based on
    the cigar string.

    Reference and query strings are expected to be trimmed for alignment clipping.

    Args:
        cigar: a string in long format (e.g. MMMMMXMMMIMMM) representing the alignment.
            Expects MIDX. If mismatches are coded as M also, will be converted to X.
        reference_seq: the reference sequence to align against.
        query_seq: the query sequence to align.
        query_qual: optional list of quality scores for the query sequence, will be
            converted to strings.

    Returns:
        a [4, len(cigar)] numpy array where the rows are:
            - reference sequence
            - cigar
            - query sequence
            - query quality
        Gaps are indicated as '-' in both the reference and query sequences.
    """

    aln_array = np.empty((4, len(cigar)), dtype=np.dtypes.StrDType)

    if reference_seq != query_seq:
        c_i = 0
        r_i = 0
        for i, c in enumerate(cigar):
            if c == "M":
                if reference_seq[r_i] != query_seq[c_i]:
                    c = "X"
                aln_array[0, i] = reference_seq[r_i]
                aln_array[1, i] = c
                aln_array[2, i] = query_seq[c_i]
                if query_qual:
                    aln_array[3, i] = query_qual[c_i]
                c_i += 1
                r_i += 1
            elif c == "I":
                aln_array[0, i] = "-"
                aln_array[1, i] = c
                aln_array[2, i] = query_seq[c_i]
                if query_qual:
                    aln_array[3, i] = query_qual[c_i]
                c_i += 1
            elif c == "D":
                aln_array[0, i] = reference_seq[r_i]
                aln_array[1, i] = c
                aln_array[2, i] = "-"
                if query_qual:
                    aln_array[3, i] = ""
                r_i += 1
            else:
                raise ValueError(f"Unknown cigar operation: {c}")

    else:
        aln_array[0, :] = list(reference_seq)
        aln_array[1, :] = list(cigar)
        aln_array[2, :] = list(query_seq)
        if query_qual:
            aln_array[3, :] = list(query_qual)

    return aln_array


def phred_to_string(scores: List[int]) -> str:
    """
    Convert a list of integer Phred quality scores into a FASTQ quality string.

    Args:
        scores: a list of integer Phred quality scores

    Returns:
        a string representing the FASTQ quality scores
    """
    return "".join(chr(s + 33) for s in scores)


def string_to_phred(qual_str: str) -> List[int]:
    """
    Convert a FASTQ quality string into a list of integer Phred quality scores.

    Args:
        qual_str: a FASTQ quality string

    Returns:
        a list of integer Phred quality scores
    """
    return [ord(c) - 33 for c in qual_str]
