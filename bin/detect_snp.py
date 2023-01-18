#!/usr/bin/env python
from collections import Counter
from dataclasses import fields
from pathlib import Path
from typing import Tuple

import numpy as np
import pysam
from vcf_tools import VCF_entry


def extract_snp_evidence(
    pileupcolumn: pysam.PileupColumn,
    assembly: str,
    add_dels: bool = False,
    high_base_quality_cutoff: int = 80,
) -> Tuple[str, Tuple[str, str], VCF_entry]:
    """
    Find SNPs in a given pileup position.

    Args:
        pileupcolumn: Pysam pileup column at a given position.
        assembly: Reference/contig name.
        add_dels: Flag to add deletions to results (Default = False).
        high_base_quality_cutoff: Cutoff to calculate ratios on high base
            quality nucleotides only (Default = 80).

    Returns:
        A tuple with reference name, reference and alternative alleles,
        and a VCF_entry object with variant information.
    """

    vcf_entry = VCF_entry(None)
    alleles = (".", ".")

    total = 0
    counted_nucs = 0
    data_present = 0
    non_ref_ratio_filtered = 0
    fwd_count = 0
    alt_base_ratio_fwd = 0
    rev_count = 0
    alt_base_ratio_rev = 0
    tot_count = 0
    tot_ratio = 0
    obs_ratio = 0
    alt_base_mean_qual = 0
    quals_mean = 0
    hc_ratio = 0
    base = None

    # all kinds of counters
    nucs_fwd = []
    nucs_rev = []
    quals = []
    hc_nucs = []
    nuc_qual = []
    alt_base_ratio_fwd = 0
    alt_base_ratio_rev = 0
    ym_ticker = []

    for pileupread in pileupcolumn.pileups:
        readpos = pileupread.query_position
        read = pileupread.alignment

        # broad try except, as this is not core functionality
        try:
            tags = read.get_tags()
            ym = [x[1] for x in tags if x[0] == "YM"][0]

        except:
            ym = 0

        # collect data on all reads spanning position
        if not pileupread.is_del and not pileupread.is_refskip:
            nuc = read.query_sequence[readpos]
            qual = read.query_qualities[readpos]

            if read.is_forward:
                nucs_fwd.append(nuc)
                nuc_qual.append((nuc, qual, "F"))

            else:
                nucs_rev.append(nuc)
                nuc_qual.append((nuc, qual, "R"))

            ym_ticker.append((nuc, ym))
            if float(qual) > high_base_quality_cutoff:
                hc_nucs.append(nuc)
            quals.append(qual)

        else:
            if add_dels:
                nuc = "-"
                if read.is_forward:
                    nucs_fwd.append(nuc)

                else:
                    nucs_rev.append(nuc)

    # Calculate metrics
    hc_count = Counter(hc_nucs)
    if len(hc_count.most_common()) > 1:
        hc_ratio = (
            1
            / sum([x[1] for x in hc_count.most_common()])
            * hc_count.most_common()[1][1]
        )

    else:
        hc_ratio = 0

    total = pileupcolumn.n
    quals_mean = np.mean(quals)

    nucs = nucs_fwd + nucs_rev
    counts_fwd, counts_rev, counts = (
        Counter(nucs_fwd),
        Counter(nucs_rev),
        Counter(nucs),
    )
    counts_fwd_mc, counts_rev_mc, counts_mc = (
        counts_fwd.most_common(),
        counts_rev.most_common(),
        counts.most_common(),
    )

    # access counts to calc totals
    counted_nucs = sum([x[1] for x in counts_mc])
    counted_nucs_fwd = sum([x[1] for x in counts_fwd_mc])
    counted_nucs_rev = sum([x[1] for x in counts_rev_mc])

    # empty position:
    if counted_nucs == 0 or total == 0:
        return (assembly, alleles, vcf_entry)

    else:
        data_present = counts_mc[0][1] / total
        non_ref_ratio_filtered = 1 - (counts_mc[0][1] / counted_nucs)

        if len(counts_mc) == 1:
            # perfect positions
            alleles = (counts_mc[0][0], ".")
            alt_base_mean_qual = 0

        else:
            # only allows 2 alleles
            # alleles = [x[0] for x in counts_mc]
            alleles = (counts_mc[0][0], counts_mc[1][0])
            base = counts_mc[0][0]
            alt_base = counts_mc[1][0]
            alt_base_mean_qual = np.mean([x[1] for x in nuc_qual if x[0] == alt_base])
            # calculate
            ticker = [x[1] for x in ym_ticker if x[0] == base]
            alt_ticker = [x[1] for x in ym_ticker if x[0] == alt_base]
            obs_ratio = 1 / (sum(ticker) + sum(alt_ticker)) * sum(alt_ticker)

        # Check support fwd
        if len(counts_fwd_mc) > 1:
            alt_base_fwd = counts_fwd_mc[1][0]
            fwd_count = counts_fwd_mc[1][1]
            alt_base_ratio_fwd = fwd_count / counted_nucs_fwd

        else:
            alt_base_fwd = None
            fwd_count = 0
            alt_base_ratio_fwd = 0

        # Check support reverse
        if len(counts_rev_mc) > 1:
            alt_base_rev = counts_rev_mc[1][0]
            rev_count = counts_rev_mc[1][1]
            alt_base_ratio_rev = rev_count / counted_nucs_rev

        else:
            alt_base_rev = None
            rev_count = 0
            alt_base_ratio_rev = 0

        base = None
        if alt_base_fwd == alt_base_rev and alt_base_fwd:
            base = alt_base_fwd

        tot_count = fwd_count + rev_count
        tot_ratio = np.mean((alt_base_ratio_fwd, alt_base_ratio_rev))

        vcf_entry.DP = total
        vcf_entry.DPQ = counted_nucs
        vcf_entry.FREQ = data_present
        vcf_entry.VAF = non_ref_ratio_filtered
        vcf_entry.FWDC = fwd_count
        vcf_entry.FWDR = alt_base_ratio_fwd
        vcf_entry.REVC = rev_count
        vcf_entry.REVR = alt_base_ratio_rev
        vcf_entry.TOTC = tot_count
        vcf_entry.TOTR = tot_ratio
        vcf_entry.SAME = 1 if base else 0
        vcf_entry.OBSR = obs_ratio
        vcf_entry.ABQ = alt_base_mean_qual
        vcf_entry.OBQ = quals_mean
        vcf_entry.HCR = hc_ratio

    return (assembly, alleles, vcf_entry)


def main(bam: Path, bed: Path, output_path: Path, pileup_depth=1_000_000):
    """
    Run SNP detection over given positions.

    Args:
        bam: Read alignments, BAM.
        bed: Genomic locations, BED.
        output_path: Output path for SNPs, VCF.
        pileup_depth: Maximum pileup depth, integer (DEFAULT=1_000_000).
    """

    from vcf_tools import create_bed_positions, initialize_output_vcf, write_vcf_entry
    from tqdm import tqdm
    import time

    # logging.debug("started main")
    # Open input files and create empty output VCF
    bam_af = pysam.AlignmentFile(bam, "rb")
    vcf = initialize_output_vcf(output_path, bam_af.references)

    # Iterate over positions in search space indicated in BED file
    for contig, pos, amplicon_ending in tqdm(create_bed_positions(bed)):
        # Check statistics for this position,
        # potentially finding a new variant
        positional_pileup = bam_af.pileup(
            contig,
            pos,
            pos + 1,
            truncate=True,
            max_depth=pileup_depth,
            min_base_quality=10,
            stepper="all",
        )

        for pileupcolumn in positional_pileup:
            result = extract_snp_evidence(pileupcolumn=pileupcolumn, assembly=contig)

            if result:
                # Reference allele is not '.', then a variant was found
                # The 'start' value is 0-based, 'stop' is 1-based
                r = vcf.new_record(
                    contig=contig, start=pos, alleles=result[1], filter="PASS"
                )

                # Write found variant as a new entry to VCF output
                for fld in fields(result[2]):
                    fld_value = getattr(result[2], fld.name)
                    if type(fld_value) in [float, np.float64, np.float32]:
                        fld_entry = str(f"{getattr(result[2], fld.name):.6f}")
                    elif type(fld_value) == int:
                        fld_entry = str(f"{getattr(result[2], fld.name)}")
                    else:
                        fld_entry = str(fld_value)

                    r.samples["Sample1"][fld.name] = fld_entry

                vcf.write(r)

            else:
                # Reference allele is '.', then no variant was found
                # Don't write anthing to VCF output file
                continue

    time.sleep(0.5)
    vcf.close()


if __name__ == "__main__":
    import argparse

    dev = False
    if not dev:
        parser = argparse.ArgumentParser(
            description="Process the information in the sequencing summary and add it to the bam."
        )

        parser.add_argument("variant_bed", type=Path)
        parser.add_argument("bam", type=Path)
        parser.add_argument("file_out", type=Path)
        args = parser.parse_args()
        # logging.info(args)

        main(args.bam, args.variant_bed, args.file_out)

    if dev:
        # PNK_01
        bed = Path("/data/projects/ROD_1125_variant_improvements/EGFR.bed")
        bam = Path(
            "/data/projects/ROD_1125_variant_improvements/ONT_20221121_EGFR/consensus_aligned/284.taged.bam"
        )
        vcf_out = Path("/data/projects/ROD_1125_variant_improvements/testsnp_EGFR.vcf")

        main(bam, bed, vcf_out)
