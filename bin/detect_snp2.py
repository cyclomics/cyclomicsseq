#!/usr/bin/env python
from collections import Counter
from pathlib import Path

import numpy as np
from vcf_tools import VCF_entry


def extract_snv_evidence(
    pileupcolumn,
    assembly,
    add_dels=False,
    high_base_quality_cutoff=80,
):
    """Returns found variants in a given pileup position"""

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
        continue

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

    vcf_entry = VCF_entry(
        DP=total,
        DPQ=counted_nucs,
        FREQ=data_present,
        VAF=non_ref_ratio_filtered,
        FWDC=fwd_count,
        FWDR=alt_base_ratio_fwd,
        REVC=rev_count,
        REVR=alt_base_ratio_rev,
        TOTC=tot_count,
        TOTR=tot_ratio,
        SAME=(1 if base else 0),
        OBSR=obs_ratio,
        ABQ=alt_base_mean_qual,
        OBQ=quals_mean,
        HCR=hc_ratio,
    )

    return (assembly, alleles, vcf_entry)


if __name__ == "__main__":
    import argparse

    dev = True
    if not dev:
        parser = argparse.ArgumentParser(
            description="Process the information in the sequencing summary and add it to the bam."
        )

        parser.add_argument("variant_bed", type=Path)
        parser.add_argument("bam", type=Path)
        parser.add_argument("file_out", type=Path)
        args = parser.parse_args()
        # logging.info(args)

        extract_snv_evidence(args.bam, args.variant_bed, args.file_out)

    if dev:
        # PNK_01
        bed = Path("tmp/PNK_Rob_custom_GRCh38.bed")
        bam = Path("tmp/PNK_01_GRCh38.p14/FAS12641.taged.bam")
        vcf_out = Path("tmp/PNK_01_GRCh38.p14/PNK_01_GRCh38.p14.snp.vcf")

        extract_snv_evidence(bam, bed, vcf_out)
