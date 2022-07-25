#!/usr/bin/env python

import logging
from collections import Counter, namedtuple
from pathlib import Path
import time

import pysam
from tqdm import tqdm
import numpy as np


def is_intable(value):
    """
    Test if the input can be converted to an int.
    """
    try:
        float(value)
        return True
    except ValueError:
        return False


def create_bed_positions(bed_file):
    with open(bed_file) as f:
        for line in f:
            L = line.strip().split("\t")
            if not is_intable(L[2]):
                logging.critical("error in bed file")
            for pos in range(int(L[1]), int(L[2])):
                yield L[0], pos


def extract_nucleotide_count(
    bam,
    assembly,
    pos,
    pileup_depth,
    minimum_base_quality=10,
    add_dels=False,
    high_base_quality_cutoff=80,
):
    tags = [
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
    Variant = namedtuple("Variant", tags)
    alleles = (".", ".")

    positional_generator = bam.pileup(
        assembly,
        pos,
        pos + 1,
        truncate=True,
        max_depth=pileup_depth,
        min_base_quality=minimum_base_quality,
    )
    for pileupcolumn in positional_generator:
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
            tags = read.get_tags()
            # broad try except, as this is not core functionality
            try:
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
                quals.append(read.query_qualities[readpos])
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
        quals_mean = np.mean(qual)

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
                alleles = (counts_mc[0][0], counts_mc[1][0])
                base = counts_mc[0][0]
                alt_base = counts_mc[1][0]
                alt_base_mean_qual = np.mean(
                    [x[1] for x in nuc_qual if x[0] == alt_base]
                )
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
            # print(test_var)

    
    # create an object to store all gathered data
    var = Variant(
        DP=total if total else 0,
        DPQ=counted_nucs if counted_nucs else 0,
        FREQ=data_present if data_present else 0,
        VAF=non_ref_ratio_filtered if non_ref_ratio_filtered else 0,
        FWDC=fwd_count if fwd_count else 0,
        FWDR=alt_base_ratio_fwd if alt_base_ratio_fwd else 0,
        REVC=rev_count if rev_count else 0,
        REVR=alt_base_ratio_rev if alt_base_ratio_rev else 0,
        TOTC=tot_count if tot_count else 0  ,
        TOTR= tot_ratio if tot_ratio else 0 ,
        SAME=(1 if base else 0),
        OBSR=obs_ratio if obs_ratio else 0,
        ABQ=alt_base_mean_qual if alt_base_mean_qual else 0,
        OBQ= quals_mean if quals_mean else 0,
        HCR=hc_ratio if hc_ratio else 0,
    )
    # print(var)
    # print("coverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
    # print(f"found positions {counted_nucs}")
    # print(f"raw counts: {counts}")
    # print(f"raw counts fwd: {counts_fwd}")
    # print(f"raw counts rev: {counts_rev}")

    # print(f"non reference_ratio : {data_present} ({data_present * 100:.2f}%)")
    # print(
    #     f"alternative ratio   : {non_ref_ratio_filtered} ({non_ref_ratio_filtered * 100:.2f}%)"
    # )

    # print(var)
    return (assembly, alleles, var)


def initialize_output_vcf(vcf_path, contigs):

    # Create a VCF header
    vcfh = pysam.VariantHeader()
    # Add a sample named "ahstram" to our VCF header
    vcfh.add_sample("Sample1")
    # Add a "FILTER" value, other than "PASS"
    vcfh.add_meta(
        "FILTER",
        items=[("ID", "RF"), ("Description", "Variant failed filter due to low RF")],
    )
    # Add a contig (chromosome 1) to our VCF header
    for contig in contigs:
        vcfh.add_meta("contig", items=[("ID", contig)])

    # Add GT to FORMAT in our VCF header
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "DP"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Depth at position"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "DPQ"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Depth at position above Q threshold"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "FREQ"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Frequency, raw frequency of non majority allele"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "VAF"),
            ("Number", 1),
            ("Type", "String"),
            (
                "Description",
                "Variant Allele Frequency, filtered frequency on base quaility of non majority allele",
            ),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "FWDC"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Forward count after Q filtering"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "FWDR"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Forward ratio after Q filtering"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "REVC"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Reverse count after Q filtering"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "REVR"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Reverse ratio after Q filtering"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "TOTC"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Total count after Q filtering"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "TOTR"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Total ratio after Q filtering"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "SAME"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Same base found forward and reverse"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "OBSR"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "observational ratio of alternative allele"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "ABQ"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "alternative base quality (mean)"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "OBQ"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "overall base quality (mean)"),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "HCR"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "high quality ratio"),
        ],
    )

    vcf = pysam.VariantFile(vcf_path, "w", header=vcfh)
    return vcf


def main(bam: Path, variants: Path, output_path, pileup_depth=1_000_000):
    logging.debug("started main")
    bam_af = pysam.AlignmentFile(bam, "rb")

    vcf = initialize_output_vcf(output_path, bam_af.references)

    for assembly, pos in tqdm(create_bed_positions(variants)):
        # print(f"{assembly}:{pos}")
        result = extract_nucleotide_count(bam_af, assembly, pos, pileup_depth)
        # print(result)
        # For a vcf file we need to convert 0 index to 1 index

        # Create a record at chr1:1000 A/T which failed filter due to "RF"
        # The 'start' value is 0-based, 'stop' is 1-based
        r = vcf.new_record(
            contig=assembly, start=pos, stop=pos + 1, alleles=result[1], filter="PASS"
        )
        for fld in result[2]._fields:
            fld_value = getattr(result[2], fld)
            if type(fld_value) == float:
                fld_entry = str(f"{getattr(result[2], fld):.6f}")
            elif type(fld_value) == int:
                fld_entry = str(f"{getattr(result[2], fld)}")
            else:
                fld_entry = str(fld_value)

            r.samples["Sample1"][fld] = fld_entry

        vcf.write(r)
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
        logging.info(args)

        main(args.bam, args.variant_bed, args.file_out)

    if dev:
        # bed = Path("dilution_series_expected_mutations.bed")
        bed = Path("/home/dami/projects/variantcalling/depth/TP53_S0.bed")
        bam = Path("/home/dami/Data/dev/000016/FAT55661.YM_gt_5.bam")
        main(bam, bed, "tmp.vcf")
