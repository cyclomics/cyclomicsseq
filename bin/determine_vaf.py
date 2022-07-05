#!/usr/bin/env python

import logging
from collections import Counter, namedtuple
from pathlib import Path

import pysam
from tqdm import tqdm


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
    ]
    Variant = namedtuple('Variant', tags)

    data_present = 0
    alt_base_fwd = None
    fwd_count = 0
    alt_base_ratio_fwd = 0
    alt_base_rev = None
    rev_count = 0
    alt_base_ratio_rev = 0
    alleles = (".", ".")
    non_ref_ratio_filtered = 1
    total = 0
    base = None
    counted_nucs = 0

    for pileupcolumn in bam.pileup(
        assembly,
        pos,
        pos + 1,
        truncate=True,
        max_depth=pileup_depth,
        min_base_quality=minimum_base_quality,
    ):
        nucs_fwd = []
        nucs_rev = []
        alt_base_ratio_fwd = 0
        alt_base_ratio_rev = 0
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                nuc = pileupread.alignment.query_sequence[pileupread.query_position]
                if pileupread.alignment.is_forward:
                    nucs_fwd.append(nuc)
                else:
                    nucs_rev.append(nuc)

            else:
                if add_dels:
                    nuc = "-"
                    if pileupread.alignment.is_forward:
                        nucs_fwd.append(nuc)
                    else:
                        nucs_rev.append(nuc)

        nucs = nucs_fwd + nucs_rev
        counts_fwd, counts_rev, counts = (
            Counter(nucs_fwd),
            Counter(nucs_rev),
            Counter(nucs),
        )

        counted_nucs = sum([x[1] for x in counts.most_common()])
        counted_nucs_fwd = sum([x[1] for x in counts_fwd.most_common()])
        counted_nucs_rev = sum([x[1] for x in counts_rev.most_common()])

        total = pileupcolumn.n
        
        # empty position:
        if counted_nucs ==0:
            continue

        else:   
            if total > 0 and len(counts.most_common()) > 0:
                data_present = counts.most_common()[0][1] / total
            else:
                data_present = 0

            if counted_nucs > 0:
                non_ref_ratio_filtered = 1 - (counts.most_common()[0][1] / counted_nucs)
            else:
                non_ref_ratio_filtered = 1

            if len(counts.most_common()) > 1:
                alleles = (counts.most_common()[0][0], counts.most_common()[1][0])
            elif len(counts.most_common()) == 1:
                # perfect positions
                alleles = (counts.most_common()[0][0], ".")
            else:
                alleles = (".", ".")

            if len(counts_fwd.most_common()) > 1:
                alt_base_fwd = counts_fwd.most_common()[1][0]
                fwd_count = counts_fwd.most_common()[1][1]
                alt_base_ratio_fwd = fwd_count / counted_nucs_fwd
            else:
                alt_base_fwd = None
                fwd_count = 0
                alt_base_ratio_fwd = 0


            if len(counts_rev.most_common()) > 1:
                alt_base_rev = counts_rev.most_common()[1][0]
                rev_count = counts_rev.most_common()[1][1]
                alt_base_ratio_rev = rev_count / counted_nucs_rev
            else:
                alt_base_rev = None
                rev_count = 0
                alt_base_ratio_rev = 0


            base = None
            if alt_base_fwd == alt_base_rev and alt_base_fwd:
                base = alt_base_fwd

        # create an object to store all gathered data
    var = Variant(DP=total, DPQ=counted_nucs, FREQ= data_present , VAF= non_ref_ratio_filtered, FWDC= fwd_count, FWDR= alt_base_ratio_fwd, REVC= rev_count, REVR= alt_base_ratio_rev, TOTC= '', TOTR= '', SAME = (1 if base else 0)) 
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
    return (
        assembly,
        alleles,
        var
    )


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
            ("Description", "Depth at position",),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "DPQ"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Depth at position above Q threshold",),
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
            ("Description", "Forward count after Q filtering",),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "FWDR"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Forward ratio after Q filtering",),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "REVC"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Reverse count after Q filtering",),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "REVR"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Reverse ratio after Q filtering",),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "TOTC"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Total count after Q filtering",),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "TOTR"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Total ratio after Q filtering",),
        ],
    )
    vcfh.add_meta(
        "FORMAT",
        items=[
            ("ID", "SAME"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "Same base found forward and reverse",),
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
                fld_entry =  str(f"{getattr(result[2], fld):.4f}")
            elif type(fld_value) == int:
                fld_entry =  str(f"{getattr(result[2], fld)}")
            else:
                fld_entry = str(fld_value)

            r.samples["Sample1"][fld] = fld_entry
        
        vcf.write(r)
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
        bed = Path(
            "/media/dami/cyclomics_003/raw_data/ONT/PNK_05_Rob/variant_locations/real_variants_chm13-v2_0_panels.bed"
        )
        bam = Path(
            "/media/dami/cyclomics_003/results/ONT/PNK_01_rob/Annotate/AnnotateBamXTags/FAS16079.annotated.bam"
        )
        main(bam, bed, "tmp.vcf")
