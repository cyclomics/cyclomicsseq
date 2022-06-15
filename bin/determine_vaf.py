from collections import Counter
from pathlib import Path
import pysam

import logging


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
    bam, assembly, pos, pileup_depth, minimum_base_quality=30, add_dels=False
):

    for pileupcolumn in bam.pileup(
        assembly,
        pos,
        pos + 1,
        truncate=True,
        max_depth=pileup_depth,
        min_base_quality=minimum_base_quality,
    ):
        nucs = []
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                nuc = pileupread.alignment.query_sequence[pileupread.query_position]
                nucs.append(nuc)
            else:
                if add_dels:
                    nucs.append("-")

        counts = Counter(nucs)
        counted_nucs = sum([x[1] for x in counts.most_common()])
        total = pileupcolumn.n
        non_ref_ratio = 1 - (counts.most_common()[0][1] / total)
        non_ref_ratio_filtered = 1 - (counts.most_common()[0][1] / counted_nucs)
        if len(counts.most_common()) > 1:
            alleles = (counts.most_common()[0][0], counts.most_common()[1][0])
        else:
            # perfect positions
            alleles = (counts.most_common()[0][0],".")
        print("coverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
        print(f"found positions {counted_nucs}")
        print(f"raw counts: {counts}")
        print(f"non reference_ratio : {non_ref_ratio} ({non_ref_ratio * 100:.2f}%)")
        print(
            f"alternative ratio   : {non_ref_ratio_filtered} ({non_ref_ratio_filtered * 100:.2f}%)"
        )
        print("")

        return (
            assembly,
            alleles,
            non_ref_ratio,
            non_ref_ratio_filtered,
            minimum_base_quality,
        )


def initialize_output_vcf(vcf_path):

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
    vcfh.add_meta("contig", items=[("ID", 1)])
    vcfh.add_meta("contig", items=[("ID", "chr17")])

    # Add GT to FORMAT in our VCF header
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

    vcf = pysam.VariantFile(vcf_path, "w", header=vcfh)
    return vcf


def main(bam: Path, variants: Path, pileup_depth=1_000_000):
    logging.debug("started main")
    bam_af = pysam.AlignmentFile(bam, "rb")

    vcf = initialize_output_vcf("example.vcf")

    for assembly, pos in create_bed_positions(variants):
        print(f"{assembly}:{pos}")
        result = extract_nucleotide_count(bam_af, assembly, pos, pileup_depth)
        print(result)
        # For a vcf file we need to convert 0 index to 1 index
        print("")

        # Create a record at chr1:1000 A/T which failed filter due to "RF"
        # The 'start' value is 0-based, 'stop' is 1-based
        r = vcf.new_record(
            contig=assembly, start=pos, stop=pos + 1, alleles=result[1], filter="PASS"
        )
        r.samples["Sample1"]["FREQ"] = f"{result[2]*100:.3f}%"
        r.samples["Sample1"]["VAF"] = f"{result[3]*100:.3f}%"

        vcf.write(r)


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

        main(args.bam, args.variants)

    if dev:
        perc_1_vars = [7675937, 7675939, 7675945, 7675984]
        perc_05_vars = [7675934, 7675948, 7675950, 7675977]
        perc_02_vars = [7675943, 7675954, 7675975, 7675986]

        # bed = Path("dilution_series_expected_mutations.bed")
        bed = Path("dilution_series_full_amplicon.bed")

        bam = Path(
            "/media/dami/cyclomics_003/results/Cyclomics/000001_v3/Annotate/AnnotateBamXTags/fastq.tag_annotated.sort.bam"
        )
        main(bam, bed)
