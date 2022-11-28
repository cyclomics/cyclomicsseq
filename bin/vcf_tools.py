import logging
from dataclasses import dataclass, fields

import numpy as np
import pysam


@dataclass
class VCF_entry:
    DP: int = 0
    DPQ: int = 0
    FREQ: float = 0.0
    VAF: float = 0.0
    FWDC: int = 0
    FWDR: float = 0.0
    REVC: int = 0
    REVR: float = 0.0
    TOTC: int = 0
    TOTR: float = 0.0
    SAME: int = 0
    OBSR: float = 0.0
    ABQ: float = 0.0
    OBQ: float = 0.0
    HCR: int = 0


def is_intable(value):
    """
    Test if the input can be converted to an int.
    """
    try:
        float(value)
        return True
    except ValueError:
        return False


def create_bed_positions(bed_file, end_warning_length=4):
    with open(bed_file) as f:
        for line in f:
            print(line)
            # bed should be tab delimited
            L = line.strip().split("\t")
            # try space splitting
            if len(L) == 1:
                L = L[0].split(" ")
                L = [x for x in L if x]

            if not is_intable(L[2]):
                logging.critical("error in bed file")
            for pos in range(int(L[1]), int(L[2])):
                close_to_end = pos + end_warning_length >= int(L[2])
                yield L[0], pos, close_to_end


def write_vcf_entry(vcf, contig, pos, vcf_entry):
    r = vcf.new_record(
        contig=contig, start=pos, stop=pos + 1, alleles=vcf_entry[1], filter="PASS"
    )

    for fld in fields(vcf_entry[2]):
        fld_value = getattr(vcf_entry[2], fld.name)
        if type(fld_value) in [float, np.float64, np.float32]:
            fld_entry = str(f"{getattr(vcf_entry[2], fld.name):.6f}")
        elif type(fld_value) == int:
            fld_entry = str(f"{getattr(vcf_entry[2], fld.name)}")
        else:
            fld_entry = str(fld_value)

        r.samples["Sample1"][fld.name] = fld_entry

    vcf.write(r)


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
