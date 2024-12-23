#!/usr/bin/env python
import logging
import time
from dataclasses import dataclass, fields
from typing import Dict

import numpy as np
import pysam


@dataclass
class VCF_entry:
    """
    Stores variant evidence annotations for a new VCF entry.
    """

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


def create_ref_mapper(reference_fasta, contig, start, stop) -> Dict[int, str]:
    """
    Code to create a dict with pos -> base for a given contig.
    """
    try:
        reference = pysam.FastaFile(reference_fasta)
    except OSError:
        try:
            time.sleep(0.2)
            reference = pysam.FastaFile(reference_fasta)
        except:
            raise OSError(f"Reference genome {reference_fasta} could not be read.")

    ref_mapper = {}
    for pos in range(start, stop + 1):
        ref_nuc = str(reference.fetch(reference=contig, start=pos, end=pos + 1)).upper()
        ref_mapper[pos] = ref_nuc

    return ref_mapper


def create_bed_lines(bed_file):
    """
    Given a bedfile yield the line as an array.
    """
    with open(bed_file) as f:
        for line in f:
            # bed should be tab delimited
            L = line.strip().split("\t")
            # try space splitting
            if len(L) == 1:
                L = L[0].split(" ")
                L = [x for x in L if x]
            yield L


def create_bed_positions(bed_file, end_warning_length=4):
    """
    Takes a BED file and yields a list of positions between the start and stop of each entry.
    """
    for L in create_bed_lines(bed_file):
        if not is_intable(L[2]):
            logging.critical("error in bed file")
        for pos in range(int(L[1]), int(L[2])):
            close_to_edge = (
                # close to start
                pos - end_warning_length <= int(L[1])
                or
                # close to end
                pos + end_warning_length >= int(L[2])
            )
            yield L[0], pos, close_to_edge


def write_vcf_entry(vcf, contig, pos, vcf_entry):
    """
    Writes variant evidence as an entry to a VCF file.

    Args:
        vcf: A pysam.VariantFile to write evidence variants to.
        contig: Contig name to write to VCF entry, column "CHROM", str.
        pos: Position in reference where the variant occurs, column "POS", int.
        vcf_entry: A tuple containing  the contig name, a tuple with reference and
            alternative alleles, and a VCF_entry object with variant annotations to add.
    """
    pos_alleles = vcf_entry[1]
    variant_format = vcf_entry[2]
    variant_info = vcf_entry[3]

    r = vcf.new_record(
        contig=contig, start=pos, stop=pos + 1, alleles=pos_alleles, filter="PASS"
    )

    for fld in fields(variant_format):
        fld_value = getattr(variant_format, fld.name)
        if type(fld_value) in [float, np.float64, np.float32]:
            fld_entry = str(f"{getattr(variant_format, fld.name):.6f}")
        elif type(fld_value) == int:
            fld_entry = str(f"{getattr(variant_format, fld.name)}")
        else:
            fld_entry = str(fld_value)

        r.samples["SAMPLE1"][fld.name] = fld_entry

    # Write info tag values to VCF output
    # These values will be used to filter the VCF
    for tag_id, value in variant_info.items():
        r.info[tag_id] = value

    vcf.write(r)


def initialize_output_vcf(vcf_path, contigs):
    """
    Initializes the VCF file where variant evidence will be written to.

    Args:
        vcf_path: Path, including file name, where the new VCF file will be created.
        contigs: Contig names, i.e. pysam references, where variants will be found,
            to be headed to VCF header.

    Returns:
        A writable pysam.VariantFile where variants evidence can be written.
    """
    # Create a VCF header
    vcfh = pysam.VariantHeader()
    # Add a sample named "ahstram" to our VCF header
    vcfh.add_sample("SAMPLE1")
    # Add a "FILTER" value, other than "PASS"
    vcfh.add_meta(
        "FILTER",
        items=[("ID", "RF"), ("Description", "Variant failed filter due to low RF")],
    )
    # Add a contig (chromosome 1) to our VCF header
    for contig in contigs:
        vcfh.add_meta("contig", items=[("ID", contig)])

    # Add to FORMAT in our VCF header
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
            ("Description", "Frequency, raw frequency of the reference allele"),
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

    # Add to INFO in our VCF header
    vcfh.add_meta(
        "INFO",
        items=[
            ("ID", "QA"),
            ("Number", 1),
            ("Type", "Float"),
            ("Description", "Alternate allele quality sum in phred"),
        ],
    )
    vcfh.add_meta(
        "INFO",
        items=[
            ("ID", "AO"),
            ("Number", 1),
            ("Type", "Integer"),
            ("Description", "Count of full observations of this alternate haplotype."),
        ],
    )
    vcfh.add_meta(
        "INFO",
        items=[
            ("ID", "DP"),
            ("Number", 1),
            ("Type", "Integer"),
            ("Description", "Total read depth at the locus"),
        ],
    )
    vcfh.add_meta(
        "INFO",
        items=[
            ("ID", "TYPE"),
            ("Number", 1),
            ("Type", "String"),
            ("Description", "The type of allele, either snp, ins, or del."),
        ],
    )
    vcfh.add_meta(
        "INFO",
        items=[
            ("ID", "SAF"),
            ("Number", 1),
            ("Type", "Integer"),
            ("Description", "Number of alternate observations on the forward strand"),
        ],
    )
    vcfh.add_meta(
        "INFO",
        items=[
            ("ID", "SAR"),
            ("Number", 1),
            ("Type", "Integer"),
            ("Description", "Number of alternate observations on the reverse strand"),
        ],
    )
    vcfh.add_meta(
        "INFO",
        items=[
            ("ID", "END"),
            ("Number", 1),
            ("Type", "Integer"),
            ("Description", "Stop position of the interval"),
        ],
    )

    vcf = pysam.VariantFile(vcf_path, "w", header=vcfh)
    return vcf
