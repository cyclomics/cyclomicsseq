#!/usr/bin/env python

from abc import abstractmethod
from collections import OrderedDict
from pathlib import Path
from typing import Any, Tuple
import pandas as pd
import io
import requests
import sys
import json


class VCF_file:
    def __init__(self, vcf_file):
        self.vcf_file = vcf_file
        self.vcf_header = ""
        self.vcf = self.read_vcf(self.vcf_file)

    @staticmethod
    def relaxed_float(x: Any) -> float:
        """Return a float, with value error catch"""
        try:
            my_float = float(x)
        except ValueError:
            my_float = float(0)
        return my_float

    def read_vcf(self, path: Path) -> pd.DataFrame:
        """Read in a VCF file and return it as Pandas DataFrame"""
        with open(path, "r") as f:
            header = []
            lines = []
            for l in f:
                if l.startswith("##"):
                    header.append(l)
                else:
                    lines.append(l)

        self.vcf_header = header
        df = pd.read_csv(
            io.StringIO("".join(lines)),
            dtype={
                "#CHROM": str,
                "POS": int,
                "ID": str,
                "REF": str,
                "ALT": str,
                "QUAL": str,
                "FILTER": str,
                "INFO": str,
            },
            sep="\t",
        ).rename(columns={"#CHROM": "CHROM"})

        if not df.empty:
            formats = df.FORMAT[0].split(":")
            for i, fmt in enumerate(formats):
                df[fmt] = df.Sample1.apply(
                    lambda x: self.relaxed_float(x.split(":")[i])
                    if (x.split(":")[i])
                    else 0
                )

        return df

    def write(self, path: Path):
        """Write output VCF file"""
        with open(path, "w") as new_vcf:
            new_vcf.writelines(self.vcf_header)

            writeable_vcf = self.vcf.rename(columns={"CHROM": "#CHROM"})
            writeable_vcf = writeable_vcf[
                [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "INFO",
                    "FORMAT",
                    "Sample1",
                ]
            ]

            new_vcf.writelines(writeable_vcf.to_csv(sep="\t", index=False))

    @staticmethod
    def get_query_allele_positions(
        ref_allele: str, alt_allele: str, start: int
    ) -> Tuple[int, int]:
        """Determine start and end positions for Ensembl query

        Positions are determined based on reference and alternative
        alleles, which inform if a variant is a SNP, an 1 or a
        deletion. Index adjustments have to be made for indels.
        """

        if len(ref_allele) > len(alt_allele):
            # this is a deletion
            # remove first base of ref which is maintained and not part of the deletion
            end = start + len(ref_allele) - 1
            start += 1
            alt_allele = "-"

        elif len(ref_allele) < len(alt_allele):
            # This is an insertion
            end = start
            # remove first base of ref which is maintained and not part of the deletion
            start += 1
            ref_allele = "-"

        elif len(ref_allele) == len(alt_allele):
            # This is a snp (ignoring insdel events)
            end = start

        # Ensembl-VEP needs REF to be the strand, which is always forward (1)
        ref_allele = "1"

        return (start, end, ref_allele, alt_allele)

    @staticmethod
    def ensembl_vep(query: str) -> dict:
        """Query Ensembl-VEP through API to annotate variants"""
        try:
            response = requests.get(
                url=query, headers={"Content-Type": "application/json"}
            )

            json_data = json.loads(response.text)[0]

        except:
            # Any query that doesn't exist in Ensembl will return an SLLError or HTTPError,
            # e.g. querying variants in a backbone sequence.
            # In that case, we return nothing
            return

        return json_data

    @staticmethod
    def get_annotation_text(vep_json: dict) -> str:
        """Parse VEP response dict and find relevant annotations"""
        # Initialize variant annotations
        variant_class = None
        consequence = None
        mutation_ids = None
        cosmic_legacy_ids = None
        gene = None
        impact = None
        biotype = None
        amino_acids = None
        canonical = None
        sift = None
        polyphen = None

        if not vep_json:
            # Ensembl-VEP query did not return any response
            # e.g. because the variant was in a backbone sequence
            annot_text = "."
            return annot_text

        # Find relevant annotations in response dict
        variant_class = vep_json.get("variant_class")
        consequence = vep_json.get("most_severe_consequence")

        # Try-Except block: cannot always find colocated variants
        # and thus COSMIC IDs
        try:
            colocated_variants = vep_json.get("colocated_variants")
            mutation_ids = []
            cosmic_legacy_ids = []
            for xref in colocated_variants:
                mutation_ids.append(xref["id"])
                if xref["allele_string"] == "COSMIC_MUTATION":
                    cosmic_legacy_ids.append(xref["var_synonyms"]["COSMIC"][0])

            # Join list of found IDs into comma-separated string
            mutation_ids = ",".join(mutation_ids)
            cosmic_legacy_ids = ",".join(cosmic_legacy_ids)
        except:
            mutation_ids = "None"
            cosmic_legacy_ids = "None"

        transcript_cons = vep_json.get("transcript_consequences")
        # Transcript consequences can differ a lot per query
        # If something is not found, will be returned as 'None'
        if transcript_cons:
            gene = transcript_cons[0].get("gene_symbol")
            impact = transcript_cons[0].get("impact")
            biotype = transcript_cons[0].get("biotype")
            amino_acids = transcript_cons[0].get("amino_acids")
            canonical = transcript_cons[0].get("canonical")

            sift_prediction = transcript_cons[0].get("sift_prediction")
            sift_score = transcript_cons[0].get("sift_score")
            if sift_prediction:
                sift = f"{sift_prediction}({sift_score})"

            polyphen_prediction = transcript_cons[0].get("polyphen_prediction")
            polyphen_score = transcript_cons[0].get("polyphen_score")
            if polyphen_prediction:
                polyphen = f"{polyphen_prediction}({polyphen_score})"

        # Merge all annotations into a string to be returned
        annot_dict = OrderedDict(
            {
                "variant_class": variant_class,
                "consequence": consequence,
                "COSMIC": mutation_ids,
                "COSMIC_legacy": cosmic_legacy_ids,
                "gene": gene,
                "impact": impact,
                "biotype": biotype,
                "amino_acids": amino_acids,
                "canonical": canonical,
                "SIFT": sift,
                "PolyPhen": polyphen,
            }
        )

        annot_text = ";".join([f"{k}={v}" for k, v in annot_dict.items()])
        # If annotation is None, don't print it
        # annot_text = ";".join([f"{k}={v}" for k, v in annot_dict.items() if v])

        return annot_text

    def annotate_vep(self, server: str):
        """Annotate a set of variants from a VCF file with Ensembl-VEP

        Input: Server URL (string), e.g. "https://rest.ensembl.org"
        """
        # Set API search options
        params = "canonical=1&variant_class=1&hgvs=1&vcf_string=1&pick=1"

        if self.vcf.empty:
            # There are no variants to annotate
            return

        annotations = []
        # Loop over variants in VCF file, annotate one at a time
        # TODO: Parallelize
        for var in self.vcf.iterrows():
            chromosome = var[1]["CHROM"]
            start = var[1]["POS"]
            ref_allele = var[1]["REF"]
            alt_allele = var[1]["ALT"]

            start, end, ref_allele, alt_allele = self.get_query_allele_positions(
                ref_allele, alt_allele, start
            )

            query = (
                f"{server}/vep/human/region/"
                f"{chromosome}:{start}-{end}:"
                f"{ref_allele}/{alt_allele}?"
                f"{params}?"
            )

            # Query Ensembl with the VEP API, returns a JSON dict
            vep_json = self.ensembl_vep(query)
            # Parse JSON response to get annotations in a string
            annotation_text = self.get_annotation_text(vep_json)
            # Add to annotations list
            annotations.append(annotation_text)

        # Write annotations to INFO column in VCF output file
        self.vcf["INFO"] = annotations


if __name__ == "__main__":
    dev = False
    if not dev:
        import argparse

        parser = argparse.ArgumentParser(description="Filter a vcf")

        parser.add_argument("variant_vcf", type=Path)
        parser.add_argument("file_out", type=Path)
        args = parser.parse_args()

        # Can be added to argparse
        server = "https://rest.ensembl.org"

        vcf = VCF_file(args.variant_vcf)
        vcf.annotate_vep(server)
        vcf.write(args.file_out)

    if dev:
        # variant_vcf = "/scratch/nxf_work/rodrigo/a4/a030f4325cb5bf3385a1ceb02f6c3a/FAS12641.merged.vcf"
        variant_vcf = "testindel_EGFR.vcf"

        server = "https://rest.ensembl.org"

        vcf = VCF_file(variant_vcf)
        vcf.annotate_vep(server)
        vcf.write("testannotate_EGFR.vcf")
