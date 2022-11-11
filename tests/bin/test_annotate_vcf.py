import unittest
from collections import defaultdict
from math import exp
from unittest.mock import MagicMock, Mock

from bin import annotate_vcf as target


class TestSmokeTest(unittest.TestCase):
    def test_Some(self):
        """
        test should always pass.
        """
        self.assertEqual("foo".upper(), "FOO")


class TestGetQueryPositions(unittest.TestCase):
    def test_deletion(self):
        """Test if deletion positions are correctly parsed before query."""
        # test_params1 = {'ref_allele':'ATG', 'alt_allele':'A', 'start':2}
        test_params1_result = tuple([2, 3, "1", "-"])
        self.assertEqual(
            target.VCF_file.get_query_allele_positions("ATG", "A", 1),
            test_params1_result,
        )

    def test_insertion(self):
        """Test if insertion positions are correctly parsed before query."""
        # test_params2 = '[ref_allele:A, alt_allele:ATG, start:2]'
        test_params2_result = tuple([2, 1, "1", "ATG"])
        self.assertEqual(
            target.VCF_file.get_query_allele_positions("A", "ATG", 1),
            test_params2_result,
        )

    def test_snp(self):
        """Test if SNP positions are correctly parsed before query."""
        # test_params3 = '[ref_allele:A, alt_allele:G, start:1]'
        test_params3_result = tuple([1, 1, "1", "G"])
        self.assertEqual(
            target.VCF_file.get_query_allele_positions("A", "G", 1), test_params3_result
        )


class TestEnsemblVep(unittest.TestCase):
    # TODO: mock response.get
    def test_api_real(self):
        """ " Test if API can get the correct response for a real variant."""
        query1 = "https://rest.ensembl.org/vep/human/region/CM000669.2:55165450-55165450:1/C?canonical=1&variant_class=1&hgvs=1&vcf_string=1&pick=1?"

        test_query1_result = {
            "start": 55165450,
            "input": "7 55165450 55165450 T/C 1",
            "seq_region_name": "7",
            "assembly_name": "GRCh38",
            "variant_class": "SNV",
            "most_severe_consequence": "intron_variant",
            "vcf_string": "7-55165450-T-C",
            "transcript_consequences": [
                {
                    "gene_id": "ENSG00000146648",
                    "impact": "MODIFIER",
                    "hgvsc": "ENST00000275493.7:c.1880+13T>C",
                    "gene_symbol_source": "HGNC",
                    "hgnc_id": "HGNC:3236",
                    "transcript_id": "ENST00000275493",
                    "consequence_terms": ["intron_variant"],
                    "biotype": "protein_coding",
                    "variant_allele": "C",
                    "gene_symbol": "EGFR",
                    "strand": 1,
                    "canonical": 1,
                }
            ],
            "strand": 1,
            "id": "7_55165450_T/C",
            "end": 55165450,
            "allele_string": "T/C",
        }
        self.assertEqual(target.VCF_file.ensembl_vep(query1), test_query1_result)


class TestGetAnnotationText(unittest.TestCase):
    def test_get_annotation_text(self):
        """Test if API response is parsed into the correct annotation text."""
        json1 = {
            "start": 55165450,
            "input": "7 55165450 55165450 T/C 1",
            "seq_region_name": "7",
            "assembly_name": "GRCh38",
            "variant_class": "SNV",
            "most_severe_consequence": "intron_variant",
            "vcf_string": "7-55165450-T-C",
            "transcript_consequences": [
                {
                    "gene_id": "ENSG00000146648",
                    "impact": "MODIFIER",
                    "hgvsc": "ENST00000275493.7:c.1880+13T>C",
                    "gene_symbol_source": "HGNC",
                    "hgnc_id": "HGNC:3236",
                    "transcript_id": "ENST00000275493",
                    "consequence_terms": ["intron_variant"],
                    "biotype": "protein_coding",
                    "variant_allele": "C",
                    "gene_symbol": "EGFR",
                    "strand": 1,
                    "canonical": 1,
                }
            ],
            "strand": 1,
            "id": "7_55165450_T/C",
            "end": 55165450,
            "allele_string": "T/C",
        }
        test_json1_result = "variant_class=SNV;consequence=intron_variant;COSMIC=None;COSMIC_legacy=None;gene=EGFR;impact=MODIFIER;biotype=protein_coding;amino_acids=None;canonical=1;SIFT=None;PolyPhen=None"
        self.assertEqual(target.VCF_file.get_annotation_text(json1), test_json1_result)
