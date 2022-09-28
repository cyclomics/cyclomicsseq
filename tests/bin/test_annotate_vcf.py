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
        """ Test if deletion positions are correctly parsed before query."""
        #test_params1 = {'ref_allele':'ATG', 'alt_allele':'A', 'start':1}
        test_params1_result = tuple([1, 1])
        self.assertEqual(target.VCF_file.get_query_positions('ATG', 'A', 1), test_params1_result)
    
    def test_insertion(self):
        """ Test if insertion positions are correctly parsed before query."""
        #test_params2 = '[ref_allele:A, alt_allele:ATG, start:2]'
        test_params2_result = tuple([2, 1])
        self.assertEqual(target.VCF_file.get_query_positions('A', 'ATG', 2), test_params2_result)

    def test_snp(self):
        """ Test if SNP positions are correctly parsed before query."""
        #test_params3 = '[ref_allele:ATG, alt_allele:GTA, start:3]'
        test_params3_result = tuple([3, 3])
        self.assertEqual(target.VCF_file.get_query_positions('ATG', 'GTA', 3), test_params3_result)

class TestEnsemblVep(unittest.TestCase):
    # TODO: mock response.get
    def test_api_real(self):
        """" Test if API can get the correct response for a real variant."""
        query1 = (
                "https://rest.ensembl.org/vep/human/region/"
                "CM000663.2:114713908-114713908:"
                "1/C?"
                "canonical=1"
                "&variant_class=1"
                "&hgvs=1"
                "&vcf_string=1"
                "&pick=1?"
            )
        
        test_query1_result = {'strand': 1, 'id': '1_114713908_T/C', 'vcf_string': '1-114713908-T-C', 'assembly_name': 'GRCh38', 'end': 114713908, 'most_severe_consequence': 'missense_variant', 'colocated_variants': [{'allele_string': 'COSMIC_MUTATION', 'id': 'COSV54736340', 'start': 114713908, 'strand': 1, 'phenotype_or_disease': 1, 'somatic': 1, 'end': 114713908, 'seq_region_name': '1', 'var_synonyms': {'COSMIC': ['COSM584']}}, {'var_synonyms': {'COSMIC': ['COSM583']}, 'seq_region_name': '1', 'somatic': 1, 'phenotype_or_disease': 1, 'strand': 1, 'end': 114713908, 'start': 114713908, 'id': 'COSV54736624', 'allele_string': 'COSMIC_MUTATION'}, {'id': 'COSV54738969', 'allele_string': 'COSMIC_MUTATION', 'strand': 1, 'somatic': 1, 'phenotype_or_disease': 1, 'end': 114713908, 'start': 114713908, 'var_synonyms': {'COSMIC': ['COSM582']}, 'seq_region_name': '1'}, {'start': 114713908, 'end': 114713908, 'strand': 1, 'phenotype_or_disease': 1, 'somatic': 1, 'seq_region_name': '1', 'var_synonyms': {'COSMIC': ['COSM5044298']}, 'allele_string': 'COSMIC_MUTATION', 'id': 'COSV54747786'}, {'end': 114713908, 'strand': 1, 'phenotype_or_disease': 1, 'start': 114713908, 'clin_sig': ['likely_pathogenic', 'pathogenic'], 'var_synonyms': {'ClinVar': ['RCV000291285', 'VCV000280409', 'RCV000444899', 'RCV000444754', 'RCV000444660', 'RCV000444278', 'RCV000444188', 'RCV000441171', 'RCV000439765', 'RCV000439526', 'RCV000438738', 'RCV000438468', 'RCV000437545', 'RCV000437312', 'RCV000436539', 'RCV000435905', 'RCV000435412', 'RCV000434604', 'RCV000433761', 'RCV000433349', 'RCV000432170', 'RCV000431592', 'RCV000431333', 'RCV000430000', 'RCV000429082', 'RCV000428903', 'RCV000428055', 'RCV000427746', 'RCV000426654', 'RCV000426122', 'RCV000424084', 'RCV000423898', 'RCV000422640', 'RCV000421496', 'RCV000421291', 'RCV000420302', 'RCV000419201', 'RCV000419053', 'RCV000418396', 'RCV000418220', 'VCV000375874', 'RCV000014914', 'VCV000013900', 'RCV000037574', 'RCV000032847', 'RCV000435687', 'RCV000432961', 'RCV000420832', 'RCV000419710', 'RCV000431883', 'RCV000430593', 'RCV000430407', 'RCV000424960', 'RCV000424721', 'RCV000424455', 'RCV000422278', 'RCV000422078', 'RCV000413804', 'RCV000445249', 'RCV000439264', 'RCV000438052', 'RCV000441317', 'RCV000440367', 'RCV000114745', 'RCV000114744', 'RCV000148032'], 'UniProt': ['VAR_006847'], 'OMIM': [164790.0002]}, 'seq_region_name': '1', 'id': 'rs11554290', 'pubmed': [24033266, 26619011, 22499344, 1654209, 2278970, 2674680, 3122217, 6587382, 8120410, 12460918, 12727991, 14508525, 16273091, 16291983, 16434492, 17699718, 18390968, 18633438, 18948947, 19075190, 19880792, 20130576, 20149136, 20179705, 20406486, 20619739, 20736745, 21107323, 21305640, 21576590, 21729679, 21829508, 22407852, 22761467, 22773810, 23392294, 23414587, 23515407, 23538902, 23569304, 23614898, 24006476, 24370118, 25032700, 25157968, 31752122, 29525983, 29721857], 'allele_string': 'T/A/C/G', 'clin_sig_allele': 'C:likely_pathogenic;A:pathogenic;G:likely_pathogenic;A:likely_pathogenic;C:pathogenic;G:pathogenic'}], 'transcript_consequences': [{'biotype': 'protein_coding', 'sift_prediction': 'tolerated', 'polyphen_score': 0.251, 'sift_score': 0.06, 'impact': 'MODERATE', 'variant_allele': 'C', 'amino_acids': 'Q/R', 'cds_start': 182, 'canonical': 1, 'cdna_end': 313, 'gene_symbol': 'NRAS', 'hgvsp': 'ENSP00000358548.4:p.Gln61Arg', 'codons': 'cAa/cGa', 'gene_symbol_source': 'HGNC', 'cds_end': 182, 'hgvsc': 'ENST00000369535.5:c.182A>G', 'cdna_start': 313, 'strand': -1, 'polyphen_prediction': 'benign', 'gene_id': 'ENSG00000213281', 'protein_start': 61, 'transcript_id': 'ENST00000369535', 'consequence_terms': ['missense_variant'], 'hgnc_id': 'HGNC:7989', 'protein_end': 61}], 'variant_class': 'SNV', 'input': '1 114713908 114713908 T/C 1', 'allele_string': 'T/C', 'start': 114713908, 'seq_region_name': '1'}
        self.assertEqual(target.VCF_file.ensembl_vep(query1), test_query1_result)

class TestGetAnnotationText(unittest.TestCase):
    def test_get_annotation_text(self):
        """ Test if API response is parsed into the correct annotation text."""
        json1 = {'strand': 1, 'id': '1_114713908_T/C', 'vcf_string': '1-114713908-T-C', 'assembly_name': 'GRCh38', 'end': 114713908, 'most_severe_consequence': 'missense_variant', 'colocated_variants': [{'allele_string': 'COSMIC_MUTATION', 'id': 'COSV54736340', 'start': 114713908, 'strand': 1, 'phenotype_or_disease': 1, 'somatic': 1, 'end': 114713908, 'seq_region_name': '1', 'var_synonyms': {'COSMIC': ['COSM584']}}, {'var_synonyms': {'COSMIC': ['COSM583']}, 'seq_region_name': '1', 'somatic': 1, 'phenotype_or_disease': 1, 'strand': 1, 'end': 114713908, 'start': 114713908, 'id': 'COSV54736624', 'allele_string': 'COSMIC_MUTATION'}, {'id': 'COSV54738969', 'allele_string': 'COSMIC_MUTATION', 'strand': 1, 'somatic': 1, 'phenotype_or_disease': 1, 'end': 114713908, 'start': 114713908, 'var_synonyms': {'COSMIC': ['COSM582']}, 'seq_region_name': '1'}, {'start': 114713908, 'end': 114713908, 'strand': 1, 'phenotype_or_disease': 1, 'somatic': 1, 'seq_region_name': '1', 'var_synonyms': {'COSMIC': ['COSM5044298']}, 'allele_string': 'COSMIC_MUTATION', 'id': 'COSV54747786'}, {'end': 114713908, 'strand': 1, 'phenotype_or_disease': 1, 'start': 114713908, 'clin_sig': ['likely_pathogenic', 'pathogenic'], 'var_synonyms': {'ClinVar': ['RCV000291285', 'VCV000280409', 'RCV000444899', 'RCV000444754', 'RCV000444660', 'RCV000444278', 'RCV000444188', 'RCV000441171', 'RCV000439765', 'RCV000439526', 'RCV000438738', 'RCV000438468', 'RCV000437545', 'RCV000437312', 'RCV000436539', 'RCV000435905', 'RCV000435412', 'RCV000434604', 'RCV000433761', 'RCV000433349', 'RCV000432170', 'RCV000431592', 'RCV000431333', 'RCV000430000', 'RCV000429082', 'RCV000428903', 'RCV000428055', 'RCV000427746', 'RCV000426654', 'RCV000426122', 'RCV000424084', 'RCV000423898', 'RCV000422640', 'RCV000421496', 'RCV000421291', 'RCV000420302', 'RCV000419201', 'RCV000419053', 'RCV000418396', 'RCV000418220', 'VCV000375874', 'RCV000014914', 'VCV000013900', 'RCV000037574', 'RCV000032847', 'RCV000435687', 'RCV000432961', 'RCV000420832', 'RCV000419710', 'RCV000431883', 'RCV000430593', 'RCV000430407', 'RCV000424960', 'RCV000424721', 'RCV000424455', 'RCV000422278', 'RCV000422078', 'RCV000413804', 'RCV000445249', 'RCV000439264', 'RCV000438052', 'RCV000441317', 'RCV000440367', 'RCV000114745', 'RCV000114744', 'RCV000148032'], 'UniProt': ['VAR_006847'], 'OMIM': [164790.0002]}, 'seq_region_name': '1', 'id': 'rs11554290', 'pubmed': [24033266, 26619011, 22499344, 1654209, 2278970, 2674680, 3122217, 6587382, 8120410, 12460918, 12727991, 14508525, 16273091, 16291983, 16434492, 17699718, 18390968, 18633438, 18948947, 19075190, 19880792, 20130576, 20149136, 20179705, 20406486, 20619739, 20736745, 21107323, 21305640, 21576590, 21729679, 21829508, 22407852, 22761467, 22773810, 23392294, 23414587, 23515407, 23538902, 23569304, 23614898, 24006476, 24370118, 25032700, 25157968, 31752122, 29525983, 29721857], 'allele_string': 'T/A/C/G', 'clin_sig_allele': 'C:likely_pathogenic;A:pathogenic;G:likely_pathogenic;A:likely_pathogenic;C:pathogenic;G:pathogenic'}], 'transcript_consequences': [{'biotype': 'protein_coding', 'sift_prediction': 'tolerated', 'polyphen_score': 0.251, 'sift_score': 0.06, 'impact': 'MODERATE', 'variant_allele': 'C', 'amino_acids': 'Q/R', 'cds_start': 182, 'canonical': 1, 'cdna_end': 313, 'gene_symbol': 'NRAS', 'hgvsp': 'ENSP00000358548.4:p.Gln61Arg', 'codons': 'cAa/cGa', 'gene_symbol_source': 'HGNC', 'cds_end': 182, 'hgvsc': 'ENST00000369535.5:c.182A>G', 'cdna_start': 313, 'strand': -1, 'polyphen_prediction': 'benign', 'gene_id': 'ENSG00000213281', 'protein_start': 61, 'transcript_id': 'ENST00000369535', 'consequence_terms': ['missense_variant'], 'hgnc_id': 'HGNC:7989', 'protein_end': 61}], 'variant_class': 'SNV', 'input': '1 114713908 114713908 T/C 1', 'allele_string': 'T/C', 'start': 114713908, 'seq_region_name': '1'}
        test_json1_result = 'variant_class=SNV;consequence=missense_variant;COSMIC=COSV54736340,COSV54736624,COSV54738969,COSV54747786,rs11554290;COSMIC_legacy=COSM584,COSM583,COSM582,COSM5044298;impact=MODERATE;biotype=protein_coding;amino_acids=Q/R;canonical=1;SIFT=tolerated(0.06);PolyPhen=benign(0.251)'
        self.assertEqual(target.VCF_file.get_annotation_text(json1), test_json1_result)
    