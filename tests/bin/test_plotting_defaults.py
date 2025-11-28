import unittest


class TestSmokeTest(unittest.TestCase):
    def test_Some(self):
        """
        Test should always pass.
        """
        self.assertEqual("foo".upper(), "FOO")


class TestNextflowParamsParser(unittest.TestCase):
    def setUp(self) -> None:
        self.test_params1 = "[a:b, c:[d:e, f:g, h:[i:j], h2:[i2:[nested:3]]], k:l,m:n]"
        self.test_params2 = """[reference:/home/dami/Data/references/Homo_sapiens/T2T/chm13v2.0.fa, 
backbone_fasta:/home/dami/Data/backbones/backbones_db_40.fasta, 
input_read_dir:/media/dami/cyclomics_003/raw_data/USEQ/DER4820, 
output_dir:testrun6, 
user_conda_location:/home/dami/Software/cycloseq/environment.yml, 
economy_mode:null, profile_selected:conda, 
cycas_location:/home/dami/Software/cycloseq/Cycas/cycas/cycas.py, 
filtering:[minimun_raw_length:1000], 
tidehunter:[headerlines:readName	repN	copyNumbareadLenestartgneendegmenconsLen	aveMatch	fullLen	subPos	consSeq, headerlinesQual:readName	repN	copyNum	readLen	start	end	consLen	aveMatch	fullLen	subPos	consSeq	quality, minimum_match_ratio:0.4, minimum_period:60, minimum_copy:3, kmer_length:16, primer_length:10],
varscan:[min_support_reads:10, min_var_freq:0.001, min_indel_freq:0.01, min_mq:20, max_depth:100000000, min_coverge:1000, min_avg_qual:15, min_bqual:15],
 minimap2parameterized:[min_chain_score:1, min_chain_count:10, min_peak_aln_score:20], 
 medaka:[depth:1, model:r104_e81_sup_g5015, chunk_len:100, chunk_ovlp:50, method:spoa, length:50, batch_size:1000], 
 read_pattern:{pass,fastq_pass}/**.{fq,fastq,fq.gz,fastq.gz}, sequencing_quality_summary:sequencing_summary*.txt, backbone:BB41, backbone_name:, region_file:auto, qc:simple, consensus_calling:cycas, alignment:minimap, variant_calling:validate, extra_haplotyping:skip, report:yes, quick_results:false, min_repeat_count:3]"""

    def test_params_parser_simple(self):
        """Test params parser with dummy example."""
        test_params1_result = {
            "a": "b",
            "c": {"d": "e", "f": "g", "h": {"i": "j"}, "h2": {"i2": {"nested": "3"}}},
            "k": "l",
            "m": "n",
        }
        # self.assertEqual(
        #     target.nextflow_params_parser(self.test_params1), test_params1_result
        # )

    def test_params_parser_real(self):
        """Test params parser with a real example."""
        test_params2_result = {
            "reference": "/home/dami/Data/references/Homo_sapiens/T2T/chm13v2.0.fa",
            "backbone_fasta": "/home/dami/Data/backbones/backbones_db_40.fasta",
            "input_read_dir": "/media/dami/cyclomics_003/raw_data/USEQ/DER4820",
            "output_dir": "testrun6",
            "user_conda_location": "/home/dami/Software/cycloseq/environment.yml",
            "economy_mode": "null",
            "profile_selected": "conda",
            "cycas_location": "/home/dami/Software/cycloseq/Cycas/cycas/cycas.py",
            "filtering": {"minimun_raw_length": "1000"},
            "tidehunter": {
                "headerlines": "readName\trepN\tcopyNumbareadLenestartgneendegmenconsLen\taveMatch\tfullLen\tsubPos\tconsSeq",
                "headerlinesQual": "readName\trepN\tcopyNum\treadLen\tstart\tend\tconsLen\taveMatch\tfullLen\tsubPos\tconsSeq\tquality",
                "minimum_match_ratio": "0.4",
                "minimum_period": "60",
                "minimum_copy": "3",
                "kmer_length": "16",
                "primer_length": "10",
            },
            "varscan": {
                "min_support_reads": "10",
                "min_var_freq": "0.001",
                "min_indel_freq": "0.01",
                "min_mq": "20",
                "max_depth": "100000000",
                "min_coverge": "1000",
                "min_avg_qual": "15",
                "min_bqual": "15",
            },
            "minimap2parameterized": {
                "min_chain_score": "1",
                "min_chain_count": "10",
                "min_peak_aln_score": "20",
            },
            "medaka": {
                "depth": "1",
                "model": "r104_e81_sup_g5015",
                "chunk_len": "100",
                "chunk_ovlp": "50",
                "method": "spoa",
                "length": "50",
                "batch_size": "1000",
            },
            "read_pattern": "{pass,fastq_pass}/**.{fq,fastq,fq.gz,fastq.gz}",
            "sequencing_quality_summary": "sequencing_summary*.txt",
            "backbone": "BB41",
            "backbone_name": "",
            "region_file": "auto",
            "qc": "simple",
            "consensus_calling": "cycas",
            "alignment": "minimap",
            "variant_calling": "validate",
            "extra_haplotyping": "skip",
            "report": "yes",
            "quick_results": "false",
            "min_repeat_count": "3",
        }
        # self.assertEqual(
        #     target.nextflow_params_parser(self.test_params2), test_params2_result
        # )
