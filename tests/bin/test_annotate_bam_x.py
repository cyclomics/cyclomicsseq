from unittest import mock
import unittest

from bin import annotate_bam_x as target


class TestSmokeTest(unittest.TestCase):
    def test_Some(self):
        """
        test should always pass.
        """
        self.assertEqual("foo".upper(), "FOO")


class TestAnnotateBamX(unittest.TestCase):
    def setUp(self) -> None:
        info, indexes = target.load_seqsum(
            "tests/informed/sequencing_summary_FAS12345_abcdef123.txt"
        )
        self.info = info
        self.indexes = indexes

    def test_header_to_tag(self):
        headers_expected = [
            "channel",
            "mux",
            "channel_mux",
            "start_time",
            "duration",
            "adapter_duration",
            "sequence_length_template",
            "mean_qscore_template",
            "median_template",
            "mad_template",
            "end_reason",
        ]
        tags = target.header_to_tag

        for h in headers_expected:
            self.assertTrue(tags[h].startswith("X"))

    def test_load_seqsum(self):
        expected_readnames = [
            "d204f0f2-b00e-4d76-9581-120fc7e93cd0",
            "8fe91acc-e8ed-4f30-94aa-ade8cb6882c3",
            "c3c10fd5-9f6b-4843-b975-cdcf387ddeda",
            "2e689dde-28f1-42f1-ace4-a2e8ea5a3db2",
        ]

        # misleading name, but tests that the object has the same elements
        # https://docs.python.org/3.2/library/unittest.html#unittest.TestCase.assertCountEqual
        self.assertCountEqual(list(self.info.keys())[:4], expected_readnames)
        for k, v in self.indexes.items():
            self.assertIsInstance(k, str)
            self.assertIsInstance(v, int)

    def test_load_seqsum_header_mapping(self):
        expected_mapping = {
            "filename_fastq": 0,
            "filename_fast5": 1,
            "parent_read_id": 2,
            "read_id": 3,
            "run_id": 4,
            "channel": 5,
            "mux": 6,
            "start_time": 7,
            "duration": 8,
            "num_events": 9,
            "passes_filtering": 10,
            "template_start": 11,
            "num_events_template": 12,
            "template_duration": 13,
            "sequence_length_template": 14,
            "mean_qscore_template": 15,
            "strand_score_template": 16,
            "median_template": 17,
            "mad_template": 18,
            "pore_type": 19,
            "experiment_id": 20,
            "sample_id": 21,
            "end_reason": 22,
        }

        self.assertEqual(self.indexes, expected_mapping)


if __name__ == "__main__":
    unittest.main()
