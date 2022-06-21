import unittest
from collections import defaultdict
from math import exp
from unittest.mock import MagicMock, Mock

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
        """
        Test that we get valid tags for occuring tags and that we dont get errors for non occuring tags.
        """
        headers_expected = [
            "channel",
            "mux",
            "start_time",
            "duration",
            "adapter_duration",
            "sequence_length_template",
            "mean_qscore_template",
            "median_template",
            "mad_template",
            "end_reason",
        ]
        tags = target.HEADER_TO_TAG

        for h in headers_expected:
            self.assertTrue(tags[h].startswith("X"))
            self.assertTrue(len(tags[h]) == 2)

        # also test that we dont get a keyerror
        self.assertEqual(
            tags["somethingthatwillneveroccurinanormalsequencingsummaryacbdef"], -1
        )

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

        def test_tag_aln_simle(self):
            """
            Test if we can add tags to an object given the required inputs.
            """
            # define testing inputs
            test_read_name = "readname123"
            fake_seqsum = {
                test_read_name: [12, "3m 12s", "we dont need this", "good_test"]
            }
            fake_seqsum_indexes = defaultdict(lambda: -1)

            fake_seqsum_indexes_entries = {"channel": 0, "duration": 1, "end_reason": 3}
            for k, v in fake_seqsum_indexes_entries.items():
                fake_seqsum_indexes[k] = v

            test_items = list(fake_seqsum_indexes.keys())

            #  expected output
            expected_tags = [
                ("already_occuring_tag", "something usefull"),
                ("XC", 12),
                ("XD", "3m 12s"),
                ("XE", "good_test"),
            ]
            # create mocked alignmet object
            mock_aln = Mock()
            mock_aln.tags = [("already_occuring_tag", "something usefull")]

            target.tag_aln(
                aln=mock_aln,
                query_name=test_read_name,
                seqsum=fake_seqsum,
                seqsum_indexes=fake_seqsum_indexes,
                items=test_items,
            )
            # Test if we only get the colums that we expect: initial plus the three that where both in seqsum and columns
            self.assertCountEqual(expected_tags, mock_aln.tags)

    def test_tag_aln_with_extra_entries(self):
        """
        Test if we can add tags to an object given the required inputs, where there is discrepancy between columns expected and presented
        """
        # define testing inputs
        test_read_name = "readname123"
        fake_seqsum = {test_read_name: [12, "3m 12s", "we dont need this", "good_test"]}
        fake_seqsum_indexes = defaultdict(lambda: -1)

        fake_seqsum_indexes_entries = {"channel": 0, "duration": 1, "end_reason": 3}
        for k, v in fake_seqsum_indexes_entries.items():
            fake_seqsum_indexes[k] = v

        test_items = list(fake_seqsum_indexes.keys())

        # add a column we dont tag
        test_items.append("something_we_dont_have_a_tag_for")
        # add a tag that is missing
        target.HEADER_TO_TAG["test_column"] = "test"
        test_items.append("test_column")

        #  expected output
        expected_tags = [
            ("already_occuring_tag", "something usefull"),
            ("XC", 12),
            ("XD", "3m 12s"),
            ("XE", "good_test"),
        ]
        # create mocked alignmet object
        mock_aln = Mock()
        mock_aln.tags = [("already_occuring_tag", "something usefull")]

        target.tag_aln(
            aln=mock_aln,
            query_name=test_read_name,
            seqsum=fake_seqsum,
            seqsum_indexes=fake_seqsum_indexes,
            items=test_items,
        )
        # Test if we only get the colums that we expect: initial plus the three that where both in seqsum and columns
        self.assertCountEqual(expected_tags, mock_aln.tags)

    def test_process_query_name(self):
        """
        Test query name splitting
        """
        input = "abc_def"
        expected_output = "abc"

        output = target.process_query_name(input)

        self.assertEqual(expected_output, output)

        input_2 = "abc_def"
        output2 = target.process_query_name(input_2, split_queryname=False)
        self.assertEqual(input_2, output2)

        input_3 = "abc_def"
        output3 = target.process_query_name(input_3, splitter="-")
        self.assertEqual(input_3, output3)

        input_4 = "abc-def-ghi"
        expected_output4 = "abc"
        output4 = target.process_query_name(input_4, splitter="-")
        self.assertEqual(expected_output4, output4)


if __name__ == "__main__":
    unittest.main()
