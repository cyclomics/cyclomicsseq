import unittest
from collections import defaultdict
from math import exp
from unittest.mock import MagicMock, Mock

from cycas.src import alignment as target

from pysam import AlignmentFile, IndexedReads


class TestSmokeTest(unittest.TestCase):
    def test_Some(self):
        """
        test should always pass.
        """
        self.assertEqual("foo".upper(), "FOO")


class TestAlignmentGroup(unittest.TestCase):
    def setUp(self) -> None:
        self.testcase = target.AlignmentGroup
        bam = AlignmentFile(
            "tests/data/FAU48563_pass_f5d283ba_42_filtered_split.NM_50_mapq_20.bam",
            "rb",
        )
        self.test_read_name = "bd257735-62dd-40bd-a302-9b348e32fe0a"
        self.bam = IndexedReads(bam)
        self.bam.build()
        return super().setUp()

    def test_initialization(self):
        """Test the construction of the object in a rough manner."""
        x = self.testcase(read_name=self.test_read_name, bam=self.bam, filters=[])

        self.assertEqual(x.read_name, self.test_read_name)
        self.assertEqual(x.read_length_raw, 3335)
        self.assertEqual(x.filters, [])
        self.assertEqual(len(x.alignments), 11)

    def test_unaligned_regions(self):
        """Test the recovery of the length of the unaligned segments in the reads"""
        x = self.testcase(read_name=self.test_read_name, bam=self.bam, filters=[])

        unaligned_regions = [223, 99, 101, 0, 0, 99, 110, 100, 115, 97, 101, 13]
        result = x.find_unaligned_regions()
        self.assertEqual(result, unaligned_regions)
