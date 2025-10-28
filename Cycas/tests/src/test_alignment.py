import unittest

from pysam import AlignmentFile, IndexedReads

from cycas.src import alignment as target


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
        self.test_read_name = "5f298acd-37f2-5f16-99b9-7830acd14f82"
        self.test_read_length = 3848
        self.test_read_alignments = 13
        self.bam = IndexedReads(bam)
        self.bam.build()
        return super().setUp()

    def test_initialization(self):
        """Test the construction of the object in a rough manner."""
        x = self.testcase(read_name=self.test_read_name, bam=self.bam, filters=[])

        self.assertEqual(x.read_name, self.test_read_name)
        self.assertEqual(x.read_length_raw, self.test_read_length)
        self.assertEqual(x.filters, [])
        self.assertEqual(len(x.alignments), self.test_read_alignments)

    def test_unaligned_regions(self):
        """Test the recovery of the length of the unaligned segments in the reads"""
        x = self.testcase(read_name=self.test_read_name, bam=self.bam, filters=[])

        unaligned_regions = [193, -24, 101, 0, 0, 724, -12, 376, -24, 101, 0, 0, 0, 45]
        result = x.find_unaligned_regions()
        self.assertEqual(result, unaligned_regions)

    def test_intermediate_read_structure(self):
        """Test the alternative representation of the read structure"""
        x = self.testcase(read_name=self.test_read_name, bam=self.bam, filters=[])

        alternative_representation = "193:U,208:BB:R:BB41C:11:193:208:0,-24:O,207:BB:R:BB41C:1:377:207:0,101:U,219:BB:R:BB41C:1:685:219:0,101:I:F:chr17:7574727:904:101:0,222:BB:R:BB41C:1:1005:222:0,724:U,203:BB:F:BB41C:1:1951:203:1,-12:O,153:BB:F:BB41C:23:2142:153:1,376:U,214:BB:R:BB41C:11:2671:214:0,-24:O,208:BB:R:BB41C:1:2861:208:0,101:U,214:BB:R:BB41C:1:3170:214:0,99:I:F:chr17:7574727:3384:99:0,219:BB:R:BB41C:1:3483:219:0,101:I:F:chr17:7574727:3702:101:0,45:U"
        result = x.create_intermediate_read_structure()
        result_length = sum([int(x.split(":")[0]) for x in result.split(",")])
        self.assertEqual(result_length, self.test_read_length)
        self.assertEqual(result, alternative_representation)
