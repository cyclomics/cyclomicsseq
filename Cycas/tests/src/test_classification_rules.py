import unittest
from unittest.mock import MagicMock, Mock

from cycas.src import classification_rules as target


class TestSmokeTest(unittest.TestCase):
    def test_Some(self):
        """
        test should always pass.
        """
        self.assertEqual("foo".upper(), "FOO")


class TestCheckMinimumUnalignedGap(unittest.TestCase):
    def setUp(self) -> None:
        self.testcase = target.CheckMinimumUnalignedGap
        self.mock_group = Mock()
        self.mock_group.find_unaligned_regions = MagicMock(
            return_value=[124, 122, 25, 72, 2]
        )
        return super().setUp()

    def test_score_higher_target(self):
        """Test value below target"""
        this_case = self.testcase(target_value=123)

        self.assertEqual(this_case.score(self.mock_group), 1)

    def test_score_lower_target(self):
        """Test if smallest gap is bigger than the target value"""
        this_case = self.testcase(target_value=24)

        self.assertEqual(this_case.score(self.mock_group), 0)

    def test_score_low_alignment_count(self):
        self.mock_group.find_unaligned_regions = MagicMock(return_value=[124, 2])
        this_case = self.testcase(target_value=24)

        self.assertEqual(this_case.score(self.mock_group), 1)

    def test_score_empty_list(self):
        self.mock_group.find_unaligned_regions = MagicMock(return_value=[])
        this_case = self.testcase(target_value=24)

        self.assertEqual(this_case.score(self.mock_group), 1)
