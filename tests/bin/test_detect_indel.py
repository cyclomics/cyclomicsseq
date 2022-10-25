import unittest
from collections import defaultdict
from math import exp
from unittest.mock import MagicMock, Mock

from bin import determine_vaf as target


class TestSmokeTest(unittest.TestCase):
    def test_Some(self):
        """
        test should always pass.
        """
        self.assertEqual("foo".upper(), "FOO")


class TestDetermineVAF(unittest.TestCase):
    def setUp(self) -> None:
        info, indexes = target.load_seqsum()
        self.info = info
        self.indexes = indexes
