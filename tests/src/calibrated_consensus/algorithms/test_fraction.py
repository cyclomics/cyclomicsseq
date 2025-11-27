import unittest

import numpy as np

from cycas.src.calibrated_consensus.algorithms import FractionCalibratedConsensus


class TestMakeConsensus(unittest.TestCase):
    def setUp(self) -> None:
        self.model = FractionCalibratedConsensus()
        return super().setUp()

    def test_make_consensus(self):
        msa = np.array(
            [
                ["_", "A", "C", "G", "-", "A", "C", "G", "T", "_", "_", "_"],
                ["C", "A", "C", "G", "-", "A", "C", "G", "T", "_", "_", "_"],
                ["C", "C", "-", "A", "T", "A", "C", "G", "T", "_", "_", "_"],
                ["C", "A", "C", "A", "T", "A", "C", "G", "T", "A", "C", "T"],
            ]
        )

        consensus_seq, consensus_qual = self.model.make_consensus(msa)

        self.assertEqual(consensus_seq, "CACNNACGTACT")
        self.assertEqual(
            consensus_qual,
            [
                "3/3",
                "3/4",
                "3/4",
                "2/4",
                "2/4",
                "4/4",
                "4/4",
                "4/4",
                "4/4",
                "1/1",
                "1/1",
                "1/1",
            ],
        )
