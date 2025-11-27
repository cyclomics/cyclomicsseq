import unittest

from cycas.src.app import _create_fastq_entry as target


class TestCreateFastqEntry(unittest.TestCase):
    def setUp(self) -> None:
        self.read_name = "test-read_0_3D"
        self.sequence = "NNATGCATGCNN"
        self.quality = [0, 0, 6, 7, 9, 45, 59, 77, 9, 6, 0, 0]
        self.quality_str = r"!!'(*N\n*'!!"
        return super().setUp()

    def test_create_fastq_entry(self):
        expected = f"@test-read_0_3D\nNNATGCATGCNN\n+\n{self.quality_str}\n"
        result = target(self.read_name, self.sequence, self.quality)
        print(result)
        self.assertEqual(result, expected)
