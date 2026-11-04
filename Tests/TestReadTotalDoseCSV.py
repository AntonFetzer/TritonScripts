"""Read.ReadTotalDoseCSV: a TotalDose CSV in the dictionary format of totalDose."""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from Dependencies.MergeTotalDose import mergeTotalDose
from Read.ReadTotalDoseCSV import readTotalDoseCSV

TABLE = """Tile, Dose [kRad], Error [kRad], Non-Zero Entries
0,5.97e-12,2.1e-13,1722
1,0.0,0.0,0
2,4.77e-12,1.7e-13,1659
"""


class ReadTotalDoseCSV(unittest.TestCase):
    def read(self, text, *args):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "TotalDose_run.csv"
            path.write_text(text, encoding="ascii")
            return readTotalDoseCSV(path, *args)

    def test_values_and_keys(self):
        d = self.read(TABLE)
        self.assertEqual(sorted(d), ["dose", "entries", "error", "non-zeros"])
        np.testing.assert_array_equal(d["dose"], [5.97e-12, 0.0, 4.77e-12])
        np.testing.assert_array_equal(d["error"], [2.1e-13, 0.0, 1.7e-13])
        np.testing.assert_array_equal(d["non-zeros"], [1722, 0, 1659])
        np.testing.assert_array_equal(d["entries"], [0, 0, 0])

    def test_scale_applies_to_dose_and_error_only(self):
        d = self.read(TABLE, 1000.0)
        np.testing.assert_allclose(d["dose"], [5.97e-9, 0.0, 4.77e-9])
        np.testing.assert_allclose(d["error"], [2.1e-10, 0.0, 1.7e-10])
        np.testing.assert_array_equal(d["non-zeros"], [1722, 0, 1659])

    def test_merges_with_a_total_dose_dictionary(self):
        # Integer entries and non-zeros, as from totalDose: mergeTotalDose adds them in place.
        other = {
            "dose": np.array([1.0, 2.0, 3.0]),
            "error": np.array([3e-12, 4e-12, 0.0]),
            "entries": np.array([100, 100, 100]),
            "non-zeros": np.array([10, 20, 30]),
        }
        total = mergeTotalDose([other, self.read(TABLE)])
        np.testing.assert_allclose(total["dose"], [1.0 + 5.97e-12, 2.0, 3.0 + 4.77e-12])
        np.testing.assert_allclose(total["error"], [np.hypot(3e-12, 2.1e-13), 4e-12, 1.7e-13])
        np.testing.assert_array_equal(total["entries"], [100, 100, 100])
        np.testing.assert_array_equal(total["non-zeros"], [1732, 20, 1689])

    def test_wrong_header(self):
        with self.assertRaises(ValueError):
            self.read("Tile,Dose\n0,1.0\n")

    def test_tiles_out_of_order(self):
        with self.assertRaises(ValueError):
            self.read("Tile, Dose [kRad], Error [kRad], Non-Zero Entries\n1,1.0,0.1,5\n0,1.0,0.1,5\n")


if __name__ == "__main__":
    unittest.main()
