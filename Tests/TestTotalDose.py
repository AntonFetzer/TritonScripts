"""Regression test for Dependencies.TotalDose.totalDose.

This pins the numerical conventions that a hand-written GRAS parser silently
gets wrong, which is the reason to call totalDose rather than reimplement it:

* Read.ReadDose converts rad to kRad.
* The per-tile average is weighted by entries, not by file. The two fixture
  files deliberately carry different entry counts, so an unweighted mean gives
  a visibly different answer and is asserted against.
* The pooled error is sqrt(sum(N_i^2 * sigma_i^2)) / sum(N_i).
* Entries and non-zero entries are summed across files.

Run from the repository root:

    python3 -m unittest discover -s Tests -p "Test*.py"
"""
import sys
import unittest
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.TotalDose import totalDose

FIXTURE_DIR = Path(__file__).resolve().parent / "Fixtures" / "TwoFileRun"

# Fixture contents, in rad as GRAS writes them.
ENTRIES_A, ENTRIES_B = 1.0e6, 3.0e6
DOSE_A = np.array([1.0e-6, 2.0e-6, 4.0e-6])
DOSE_B = np.array([1.2e-6, 2.1e-6, 3.9e-6])
ERROR_A = np.array([2.0e-7, 2.0e-7, 2.0e-7])
ERROR_B = np.array([2.0e-7, 2.0e-7, 2.0e-7])
NON_ZEROS_A = np.array([100, 200, 300])
NON_ZEROS_B = np.array([400, 500, 600])
RAD_TO_KRAD = 1.0e-3


class TestTotalDose(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not FIXTURE_DIR.is_dir():
            raise unittest.SkipTest("fixture directory missing: {}".format(FIXTURE_DIR))
        cls.result = totalDose(str(FIXTURE_DIR))

    def test_dose_is_weighted_by_entries_and_in_krad(self):
        total_entries = ENTRIES_A + ENTRIES_B
        expected = (DOSE_A * ENTRIES_A + DOSE_B * ENTRIES_B) / total_entries * RAD_TO_KRAD
        np.testing.assert_allclose(self.result["dose"], expected, rtol=1e-12)

    def test_unweighted_mean_would_be_wrong(self):
        """Guard the specific mistake this test exists to prevent."""
        unweighted = (DOSE_A + DOSE_B) / 2.0 * RAD_TO_KRAD
        difference = np.abs(self.result["dose"] - unweighted) / self.result["dose"]
        self.assertTrue(
            np.all(difference > 1.0e-3),
            "fixtures no longer distinguish entry weighting from a plain mean",
        )

    def test_pooled_error(self):
        total_entries = ENTRIES_A + ENTRIES_B
        expected = np.sqrt(
            (ENTRIES_A * ERROR_A) ** 2 + (ENTRIES_B * ERROR_B) ** 2
        ) / total_entries * RAD_TO_KRAD
        np.testing.assert_allclose(self.result["error"], expected, rtol=1e-12)

    def test_entries_and_non_zeros_are_summed(self):
        np.testing.assert_allclose(
            self.result["entries"], np.full(3, ENTRIES_A + ENTRIES_B))
        np.testing.assert_array_equal(
            self.result["non-zeros"], NON_ZEROS_A + NON_ZEROS_B)

    def test_returns_documented_keys(self):
        for key in ("dose", "error", "entries", "non-zeros"):
            self.assertIn(key, self.result)


if __name__ == "__main__":
    unittest.main()
