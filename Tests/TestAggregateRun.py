"""Regression test for Dependencies.AggregateRun.

aggregateRun exists to stop a driver reporting a plausible number from a run
that is incomplete or does not match its tally list, so what is pinned here is
that each guard actually fires, and that the derived quantities it adds are
consistent with the arrays totalDose returns.

Reuses the TwoFileRun fixture: two files, three tiles, deliberately different
entry counts. See TestTotalDose.py for the numerical conventions themselves.

Run from the repository root:

    python3 -m unittest discover -s Tests -p "Test*.py"
"""
import sys
import unittest
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.AggregateRun import aggregateRun, tileRecords

FIXTURE_DIR = Path(__file__).resolve().parent / "Fixtures" / "TwoFileRun"
FIXTURE_FILES = 2
FIXTURE_TILES = 3


class TestAggregateRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not FIXTURE_DIR.is_dir():
            raise unittest.SkipTest(f"fixture directory missing: {FIXTURE_DIR}")
        cls.result = aggregateRun(
            FIXTURE_DIR, expectedFiles=FIXTURE_FILES, tileCount=FIXTURE_TILES)

    def test_passes_totalDose_keys_through(self):
        for key in ("dose", "error", "entries", "non-zeros"):
            self.assertIn(key, self.result)
            self.assertEqual(len(self.result[key]), FIXTURE_TILES)

    def test_reports_the_file_count(self):
        self.assertEqual(self.result["files"], FIXTURE_FILES)

    def test_relative_error_matches_the_arrays(self):
        expected = 100.0 * np.asarray(self.result["error"]) / np.asarray(
            self.result["dose"])
        np.testing.assert_allclose(
            self.result["relative_error_percent"], expected, rtol=1e-12)

    def test_hit_fraction_matches_the_arrays(self):
        expected = np.asarray(self.result["non-zeros"], dtype=float) / np.asarray(
            self.result["entries"], dtype=float)
        np.testing.assert_allclose(
            self.result["hit_fraction"], expected, rtol=1e-12)

    def test_wrong_file_count_is_refused(self):
        """An incomplete run must not aggregate silently."""
        with self.assertRaises(RuntimeError):
            aggregateRun(FIXTURE_DIR, expectedFiles=FIXTURE_FILES + 1)

    def test_wrong_tile_count_is_refused(self):
        """A driver whose tally list disagrees with the detector macro must fail."""
        with self.assertRaises(RuntimeError):
            aggregateRun(FIXTURE_DIR, tileCount=FIXTURE_TILES + 1)

    def test_missing_directory_is_refused(self):
        with self.assertRaises(FileNotFoundError):
            aggregateRun(FIXTURE_DIR / "does-not-exist")

    def test_tileRecords_follows_the_given_order(self):
        """Tally order comes from the detector macro, not from the GRAS output."""
        tiles = [(2, "C_PV", "third"), (0, "A_PV", "first")]
        records = tileRecords(self.result, tiles)
        self.assertEqual([record["volume"] for record in records],
                         ["C_PV", "A_PV"])
        self.assertEqual([record["index"] for record in records], [2, 0])
        self.assertAlmostEqual(records[0]["dose"], float(self.result["dose"][2]))
        self.assertAlmostEqual(records[1]["dose"], float(self.result["dose"][0]))


if __name__ == "__main__":
    unittest.main()
