"""Regression test for Dependencies.HistogramPlots.

The drawing helpers replaced six hand-rolled copies of the same idiom, so what
is pinned here is the part that silently went wrong in the copies: drawing a
histogram whose arrays do not line up used to raise deep inside matplotlib or
produce a misaligned plot. checkHistogram must refuse it up front instead.

Run from the repository root:

    python3 -m unittest discover -s Tests -p "Test*.py"
"""
import sys
import unittest
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.HistogramPlots import checkHistogram, nonZeroRange


def makeHistogram(bins=4):
    """A minimal histogram obeying the repository's reader contract."""
    lower = np.arange(bins, dtype=float)
    return {
        "lower": lower,
        "upper": lower + 1.0,
        "mean": lower + 0.5,
        "value": np.linspace(1.0, 2.0, bins),
        "error": np.full(bins, 0.1),
        "entries": np.arange(10, 10 + bins, dtype=float),
    }


class TestCheckHistogram(unittest.TestCase):

    def test_rejects_mismatched_lengths(self):
        """The failure the copied code turned into a misaligned plot."""
        histogram = makeHistogram()
        histogram["value"] = histogram["value"][:-1]
        with self.assertRaises(ValueError):
            checkHistogram(histogram)


class TestNonZeroRange(unittest.TestCase):
    def test_spans_the_filled_bins_only(self):
        """Empty bins at both ends are valid data but must not set the limits."""
        histogram = makeHistogram(bins=6)
        histogram["value"] = np.array([0.0, 0.0, 3.0, 0.0, 5.0, 0.0])
        self.assertEqual(nonZeroRange([histogram]), (2.0, 5.0))

    def test_all_empty_returns_none(self):
        histogram = makeHistogram()
        histogram["value"] = np.zeros(4)
        self.assertIsNone(nonZeroRange([histogram]))


if __name__ == "__main__":
    unittest.main()
