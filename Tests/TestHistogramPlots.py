"""Regression test for Dependencies.HistogramPlots.

The drawing helpers replaced six hand-rolled copies of the same idiom, so what
is pinned here is the part that silently went wrong in the copies: drawing a
histogram whose arrays do not line up, or which is missing a field, used to
raise deep inside matplotlib or produce a misaligned plot. checkHistogram must
refuse it up front instead.

Run from the repository root:

    python3 -m unittest discover -s Tests -p "Test*.py"
"""
import sys
import unittest
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # never open a window from the test suite
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.HistogramPlots import (
    REQUIRED_FIELDS, checkHistogram, colour, errorbarHistogram, logLogGrid,
    savePdf, stepHistogram)


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
    def test_accepts_a_valid_histogram(self):
        histogram = makeHistogram()
        self.assertIs(checkHistogram(histogram), histogram)

    def test_rejects_a_non_dict(self):
        with self.assertRaises(TypeError):
            checkHistogram([1, 2, 3])

    def test_rejects_a_missing_required_field(self):
        for field in REQUIRED_FIELDS:
            histogram = makeHistogram()
            del histogram[field]
            with self.assertRaises(KeyError):
                checkHistogram(histogram)

    def test_rejects_a_missing_requested_field(self):
        with self.assertRaises(KeyError):
            checkHistogram(makeHistogram(), field="weights")

    def test_rejects_mismatched_lengths(self):
        """The failure the copied code turned into a misaligned plot."""
        histogram = makeHistogram()
        histogram["value"] = histogram["value"][:-1]
        with self.assertRaises(ValueError):
            checkHistogram(histogram)

    def test_rejects_non_finite_edges(self):
        histogram = makeHistogram()
        histogram["upper"] = np.array([1.0, 2.0, np.inf, 4.0])
        with self.assertRaises(ValueError):
            checkHistogram(histogram)


class TestDrawing(unittest.TestCase):
    def setUp(self):
        self.figure = plt.figure()

    def tearDown(self):
        plt.close(self.figure)

    def test_step_draws_a_bar_per_bin_and_one_outline(self):
        axes = stepHistogram(makeHistogram(bins=4), label="run A")
        self.assertEqual(len(axes.containers), 1)
        self.assertEqual(len(axes.containers[0]), 4)
        self.assertEqual(axes.get_legend_handles_labels()[1], ["run A"])

    def test_errorbar_can_omit_the_fill(self):
        axes = errorbarHistogram(makeHistogram(), label="run B", fill=False)
        self.assertTrue(any(isinstance(container, matplotlib.container.ErrorbarContainer)
                            for container in axes.containers))

    def test_entries_field_is_drawable(self):
        histogram = makeHistogram()
        axes = stepHistogram(histogram, field="entries")
        heights = [patch.get_height() for patch in axes.containers[0]]
        np.testing.assert_allclose(heights, histogram["entries"])

    def test_loglog_grid_sets_scales_without_an_empty_legend(self):
        """Unlabelled series must not produce an empty legend box."""
        stepHistogram(makeHistogram(), label=None)
        axes = logLogGrid(title="t", xlabel="x", ylabel="y")
        self.assertEqual(axes.get_xscale(), "log")
        self.assertEqual(axes.get_yscale(), "log")
        self.assertIsNone(axes.get_legend())

    def test_linear_scale_override(self):
        """Relative-error axes are linear from zero by repository convention."""
        stepHistogram(makeHistogram())
        axes = logLogGrid(yscale="linear", xscale="linear")
        self.assertEqual(axes.get_yscale(), "linear")

    def test_colour_cycles(self):
        self.assertEqual(colour(0), colour(10))
        self.assertNotEqual(colour(0), colour(1))

    def test_savePdf_writes_a_file(self):
        import tempfile
        stepHistogram(makeHistogram())
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "figure.pdf"
            self.assertEqual(savePdf(target), str(target))
            self.assertGreater(target.stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
