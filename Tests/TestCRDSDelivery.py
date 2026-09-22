"""Regression test for Dependencies.CRDSDelivery.

The delivered value of a part is a sum over exposures of a coefficient times
the local fluence. What is pinned here is the error model, because it is the
part a hand calculation gets wrong: a coefficient reused by two exposures is
ONE measurement, so its error scales with the summed fluence rather than
adding in quadrature, while coefficients from different runs are independent.
The map sampling and the board-shift range are pinned on a synthetic linear
map written to a temporary directory, where bilinear interpolation is exact.

Run from the repository root:

    python3 -m unittest discover -s Tests -p "Test*.py"
"""
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.CRDSDelivery import combine, deliveredRecords, sampleFluence


class TestCombine(unittest.TestCase):
    def test_one_run_for_both_exposures_is_fully_correlated(self):
        delivered, error = combine({"wide": (2.0, 0.1)}, ["wide", "wide"], [3.0, 5.0])
        self.assertAlmostEqual(delivered, 16.0)
        self.assertAlmostEqual(error, 0.1 * 8.0)

    def test_different_runs_add_in_quadrature(self):
        delivered, error = combine({"narrow": (2.0, 0.1), "wide": (4.0, 0.3)},
                                   ["narrow", "wide"], [3.0, 5.0])
        self.assertAlmostEqual(delivered, 2.0 * 3.0 + 4.0 * 5.0)
        self.assertAlmostEqual(error, math.hypot(0.1 * 3.0, 0.3 * 5.0))

    def test_mismatched_lengths_are_refused(self):
        with self.assertRaises(ValueError):
            combine({"wide": (1.0, 0.1)}, ["wide"], [1.0, 2.0])


class TestDeliveredRecords(unittest.TestCase):
    """Two exposures on linear maps phi = a + b*y, sampled at a part at y = 10."""

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        root = Path(self.directory.name)
        edges = np.arange(-50.0, 51.0, 1.0)
        centres = (edges[1:] + edges[:-1]) / 2
        # values[row, column], rows along y as the accelerator-log maps store them
        yy = np.repeat(centres[:, None], len(centres), axis=1)
        for name, (a, b, average) in {"narrow": (100.0, -2.0, 90.0),
                                      "wide": (50.0, 0.0, 50.0)}.items():
            (root / name).mkdir()
            np.savez(root / name / "heatmap_isocentre.npz",
                     fluence_protons_cm2=a + b * yy, x_edges_mm=edges, y_edges_mm=edges)
            (root / name / "summary.json").write_text(
                json.dumps({"crds_mean_fluence_raw_protons_cm2": average}))
        self.delivery = {
            "maps": root,
            "parts": {"LED": (0.0, 10.0)},
            "partOfQuantity": {"DDD": "LED"},
            "alignmentShift_mm": 1.0,
            "exposures": [
                {"label": "first", "map": "narrow", "coefficientRun": {"LED": "wide"}},
                {"label": "second", "map": "wide", "coefficientRun": {"LED": "wide"}},
            ],
        }

    def tearDown(self):
        self.directory.cleanup()

    def test_linear_map_is_sampled_exactly(self):
        path = Path(self.directory.name) / "narrow" / "heatmap_isocentre.npz"
        self.assertAlmostEqual(sampleFluence(path, 0.0, 10.0), 80.0)
        self.assertAlmostEqual(sampleFluence(path, 3.3, 10.5), 79.0)

    def test_delivered_uses_local_fluence_and_reused_coefficient(self):
        record, = deliveredRecords(self.delivery, {("DDD", "active"): {"wide": (2.0, 0.1)}})
        self.assertAlmostEqual(record["delivered"], 2.0 * (80.0 + 50.0))
        self.assertAlmostEqual(record["error"], 0.1 * (80.0 + 50.0))
        self.assertAlmostEqual(record["fieldAverageDelivered"], 2.0 * (90.0 + 50.0))
        self.assertEqual([e["localFluence"] for e in record["exposures"]], [80.0, 50.0])

    def test_alignment_range_follows_the_gradient(self):
        record, = deliveredRecords(self.delivery, {("DDD", "active"): {"wide": (2.0, 0.1)}})
        # a 1 mm shift in y moves the narrow map by 2 per mm; x does nothing
        self.assertAlmostEqual(record["alignmentLow"], 2.0 * (78.0 + 50.0))
        self.assertAlmostEqual(record["alignmentHigh"], 2.0 * (82.0 + 50.0))


if __name__ == "__main__":
    unittest.main()
