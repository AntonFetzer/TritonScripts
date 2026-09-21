import tempfile
import unittest
from pathlib import Path

import numpy as np

from Read.ReadMCC import readMCC, readScan
from Dependencies.ScanMetrics import (cellBounds, depthDoseMetrics, fallingCrossing,
                                      levelCrossing, profileMetrics, scanMetrics)

FIXTURE = Path(__file__).resolve().parent / "Fixtures" / "TwoScanProfiles.mcc"


class TestReadMCC(unittest.TestCase):
    def test_reads_every_scan(self):
        scans = readMCC(FIXTURE)
        self.assertEqual([scan["curve_type"] for scan in scans],
                         ["INPLANE_PROFILE", "CROSSPLANE_PROFILE", "PDD"])

    def test_selects_one_scan_by_curve_type(self):
        self.assertEqual(readScan(FIXTURE, "PDD")["curve_type"], "PDD")
        with self.assertRaises(ValueError):
            readScan(FIXTURE)
        with self.assertRaises(ValueError):
            readScan(FIXTURE, "DIAGONAL_PROFILE")

    def test_scan_fields(self):
        scan = readMCC(FIXTURE)[0]
        self.assertEqual(scan["depth"], 10.0)
        self.assertEqual(scan["unit"], "A.U.")
        self.assertEqual(scan["metadata"]["SCAN_CURVETYPE"], "INPLANE_PROFILE")
        np.testing.assert_allclose(scan["position"], [-20.0, 0.0, 20.0])
        np.testing.assert_allclose(scan["value"], [0.5, 1.0, 0.5])
        np.testing.assert_allclose(scan["reference"], [19.207] * 3)

    def test_two_column_scan_has_no_reference(self):
        self.assertIsNone(readMCC(FIXTURE)[1]["reference"])

    def test_rejects_file_without_scans(self):
        with tempfile.TemporaryDirectory() as folder:
            empty = Path(folder) / "Empty.mcc"
            empty.write_text("BEGIN_SCAN_DATA\nEND_SCAN_DATA\n")
            with self.assertRaises(ValueError):
                readMCC(empty)

    def test_rejects_ragged_data_block(self):
        with tempfile.TemporaryDirectory() as folder:
            ragged = Path(folder) / "Ragged.mcc"
            ragged.write_text(
                "BEGIN_SCAN 1\nSCAN_CURVETYPE=INPLANE_PROFILE\nBEGIN_DATA\n"
                "0.0 1.0 2.0\n1.0 1.0\nEND_DATA\nEND_SCAN 1\n")
            with self.assertRaises(ValueError):
                readMCC(ragged)


class TestDepthDoseMetrics(unittest.TestCase):
    def setUp(self):
        self.scan = readScan(FIXTURE, "PDD")

    def test_normalises_to_maximum(self):
        metrics = depthDoseMetrics(self.scan)
        self.assertAlmostEqual(metrics["reference_value"], 1.0)
        self.assertAlmostEqual(metrics["dmax_sampled_mm"], 10.0)
        self.assertAlmostEqual(metrics["surface_percent"], 50.0)

    def test_ranges_use_the_falling_edge_only(self):
        metrics = depthDoseMetrics(self.scan)
        self.assertAlmostEqual(metrics["ranges"][50.0], 20.0)
        self.assertAlmostEqual(metrics["ranges"][80.0], 14.0)
        self.assertAlmostEqual(metrics["tail_percent"], 10.0)
        self.assertFalse(metrics["truncated"])

    def test_falling_crossing_ignores_the_rising_edge(self):
        depth = np.array([0.0, 10.0, 20.0])
        values = np.array([50.0, 100.0, 50.0])
        self.assertAlmostEqual(fallingCrossing(depth, values, 75.0), 15.0)


class TestScanMetricsDispatch(unittest.TestCase):
    def test_dispatches_on_curve_type(self):
        kinds = [scanMetrics(scan)["kind"] for scan in readMCC(FIXTURE)]
        self.assertEqual(kinds, ["profile", "profile", "depth_dose"])

    def test_rejects_unknown_curve_type(self):
        scan = dict(readMCC(FIXTURE)[0], curve_type="DIAGONAL_UNSUPPORTED")
        with self.assertRaises(ValueError):
            scanMetrics(scan)


class TestCellBounds(unittest.TestCase):
    def test_first_boundary_is_the_surface(self):
        bounds = cellBounds(np.array([0.0, 1.5, 3.0]))
        np.testing.assert_allclose(bounds, [0.0, 0.75, 2.25, 3.75])

    def test_rejects_single_depth(self):
        with self.assertRaises(ValueError):
            cellBounds(np.array([1.0]))


class TestProfileMetrics(unittest.TestCase):
    def test_normalises_to_central_axis_not_maximum(self):
        scans = readMCC(FIXTURE)
        metrics = profileMetrics(scans[1])
        self.assertAlmostEqual(metrics["reference_value"], 0.8)
        self.assertAlmostEqual(metrics["normalized"].max(), 100.0)

    def test_symmetric_triangle_widths(self):
        metrics = profileMetrics(readMCC(FIXTURE)[0])
        self.assertAlmostEqual(metrics["widths"][50.0], 40.0)
        self.assertAlmostEqual(metrics["offsets"][50.0], 0.0)
        self.assertAlmostEqual(metrics["asymmetry_percentage_points"], 0.0)

    def test_unreached_level_is_nan_not_extrapolated(self):
        metrics = profileMetrics(readMCC(FIXTURE)[0])
        self.assertTrue(np.isnan(metrics["widths"][20.0]))
        self.assertTrue(metrics["truncated"])

    def test_level_crossing_sides(self):
        position = np.array([-20.0, 0.0, 20.0])
        normalized = np.array([50.0, 100.0, 50.0])
        self.assertAlmostEqual(levelCrossing(position, normalized, 75.0, "left"), -10.0)
        self.assertAlmostEqual(levelCrossing(position, normalized, 75.0, "right"), 10.0)


if __name__ == "__main__":
    unittest.main()
