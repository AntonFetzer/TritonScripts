import tempfile
import unittest
from pathlib import Path

import numpy as np

from Read.ReadMCC import readMCC, readScan
from Dependencies.ScanMetrics import depthDoseMetrics, profileMetrics, scanMetrics

FIXTURE = Path(__file__).resolve().parent / "Fixtures" / "TwoScanProfiles.mcc"


class TestReadMCC(unittest.TestCase):
    def test_reads_every_scan(self):
        scans = readMCC(FIXTURE)
        self.assertEqual([scan["curve_type"] for scan in scans],
                         ["INPLANE_PROFILE", "CROSSPLANE_PROFILE", "PDD"])

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


class TestScanMetricsDispatch(unittest.TestCase):
    def test_dispatches_on_curve_type(self):
        kinds = [scanMetrics(scan)["kind"] for scan in readMCC(FIXTURE)]
        self.assertEqual(kinds, ["profile", "profile", "depth_dose"])


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


if __name__ == "__main__":
    unittest.main()
