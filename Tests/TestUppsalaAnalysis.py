import csv
import json
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path

import numpy as np

from UppsalaAcceleratorLogs import uppsala_analysis as ua


STREAM_HEADER = [
    "SUBMAP_NUMBER", "TIME", "X_WIDTH(mm)", "Y_WIDTH(mm)", "X_POSITION(mm)",
    "Y_POSITION(mm)", "X_DOSE(C)", "Y_DOSE(C)", "X_DOSE_RATE(A)",
    "Y_DOSE_RATE(A)", "X_WIDTH_IC1(mm)", "Y_WIDTH_IC1(mm)",
    "X_POSITION_IC1(mm)", "Y_POSITION_IC1(mm)", "DOSE_IC1_X(C)",
    "DOSE_IC1_Y(C)", "DOSE_RATE_IC1_X(A)", "DOSE_RATE_IC1_Y(A)",
    "X_CURRENT_PRIM(V)", "Y_CURRENT_PRIM(V)", "X_VOLTAGE_PRIM(V)",
    "Y_VOLTAGE_PRIM(V)", "X_CURRENT_SEC(V)", "Y_CURRENT_SEC(V)",
    "X_VOLTAGE_SEC(V)", "Y_VOLTAGE_SEC(V)", "DOSE_PRIM(C)", "DOSE_SEC(C)",
    "DOSE_RATE_PRIM(A)", "DOSE_RATE_SEC(A)", "BEAMCURRENT(V)", "X_FIELD(G)",
    "Y_FIELD(G)", "QUAD_PRIMARY_CURRENT_X(A)", "QUAD_PRIMARY_CURRENT_Y(A)",
    "QUAD_PRIMARY_VOLTAGE_X(V)", "QUAD_PRIMARY_VOLTAGE_Y(V)",
    "QUAD_REDUNDANT_CURRENT_X(A)", "QUAD_REDUNDANT_CURRENT_Y(A)",
    "QUAD_REDUNDANT_VOLTAGE_X(V)", "QUAD_REDUNDANT_VOLTAGE_Y(V)",
]


class ParserTests(unittest.TestCase):
    def test_record_stream_and_heatmap_conservation(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sample.map_record_001_part_01.csv"
            rows = [
                ["#Version:11.0.0"],
                ["beamline id = 4"],
                ["layer id = 1"],
                ["acquisition period = 250"],
                [],
                STREAM_HEADER,
            ]
            for idx, charge in enumerate((1e-11, 2e-11, -1e-12, -10000.0)):
                row = ["0"] * 41
                row[0] = str(idx + 1)
                row[1] = f"2025-06-28T12:00:0{idx}.000+02:00"
                row[2], row[3] = "5.0", "6.0"
                row[4], row[5] = "1.0", "-2.0"
                row[6], row[7] = str(charge), str(charge)
                row[26], row[27] = str(charge), str(charge * 1.01)
                rows.append(row)
            rows.append(["2025-06-28T12:01:00.000+02:00", "Layer completed"])
            with path.open("w", newline="", encoding="utf-8") as handle:
                csv.writer(handle).writerows(rows)
            entry = {
                "id": "sample",
                "record_path": str(path),
                "spec_path": None,
                "size_bytes": path.stat().st_size,
                "session": "test",
                "kind": "part",
            }
            grid = {
                "minimum_mm": -20.0,
                "maximum_mm": 20.0,
                "pixel_mm": 1.0,
                "width_quantization_mm": 0.25,
                "width_interpretation": "sigma",
            }
            summary, arrays = ua.parse_record(entry, grid)
            self.assertEqual(summary["acquisition_rows"], 4)
            self.assertEqual(len(arrays["seconds"]), 3)
            self.assertEqual(summary["sentinel_charge_counts"]["primary"], 1)
            self.assertAlmostEqual(summary["charge_primary_raw_c"], 2.9e-11, places=20)
            self.assertAlmostEqual(summary["charge_primary_positive_c"], 3e-11, places=20)
            self.assertAlmostEqual(float(arrays["heatmap_charge_c"].sum()), 3e-11, places=17)
            self.assertLess(abs(summary["heatmap_conservation_relative_error"]), 1e-12)

    def test_repeated_stream_retains_only_last_section(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "repeated.map_record_001_part_01.csv"
            rows = [["#Version:11.0.0"], STREAM_HEADER]
            first = ["0"] * 41
            first[0] = "1"
            first[1] = "2025-06-28T12:00:00.000+02:00"
            first[2:8] = ["5", "6", "1", "-2", "1e-11", "1e-11"]
            first[26], first[27] = "1e-11", "1.01e-11"
            rows.append(first)
            rows.extend([["2025-06-28T12:00:01+02:00", "Failure"], STREAM_HEADER])
            second = first.copy()
            second[26], second[27] = "2e-11", "2.02e-11"
            rows.append(second)
            with path.open("w", newline="", encoding="utf-8") as handle:
                csv.writer(handle).writerows(rows)
            entry = {
                "id": "repeated",
                "record_path": str(path),
                "spec_path": None,
                "size_bytes": path.stat().st_size,
                "session": "test",
                "kind": "part",
            }
            grid = {
                "minimum_mm": -20.0,
                "maximum_mm": 20.0,
                "pixel_mm": 1.0,
                "width_quantization_mm": 0.25,
                "width_interpretation": "sigma",
            }
            summary, arrays = ua.parse_record(entry, grid)
            self.assertEqual(summary["stream_section_count"], 2)
            self.assertEqual(summary["ignored_repeated_stream_sections"], 1)
            self.assertEqual(summary["acquisition_rows"], 1)
            self.assertAlmostEqual(summary["charge_primary_raw_c"], 2e-11)
            self.assertAlmostEqual(arrays["series"][0, 0], 2e-11)

    def test_conversion_interpolation(self):
        conversion = [
            {"energy_mev": 60.0, "protons_per_mu": 5.40e7},
            {"energy_mev": 65.0, "protons_per_mu": 5.71e7},
        ]
        self.assertAlmostEqual(ua.interpolate_ppmu(64.0, conversion), 5.648e7)

    def test_facility_64_mev_uses_65_mev_conversion_row(self):
        conversion = [
            {"energy_mev": 60.0, "protons_per_mu": 5.40e7},
            {"energy_mev": 65.0, "protons_per_mu": 5.71e7},
        ]
        ppmu, calibration_energy = ua.facility_ppmu(64.0, conversion)
        self.assertEqual(calibration_energy, 65.0)
        self.assertEqual(ppmu, 5.71e7)

    def test_monitor_mu_uses_k_factor_not_map_specification_sum(self):
        k_factor = 0.9784840376859955
        prescribed_mu = 50000.0
        delivered_charge = (
            prescribed_mu * ua.MONITOR_CHARGE_C_PER_MU * k_factor
        )
        self.assertAlmostEqual(
            ua.monitor_mu_from_charge(delivered_charge, k_factor),
            prescribed_mu,
        )
        # Canceled/restarted specifications must not enter the conversion.
        repeated_nominal_spec_charge = delivered_charge * 3
        self.assertNotAlmostEqual(
            delivered_charge / repeated_nominal_spec_charge * prescribed_mu,
            prescribed_mu,
        )

    def test_monitor_protons_per_c(self):
        k_factor = 0.9932203389830508
        ppmu = 5.71e7
        charge = 25000 * ua.MONITOR_CHARGE_C_PER_MU * k_factor
        protons = charge * ua.monitor_protons_per_c(ppmu, k_factor)
        self.assertAlmostEqual(protons, 25000 * ppmu)

    def test_occupied_axis_limits(self):
        values = np.zeros((20, 20), dtype=float)
        values[7:10, 4:12] = 10.0
        values[0, 0] = 1e-10  # numerical tail below the relative threshold
        grid = {
            "minimum_mm": -10.0,
            "maximum_mm": 10.0,
            "pixel_mm": 1.0,
        }
        x_limits, y_limits = ua.occupied_axis_limits(
            values, grid, relative_threshold=1e-6, padding_mm=2.0
        )
        self.assertEqual(x_limits, (-8.0, 4.0))
        self.assertEqual(y_limits, (-5.0, 2.0))

    def test_radex_area(self):
        self.assertAlmostEqual(ua.RADEX_INSTRUMENT_AREA_CM2, 102.872)

    def test_radex_time_rows_use_fixed_instrument_area(self):
        products = {
            "exp_a": {"time_rows": {10: {
                "protons_delivered_raw": 1028.72,
                "protons_delivered_positive": 925.848}}},
            "exp_b": {"time_rows": {20: {
                "protons_delivered_raw": 2057.44,
                "protons_delivered_positive": 2057.44}}},
        }
        rows = ua.build_radex_time_rows(products, ["exp_a", "exp_b"])
        self.assertAlmostEqual(rows[10]["radex_mean_fluence_raw_protons_cm2"], 10.0)
        self.assertAlmostEqual(rows[20]["radex_mean_fluence_raw_protons_cm2"], 20.0)
        self.assertAlmostEqual(rows[10]["radex_instrument_area_cm2"], 102.872)
        self.assertAlmostEqual(rows[20]["radex_instrument_area_cm2"], 102.872)



if __name__ == "__main__":
    unittest.main()
