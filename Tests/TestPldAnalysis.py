import csv
import tempfile
import unittest
from pathlib import Path

from UppsalaAcceleratorLogs import pld_analysis as pa


class PldParserTests(unittest.TestCase):
    def write_plan(self, directory: str) -> Path:
        path = Path(directory) / "test.pld"
        rows = [
            ["Beam", "id", "name", "initial", "first", "plan", "beam", "100", "2", "1"],
            ["Layer", "Spot1", "85", "2", "4", "10"],
            ["Element", "-1", "-2", "0", "0"],
            ["Element", "-1", "-2", "1", "0"],
            ["Element", "1", "2", "0", "0"],
            ["Element", "1", "2", "1", "0"],
        ]
        with path.open("w", encoding="utf-8", newline="") as handle:
            csv.writer(handle).writerows(rows)
        return path

    def test_parse_collapses_element_pairs_and_swaps_machine_axes(self):
        with tempfile.TemporaryDirectory() as tmp:
            plan = pa.parse_pld(self.write_plan(tmp))
        self.assertEqual(plan.total_mu, 100.0)
        self.assertEqual(plan.parsed_weight_sum, 2.0)
        self.assertEqual(len(plan.spots), 2)
        self.assertEqual(plan.spots[0].machine_x_mm, -2.0)
        self.assertEqual(plan.spots[0].machine_y_mm, -1.0)
        self.assertEqual(plan.spots[0].paintings, 10)

    def test_reconstruction_conserves_planned_protons(self):
        conversion = [
            {"energy_mev": 80.0, "protons_per_mu": 6.65e7},
            {"energy_mev": 85.0, "protons_per_mu": 6.97e7},
        ]
        grid = {
            "minimum_mm": -30.0,
            "maximum_mm": 30.0,
            "pixel_mm": 1.0,
            "width_interpretation": "sigma",
        }
        with tempfile.TemporaryDirectory() as tmp:
            plan = pa.parse_pld(self.write_plan(tmp))
            summary, arrays, spots = pa.reconstruct_plan(plan, conversion, grid)
        self.assertEqual(len(spots), 2)
        self.assertAlmostEqual(summary["planned_protons"], 100.0 * 6.97e7)
        self.assertLess(abs(summary["heatmap_conservation_relative_error"]), 1e-12)
        self.assertAlmostEqual(
            float(arrays["planned_protons_per_pixel"].sum()),
            summary["planned_protons"],
            delta=summary["planned_protons"] * 1e-6,
        )

    def test_parser_rejects_mismatched_pair_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write_plan(tmp)
            text = path.read_text(encoding="utf-8").replace(
                "Element,-1,-2,1,0", "Element,-1,-3,1,0"
            )
            path.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "spot coordinates differ"):
                pa.parse_pld(path)


if __name__ == "__main__":
    unittest.main()
