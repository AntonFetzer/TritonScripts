import unittest

import numpy as np

from UppsalaAcceleratorLogs import radex_isocentre_outputs as rio


class RadexIsocentreOutputTests(unittest.TestCase):
    def test_bilinear_sample_reproduces_linear_field(self):
        edges = np.linspace(-2.0, 2.0, 5)
        centers = (edges[:-1] + edges[1:]) / 2.0
        x, y = np.meshgrid(centers, centers)
        values = 2.0 * x + 3.0 * y
        self.assertAlmostEqual(
            rio.bilinear_sample(values, edges, edges, 0.25, -0.25),
            -0.25,
        )

    def test_radfet_locations_match_documented_geometry(self):
        self.assertEqual(len(rio.RADEX_RADFETS), 12)
        self.assertEqual(rio.RADEX_RADFETS[0], (0, "VT01, 100% Pb", 60.0, 16.7))
        self.assertEqual(
            rio.RADEX_RADFETS[10],
            (10, "VT01, readout PCB back", 28.0, -20.0),
        )
        self.assertEqual(rio.RADEX_RADFETS[11], (11, "VT05, readout PCB back", 0.0, 0.0))

    def test_transform_constants_are_internally_consistent(self):
        transform = rio.RADEX_CHAMBER_TO_ISOCENTRE
        area = (
            transform["scale_x_chamber_to_isocentre"]
            * transform["scale_y_chamber_to_isocentre"]
        )
        self.assertAlmostEqual(area, transform["area_scale_chamber_to_isocentre"])
        self.assertAlmostEqual(
            1.0 / area,
            transform["fluence_density_jacobian_chamber_to_isocentre"],
        )


    def test_facility_spot_sigma_interpolates_linearly_without_clamping(self):
        self.assertEqual(rio.facility_spot_sigma(85.0), (5.4, 5.4))
        sigma_x, sigma_y = rio.facility_spot_sigma(64.0)
        self.assertAlmostEqual(sigma_x, 6.68)
        self.assertAlmostEqual(sigma_y, 6.68)
        with self.assertRaises(ValueError):
            rio.facility_spot_sigma(90.0)

    def test_forward_transform_conserves_protons_and_moves_points(self):
        edges = np.linspace(-50.0, 50.0, 101)
        values = np.zeros((100, 100))
        values[50, 60] = 3.0  # chamber pixel centred on x=10.5, y=0.5 mm
        values[40, 45] = 1.0
        transform = dict(rio.RADEX_CHAMBER_TO_ISOCENTRE)
        mapped, outside = rio.transform_positions_to_isocentre(
            values, edges, edges, edges, edges, transform
        )
        self.assertEqual(outside, 0.0)
        self.assertAlmostEqual(mapped.sum(), values.sum(), places=12)
        centres = (edges[:-1] + edges[1:]) / 2.0
        x, y = np.meshgrid(centres, centres)
        weight = mapped.copy()
        weight[y < -5.0] = 0.0
        x_mean = float((weight * x).sum() / weight.sum())
        expected = transform["scale_x_chamber_to_isocentre"] * 10.5 + transform[
            "offset_x_isocentre_mm"
        ]
        self.assertAlmostEqual(x_mean, expected, places=9)

    def test_spot_lines_recover_a_comb_lattice(self):
        edges = np.arange(-30.0, 30.0001, 0.05)
        centres = (edges[:-1] + edges[1:]) / 2.0
        pitch = 1.35
        lines = -20.25 + pitch * np.arange(31)
        projection = np.zeros_like(centres)
        for index, line in enumerate(lines):
            weight = 0.5 if index in (0, 30) else 1.0
            projection += weight * np.exp(-0.5 * ((centres - line) / 0.1) ** 2)
        positions, weights = rio.spot_lines(projection, edges)
        self.assertEqual(len(positions), 31)
        self.assertTrue(np.allclose(positions, lines, atol=0.01))
        self.assertAlmostEqual(float(weights[0] / np.median(weights)), 0.5, places=2)
        statistics = rio.lattice_statistics(positions, weights, rio.PLD_SPOT_PITCH_MM)
        self.assertEqual(statistics["lines_used"], 31)
        self.assertAlmostEqual(statistics["pitch_median_mm"], pitch, places=3)
        self.assertAlmostEqual(statistics["span_mm"], 30 * pitch, places=2)
        self.assertAlmostEqual(
            statistics["implied_scale_chamber_to_isocentre"], 2.0 / pitch, places=3
        )

    def test_row_profile_samples_a_linear_field(self):
        edges = np.linspace(-10.0, 10.0, 21)
        centres = (edges[:-1] + edges[1:]) / 2.0
        x, y = np.meshgrid(centres, centres)
        values = 1.0 + 0.1 * x
        profile = rio.row_profile(values, edges, edges, 0.0, np.array([-5.0, 0.0, 5.0]))
        self.assertTrue(np.allclose(profile, [0.5, 1.0, 1.5]))


if __name__ == "__main__":
    unittest.main()
