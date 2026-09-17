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


if __name__ == "__main__":
    unittest.main()
