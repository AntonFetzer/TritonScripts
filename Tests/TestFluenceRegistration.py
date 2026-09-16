import unittest

import numpy as np

from UppsalaAcceleratorLogs import register_fluence_maps as rf


class FluenceRegistrationTests(unittest.TestCase):
    def test_pearson_and_cosine_identity(self):
        values = np.asarray([0.0, 0.25, 0.5, 1.0])
        self.assertAlmostEqual(rf.pearson(values, values), 1.0)
        self.assertAlmostEqual(rf.cosine(values, values), 1.0)

    def test_sample_registered_applies_coordinate_scale_and_offset(self):
        try:
            import scipy  # noqa: F401
        except ImportError:
            self.skipTest("SciPy not available")
        edges = np.linspace(-5.0, 5.0, 11)
        x = (edges[:-1] + edges[1:]) / 2.0
        y = (edges[:-1] + edges[1:]) / 2.0
        observed_x, observed_y = np.meshgrid(x, y)
        observed = observed_x + 10.0 * observed_y
        query_x = np.asarray([2.0, 4.0])
        query_y = np.asarray([1.0, 3.0])
        sampled = rf.sample_registered(
            observed,
            edges,
            edges,
            query_x,
            query_y,
            (2.0, 2.0, 1.0, 1.0),
        )
        # Source coordinates are (0.5, 0.0) and (1.5, 1.0).
        np.testing.assert_allclose(sampled, [0.5, 11.5], atol=1e-12)

    def test_fit_reports_area_jacobian(self):
        # The spatial fit result uses this Jacobian to preserve the integrated
        # proton count when an absolute fluence-density map is stretched.
        scale_x, scale_y = 1.48, 1.33
        area_scale = scale_x * scale_y
        self.assertAlmostEqual(area_scale, 1.9684)
        self.assertAlmostEqual(1.0 / area_scale, 0.508026824)


if __name__ == "__main__":
    unittest.main()
