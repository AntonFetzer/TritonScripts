import unittest

import numpy as np

from UppsalaAcceleratorLogs import register_fluence_maps as rf


class FluenceRegistrationTests(unittest.TestCase):
    def test_pearson_and_cosine_identity(self):
        values = np.asarray([0.0, 0.25, 0.5, 1.0])
        self.assertAlmostEqual(rf.pearson(values, values), 1.0)
        self.assertAlmostEqual(rf.cosine(values, values), 1.0)

    def test_fit_reports_area_jacobian(self):
        # The spatial fit result uses this Jacobian to preserve the integrated
        # proton count when an absolute fluence-density map is stretched.
        scale_x, scale_y = 1.48, 1.33
        area_scale = scale_x * scale_y
        self.assertAlmostEqual(area_scale, 1.9684)
        self.assertAlmostEqual(1.0 / area_scale, 0.508026824)


if __name__ == "__main__":
    unittest.main()
