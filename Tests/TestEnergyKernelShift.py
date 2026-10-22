import unittest

import numpy as np

from Dependencies.EnergyKernel import fitShiftAmplitudes, shiftResponse


class TestShiftResponse(unittest.TestCase):
    def test_pure_power_law_and_knee(self):
        self.assertAlmostEqual(shiftResponse(100.0, 0.5), 10.0)
        self.assertEqual(shiftResponse(-5.0, 0.5), 0.0)
        # linear below the knee: slope n D0^(n-1)
        d0, n = 360.0, 0.51
        small = 1e-3
        self.assertAlmostEqual(shiftResponse(small, n, d0) / small, n * d0 ** (n - 1), places=6)


class TestFitShiftAmplitudes(unittest.TestCase):
    def setUp(self):
        self.dose = np.array([16000.0, 20000.0, 14000.0, 41000.0, 0.0])
        self.exponent = np.array([0.51, 0.51, 0.51, 0.64, 0.64])
        self.offset = np.array([360.0, 360.0, 360.0, 0.0, 0.0])
        self.groups = {"VT01": np.array([0, 1, 2]), "VT05": np.array([3, 4])}

    def test_recovers_amplitudes(self):
        truth = {"VT01": 0.012, "VT05": 0.00025}
        shift = np.array([truth["VT01" if i < 3 else "VT05"] * shiftResponse(d, n, o)
                          for i, (d, n, o) in enumerate(zip(self.dose, self.exponent, self.offset))])
        amplitudes, chi2, pulls, predicted = fitShiftAmplitudes(
            shift, np.full(5, 0.005), self.exponent, self.offset, 0.15,
            self.dose, 0.02 * self.dose, self.groups)
        self.assertAlmostEqual(amplitudes["VT01"], truth["VT01"])
        self.assertAlmostEqual(amplitudes["VT05"], truth["VT05"])
        self.assertAlmostEqual(chi2, 0.0)
        np.testing.assert_allclose(predicted, shift)
        self.assertEqual(predicted[4], 0.0)

    def test_previous_dose_and_pull_sign(self):
        shift = np.array([1.5, 1.6, 1.4, 0.26, 0.0])
        a0, chi0, _, _ = fitShiftAmplitudes(shift, np.full(5, 0.005), self.exponent, self.offset, 0.15,
                                           self.dose, 0.0 * self.dose, self.groups)
        a1, _, pulls, predicted = fitShiftAmplitudes(shift, np.full(5, 0.005), self.exponent, self.offset,
                                                     0.15, self.dose, 0.0 * self.dose, self.groups,
                                                     previous=np.full(5, 5000.0))
        self.assertLess(a1["VT01"], a0["VT01"])
        self.assertTrue(np.all(np.sign(pulls[:3]) == np.sign(predicted[:3] - shift[:3])))


if __name__ == "__main__":
    unittest.main()
