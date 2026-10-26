import unittest

import numpy as np

from Dependencies.EnergyKernel import fitShiftAmplitudes, foldKernel, shiftResponse, tabulatedResponse


class TestShiftResponse(unittest.TestCase):
    def test_pure_power_law_and_knee(self):
        self.assertAlmostEqual(shiftResponse(100.0, 0.5), 10.0)
        self.assertEqual(shiftResponse(-5.0, 0.5), 0.0)
        # linear below the knee: slope n D0^(n-1)
        d0, n = 360.0, 0.51
        small = 1e-3
        self.assertAlmostEqual(shiftResponse(small, n, d0) / small, n * d0 ** (n - 1), places=6)


class TestTabulatedResponse(unittest.TestCase):
    def test_interpolates_and_extends_in_log_log(self):
        table = (np.array([1000.0, 10000.0, 100000.0]), np.array([0.4, 2.0, 5.6]))
        self.assertAlmostEqual(tabulatedResponse(10000.0, table), 2.0)
        self.assertAlmostEqual(tabulatedResponse(np.sqrt(1e3 * 1e4), table), np.sqrt(0.4 * 2.0))
        self.assertAlmostEqual(tabulatedResponse(100.0, table), 0.4 * (0.1) ** np.log10(5.0))
        self.assertEqual(tabulatedResponse(0.0, table), 0.0)

    def test_amplitude_is_ratio_to_table(self):
        grid = np.geomspace(10.0, 1e6, 200)
        table = (grid, 0.3 * shiftResponse(grid, 0.45, 200.0))
        dose = np.array([16000.0, 20000.0, 14000.0])
        shift = 0.52 * tabulatedResponse(dose, table)
        amplitudes, chi2, _, predicted = fitShiftAmplitudes(
            shift, np.full(3, 0.005), np.full(3, 0.51), np.full(3, 360.0), 0.15,
            dose, 0.02 * dose, {"VT01": np.arange(3)}, response=[table] * 3)
        self.assertAlmostEqual(amplitudes["VT01"], 0.52, places=6)
        self.assertAlmostEqual(chi2, 0.0, places=6)
        np.testing.assert_allclose(predicted, shift)


class TestFoldKernel(unittest.TestCase):
    def test_gaussian_and_upper_cut(self):
        # A kernel linear in energy returns the spectrum's mean energy.
        energies = np.arange(55.0, 70.0001, 0.25)
        kernel = np.column_stack([energies, np.ones_like(energies)])
        mu, sigma = 62.5, 0.3
        full, _ = foldKernel(energies, kernel, 0.01 * kernel, [mu], sigma)
        self.assertAlmostEqual(full[0, 0], mu, places=6)
        self.assertAlmostEqual(full[0, 1], 1.0)
        same, _ = foldKernel(energies, kernel, 0.01 * kernel, [mu], sigma, upperCut=None)
        np.testing.assert_array_equal(same, full)
        # Cut off above mu + 1 sigma: mean of the one-sided truncated normal,
        # mu - sigma phi(1) / Phi(1) = mu - 0.2876 sigma, to the 0.025 MeV grid's
        # discretisation (the cut point enters with full weight).
        cut, _ = foldKernel(energies, kernel, 0.01 * kernel, [mu], sigma, upperCut=1.0)
        self.assertAlmostEqual(cut[0, 0], mu - 0.2876 * sigma, delta=0.02 * sigma)
        self.assertAlmostEqual(cut[0, 1], 1.0)


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
