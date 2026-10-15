"""Regression test for totalDose(..., replica_error=True).

GRAS 5.02.01 quotes wrong errors for importance-weighted sources
(GRASStatDouble::fill treats the event weight as a frequency weight), so the
error of such runs has to come from the scatter between replica files. The
ReplicaRun fixture has five replicas with unequal entry counts whose doses
scatter far less than their quoted errors, the signature of that bug.

This pins:

* The replica error is sqrt(sum(N_i * (y_i - mean)^2) / ((k - 1) * sum(N_i))),
  the error of the entry-weighted mean, with weight N_i rather than N_i^2.
* Quoted errors far above the scatter trigger a warning, in either mode, but
  the chance low scatter of two consistent files does not.
* Fewer files than min_replicas is refused.
* The default mode still returns the pooled quoted error and the same keys.

Run from the repository root:

    python3 -m unittest discover -s Tests -p "Test*.py"
"""
import contextlib
import io
import sys
import unittest
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.TotalDose import totalDose

FIXTURES = Path(__file__).resolve().parent / "Fixtures"
REPLICA_DIR = FIXTURES / "ReplicaRun"
TWO_FILE_DIR = FIXTURES / "TwoFileRun"

# Fixture contents, in rad as GRAS writes them; rows are files, columns tiles.
ENTRIES = np.array([1e6, 2e6, 1e6, 3e6, 1e6])
DOSE = np.array([
    [1.000e-6, 2.00e-6, 4.00e-6],
    [1.010e-6, 1.98e-6, 4.02e-6],
    [0.990e-6, 2.03e-6, 3.97e-6],
    [1.005e-6, 2.01e-6, 4.01e-6],
    [0.995e-6, 1.99e-6, 3.99e-6],
])
QUOTED_ERROR = 2.0e-7
RAD_TO_KRAD = 1.0e-3


def run(**kwargs):
    """Call totalDose and return its result together with what it printed."""
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        result = totalDose(str(REPLICA_DIR), **kwargs)
    return result, out.getvalue()


class TestReplicaError(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not REPLICA_DIR.is_dir():
            raise unittest.SkipTest("fixture directory missing: {}".format(REPLICA_DIR))
        cls.replica, cls.replica_log = run(replica_error=True)
        cls.quoted, cls.quoted_log = run()

    def test_replica_error_of_entry_weighted_mean(self):
        n = ENTRIES[:, None]
        mean = np.sum(n * DOSE, axis=0) / ENTRIES.sum()
        expected = np.sqrt(
            np.sum(n * (DOSE - mean) ** 2, axis=0) / ((len(ENTRIES) - 1) * ENTRIES.sum())
        ) * RAD_TO_KRAD
        np.testing.assert_allclose(self.replica["error"], expected, rtol=1e-12)
        np.testing.assert_allclose(self.replica["dose"], mean * RAD_TO_KRAD, rtol=1e-12)

    def test_default_mode_unchanged(self):
        expected = QUOTED_ERROR * np.sqrt(np.sum(ENTRIES ** 2)) / ENTRIES.sum() * RAD_TO_KRAD
        np.testing.assert_allclose(self.quoted["error"], np.full(3, expected), rtol=1e-12)
        np.testing.assert_array_equal(self.quoted["dose"], self.replica["dose"])
        np.testing.assert_array_equal(self.quoted["entries"], self.replica["entries"])
        np.testing.assert_array_equal(self.quoted["non-zeros"], self.replica["non-zeros"])
        self.assertEqual(set(self.quoted), {"dose", "error", "entries", "non-zeros"})
        self.assertEqual(set(self.replica), set(self.quoted))

    def test_overestimated_quoted_errors_are_flagged(self):
        self.assertIn("quoted statistical errors are overestimated", self.quoted_log)
        self.assertIn("quoted statistical errors are overestimated", self.replica_log)

    def test_consistent_files_are_not_flagged(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            totalDose(str(TWO_FILE_DIR))
        self.assertNotIn("overestimated", out.getvalue())

    def test_too_few_replicas_refused(self):
        with contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaises(ValueError):
                totalDose(str(TWO_FILE_DIR), replica_error=True)
            with self.assertRaises(ValueError):
                totalDose(str(REPLICA_DIR), replica_error=True, min_replicas=6)


if __name__ == "__main__":
    unittest.main()
