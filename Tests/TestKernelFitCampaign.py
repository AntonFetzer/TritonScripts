"""Regression test for the campaign file of RadEx.KernelFit with a table that carries its own errors.

The HUS electron campaign keeps both phases in one table, selected by `rows`, gives
each channel's dose errors as columns, reads the phase fluence from a column and
names its device-type column. What is pinned here is that those rows and columns
are the ones used, that sigma^2 = noise^2 + conversion^2 + (systematic x dose)^2
in the table's unit, that --systematic overrides the file, and that the scale
groups leave excluded channels out.

Run from the repository root:

    python3 -m unittest discover -s Tests -p "Test*.py"
"""
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from RadEx import KernelFit as kf

TABLE = """phase,channel,radfet_type,dose_fit_krad,noise_krad,conversion_uncertainty_krad,fluence
late,0,VT01,10.0,0.3,0.4,2e12
late,1,VT01,20.0,0.0,,2e12
late,2,VT05,5.0,0.1,0.0,2e12
first,0,VT01,1.0,0.1,0.1,1e12
first,1,VT01,2.0,0.1,0.1,1e12
first,2,VT05,3.0,0.1,0.1,1e12
"""


def phase(rows):
    return {"rows": {"phase": rows}, "dose_column": "dose_fit_krad", "noise_column": "noise_krad",
            "conversion_column": "conversion_uncertainty_krad", "log_fluence_column": "fluence"}


class KernelFitCampaignTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        folder = Path(self.tmp.name)
        (folder / "doses.csv").write_text(TABLE)
        self.path = folder / "campaign.json"
        self.path.write_text(json.dumps({
            "dose": "doses.csv", "device_column": "radfet_type", "fit_space": "dose",
            "systematic_fraction": 0.1, "phases": {"late": phase("late"), "first": phase("first")},
            "layers": {"oxide": 0, "die": 3}, "labels": {"0": "a", "1": "b", "2": "c"}}))

    def tearDown(self):
        self.tmp.cleanup()

    def test_rows_columns_and_sigma(self):
        campaign = kf.load_campaign(self.path)
        self.assertEqual(campaign["channels"], [0, 1, 2])
        self.assertEqual(campaign["device"], {0: "VT01", 1: "VT01", 2: "VT05"})
        self.assertEqual(campaign["phases"]["late"]["log_fluence"], 2e12)
        self.assertEqual(campaign["phases"]["first"]["log_fluence"], 1e12)
        meas = kf.measurement(campaign, "late")
        np.testing.assert_allclose(meas["value"], [10.0, 20.0, 5.0])
        expected = np.sqrt(np.array([0.3 ** 2 + 0.4 ** 2, 0.0, 0.1 ** 2]) + (0.1 * np.array([10.0, 20.0, 5.0])) ** 2)
        np.testing.assert_allclose(meas["sigma"], expected)
        np.testing.assert_allclose(kf.measurement(campaign, "first")["value"], [1.0, 2.0, 3.0])

    def test_systematic_override_and_passthrough(self):
        campaign = kf.load_campaign(self.path, systematic=0.2)
        np.testing.assert_allclose(kf.measurement(campaign, "late")["sigma"][1], 4.0)
        self.assertEqual(campaign["labels"], {0: "a", 1: "b", 2: "c"})
        self.assertEqual(kf.layer_slice(campaign, "die"), slice(3, 6))
        np.testing.assert_allclose(kf.channel_weights(campaign, "late", campaign["channels"]), 1.0)

    def test_scale_groups_exclude(self):
        campaign = kf.load_campaign(self.path)
        groups = kf.scale_groups(campaign, campaign["channels"], "per-type", exclude=(1,))
        self.assertEqual(list(groups), ["VT01", "VT05"])
        np.testing.assert_array_equal(groups["VT01"], [0])
        np.testing.assert_array_equal(kf.scale_groups(campaign, campaign["channels"], "shared", (2,))["all"], [0, 1])


if __name__ == "__main__":
    unittest.main()
