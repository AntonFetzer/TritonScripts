"""Regression test for Dependencies.SEERates and the result-file handling of the LET and fluence readers.

These steps used to exist as separate copies in the three SEU rate scripts and
the three rate plots. What is pinned here is what those copies did differently
or could get silently wrong: which runs are divided by the mission duration, how
a zero-mean bin is folded, the order of the thicknesses, and that a missing or
incomplete result file is skipped, without deleting it or stopping the pass.

Run from the repository root:

    python3 -m unittest discover -s Tests -p "Test*.py"
"""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import matplotlib
matplotlib.use("Agg")  # never open a window from the test suite
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.HistogramPlots import histogramWithCrossSection
from Dependencies.MergeHistograms import completeResultFiles
from Dependencies.SEERates import (
    LET_CROSS_SECTIONS, MISSION_SECONDS, PROTON_CROSS_SECTIONS, RATE_CSV_HEADER, crossSection, foldCrossSection,
    normaliseToRate, readRateCSVs, runFolders, seeRate, thicknessFolders, usableHistogram, writeRateCSV)
from Dependencies.TotalFluenceHistograms import totalFluenceHistos
from Dependencies.TotalLETHistos import totalLETHistos


def histogram():
    """Two bins, the first one with a zero mean."""
    return {
        'lower': np.array([1.0, 2.0]), 'upper': np.array([2.0, 3.0]), 'mean': np.array([0.0, 2.5]),
        'value': np.array([10.0, 20.0]), 'error': np.array([1.0, 2.0]), 'entries': np.array([3, 4]),
    }


class TestCrossSection(unittest.TestCase):

    def test_weibull_is_zero_at_and_below_the_onset(self):
        Parameters = {"L0": 1.0, "W": 2.0, "S": 1.5, "A0": 1e-8}
        Result = crossSection(np.array([0.5, 1.0, 3.0]), Parameters)
        self.assertEqual(Result[0], 0.0)
        self.assertEqual(Result[1], 0.0)
        self.assertAlmostEqual(Result[2] / 1e-8, 1 - np.exp(-1.0))

    def test_measured_points_interpolate_log_log(self):
        Parameters = {"E_pts": np.array([1.0, 100.0]), "S_pts": np.array([1e-10, 1e-14])}
        Result = crossSection(np.array([0.5, 1.0, 10.0, 100.0, 1000.0]), Parameters)
        self.assertEqual(Result[0], 0.0)  # zero below the lowest point
        self.assertAlmostEqual(Result[1] / 1e-10, 1.0)
        self.assertAlmostEqual(Result[2] / 1e-12, 1.0)
        self.assertAlmostEqual(Result[3] / 1e-14, 1.0)
        self.assertAlmostEqual(Result[4] / 1e-14, 1.0)  # constant above the highest point

    def test_device_parameters_are_in_cm2_per_bit(self):
        LET = {X["name"]: X for X in LET_CROSS_SECTIONS}
        Proton = {X["name"]: X for X in PROTON_CROSS_SECTIONS}
        self.assertEqual(LET["NanoXplore NG-Medium"]["A0"], 5.2e-9)
        self.assertEqual(LET["Cypress CY62167GE30-45ZXI"]["A0"], 2.6e-7)
        self.assertEqual(Proton["NanoXplore NG-Medium"]["A0"], 4.85e-16)
        self.assertAlmostEqual(np.max(Proton["Cypress CY62167GE30-45ZXI"]["S_pts"]) / 1.2757607330635113e-9, 1.0)


class TestNormaliseToRate(unittest.TestCase):

    def test_only_mission_fluence_runs_are_divided(self):
        for Folder in ("GEO-solar-proton", "LEO-electron", "GEO-cosmic-iron", "GEO-SAPPHIRE-Solar-Iron"):
            self.assertAlmostEqual(normaliseToRate(histogram(), Folder)['value'][1] * MISSION_SECONDS, 20.0, msg=Folder)
        for Folder in ("GEO-SolarProton-5minPeakFlux", "GEO-CREME96-Solar-Iron", "LEO-SAPPHIRE-Solar-Iron", "VAP-GCR-Iron", "Carrington-SEP-Expected-EVT2026"):
            self.assertEqual(normaliseToRate(histogram(), Folder)['value'][1], 20.0, msg=Folder)

    def test_input_is_not_changed(self):
        Hist = histogram()
        normaliseToRate(Hist, "GEO-solar-proton")
        self.assertEqual(Hist['value'][1], 20.0)


class TestFoldAndRate(unittest.TestCase):

    def test_zero_mean_uses_the_lower_edge_and_input_is_not_changed(self):
        Hist = histogram()
        SEEHist = foldCrossSection(Hist, {"L0": 0.5, "W": 1.0, "S": 1.0, "A0": 1.0})
        self.assertEqual(SEEHist['mean'][0], 1.0)
        self.assertGreater(SEEHist['value'][0], 0.0)
        self.assertEqual(Hist['mean'][0], 0.0)
        self.assertEqual(Hist['value'][0], 10.0)

    def test_rate_is_not_energy_weighted(self):
        """The rate is value x cross section: the same flux gives the same rate wherever the bins lie."""
        Hist = histogram()
        Hist['mean'] = np.array([1.5, 2.5])
        Flat = {"E_pts": np.array([0.1, 1000.0]), "S_pts": np.array([1e-10, 1e-10])}
        self.assertAlmostEqual(seeRate(foldCrossSection(Hist, Flat))[0] / 1e-10, np.sum(Hist['value']))

    def test_rate_error_and_contributing_entries(self):
        SEEHist = {'value': np.array([0.0, 1.0, 3.0]), 'error': np.array([0.0, 3.0, 4.0]), 'entries': np.array([10, 20, 30])}
        Rate, Error, Relative, Entries = seeRate(SEEHist)
        self.assertEqual((Rate, Error, Relative, Entries), (4.0, 5.0, 1.25, 50))

    def test_zero_rate_has_no_relative_error(self):
        SEEHist = {'value': np.zeros(2), 'error': np.zeros(2), 'entries': np.array([1, 2])}
        Rate, _, Relative, Entries = seeRate(SEEHist)
        self.assertEqual(Rate, 0.0)
        self.assertTrue(np.isnan(Relative))
        self.assertEqual(Entries, 0)

    def test_inconsistent_or_missing_histograms_are_not_usable(self):
        Bad = histogram()
        Bad['error'] = np.array([1.0])
        with mock.patch("builtins.print"):
            self.assertTrue(usableHistogram(histogram(), "Res/"))
            self.assertFalse(usableHistogram(None, "Res/"))
            self.assertFalse(usableHistogram(Bad, "Res/"))


class TestFolders(unittest.TestCase):

    def test_run_folders_default_skips_what_is_not_a_run(self):
        with tempfile.TemporaryDirectory() as Root:
            for Name in ("RunB", "RunA", "OldResults", "NoSEEs"):
                (Path(Root) / Name).mkdir()
            (Path(Root) / "notes.md").write_text("x")
            self.assertEqual(runFolders(Root + "/", []), (Root + "/", ["RunA", "RunB"]))

    def test_run_folders_from_arguments(self):
        with tempfile.TemporaryDirectory() as Root:
            Directory, FolderList = runFolders("/unused/", [Root + "/RunA/", Root + "/RunB"])
            self.assertEqual(Directory, Root + "/")
            self.assertEqual(FolderList, ["RunA", "RunB"])

    def test_run_folders_must_share_a_directory(self):
        with self.assertRaises(SystemExit):
            runFolders("/unused/", ["/a/RunA", "/b/RunB"])

    def test_thicknesses_are_in_natural_order(self):
        with tempfile.TemporaryDirectory() as Root:
            for Name in ("16mm", "1mm", "2mm", "8mm", "4mm"):
                (Path(Root) / "Run" / Name).mkdir(parents=True)
            (Path(Root) / "Run" / "Run SEE-Rate.csv").write_text("x")
            self.assertEqual(thicknessFolders(Root, "Run"), ["1mm", "2mm", "4mm", "8mm", "16mm"])


class TestRateCSVs(unittest.TestCase):

    Rows = [("1mm", 1.5, 0.5, 1 / 3, 10), ("16mm", 0.25, 0.125, 0.5, 2)]

    def test_round_trip_keeps_values_and_order(self):
        with tempfile.TemporaryDirectory() as Root:
            writeRateCSV(Root + "/RunA Dev SEE-Rate.csv", "RunA", "Dev", self.Rows)
            self.assertEqual(Path(Root + "/RunA Dev SEE-Rate.csv").read_text().splitlines()[0], RATE_CSV_HEADER)
            Dict = readRateCSVs([Root + "/* Dev SEE-Rate.csv"], "Dev")
        self.assertEqual(list(Dict["RunA"]['Shielding']), ["1mm", "16mm"])
        self.assertEqual(list(Dict["RunA"]['SEE_Rate']), [1.5, 0.25])
        self.assertEqual(list(Dict["RunA"]['SEE_Error']), [0.5, 0.125])
        self.assertEqual(list(Dict["RunA"]['Entries_Contributing_To_SEE']), [10.0, 2.0])
        self.assertIsInstance(Dict["RunA"]['SEE_Rate'], np.ndarray)

    def test_several_patterns_are_collected(self):
        with tempfile.TemporaryDirectory() as Root:
            writeRateCSV(Root + "/RunA Dev SEE-Rate.csv", "RunA", "Dev", self.Rows)
            writeRateCSV(Root + "/Total Dev SEE-Rate.csv", "Total", "Dev", self.Rows)
            Dict = readRateCSVs([Root + "/RunA * SEE-Rate.csv", Root + "/Total * SEE-Rate.csv"], "Dev")
        self.assertEqual(sorted(Dict), ["RunA", "Total"])

    def test_wrong_cross_section_is_refused(self):
        with tempfile.TemporaryDirectory() as Root:
            writeRateCSV(Root + "/RunA Dev SEE-Rate.csv", "RunA", "Other", self.Rows)
            with self.assertRaises(ValueError):
                readRateCSVs([Root + "/* SEE-Rate.csv"], "Dev")

    def test_same_dataset_from_two_files_is_refused(self):
        with tempfile.TemporaryDirectory() as Root:
            writeRateCSV(Root + "/One SEE-Rate.csv", "RunA", "Dev", self.Rows)
            writeRateCSV(Root + "/Two SEE-Rate.csv", "RunA", "Dev", self.Rows)
            with self.assertRaises(ValueError):
                readRateCSVs([Root + "/* SEE-Rate.csv"], "Dev")

    def test_no_rows_writes_no_file(self):
        with tempfile.TemporaryDirectory() as Root, mock.patch("builtins.print"):
            writeRateCSV(Root + "/RunA Dev SEE-Rate.csv", "RunA", "Dev", [])
            self.assertFalse(Path(Root + "/RunA Dev SEE-Rate.csv").exists())


class TestResultFiles(unittest.TestCase):

    def test_incomplete_files_are_skipped_and_kept(self):
        with tempfile.TemporaryDirectory() as Root, mock.patch("builtins.print"):
            (Path(Root) / "LET_1.csv").write_text("data\n'End of File'\n")
            (Path(Root) / "LET_2.csv").write_text("truncated")
            (Path(Root) / "Fluence_1.csv").write_text("data\n'End of File'\n")
            self.assertEqual(completeResultFiles(Root, "LET"), ["LET_1.csv"])
            self.assertTrue((Path(Root) / "LET_2.csv").exists())

    def test_missing_folder_has_no_files(self):
        self.assertEqual(completeResultFiles("/no/such/folder", "LET"), [])

    def test_readers_return_nothing_without_complete_files(self):
        with tempfile.TemporaryDirectory() as Root, mock.patch("builtins.print"):
            (Path(Root) / "LET_1.csv").write_text("truncated")
            (Path(Root) / "Fluence_1.csv").write_text("truncated")
            self.assertEqual(totalLETHistos(Root), (None, None))
            self.assertEqual(totalFluenceHistos(Root), (None, None))
            self.assertTrue((Path(Root) / "LET_1.csv").exists())
            self.assertTrue((Path(Root) / "Fluence_1.csv").exists())


class TestHistogramWithCrossSection(unittest.TestCase):

    def test_writes_the_pdf_and_closes_the_figure(self):
        Hist = histogram()
        Hist['mean'] = np.array([1.5, 2.5])  # positive for the log axes
        Before = plt.get_fignums()
        with tempfile.TemporaryDirectory() as Root:
            Path_ = Path(Root) / "hist.pdf"
            Written = histogramWithCrossSection(Hist, np.full(2, 1e-9), Path_, "title", "x", "y", "hist",
                                                "sigma", "sigma axis", legendLoc="upper left", capsize=2)
            self.assertTrue(Path_.exists())
        self.assertEqual(Written, str(Path_))
        self.assertEqual(plt.get_fignums(), Before)


if __name__ == "__main__":
    unittest.main()
