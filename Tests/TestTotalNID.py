"""Regression test for Read.ReadNID and Dependencies.TotalNID.

What is pinned here is what a fresh displacement-damage aggregation gets wrong:
averaging by file instead of by entries, pooling errors as a plain quadrature
sum, and accepting a result directory that is short a few files because the
job is still running.

Run from the repository root:

    python3 -m unittest discover -s Tests -p "Test*.py"
"""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.TotalNID import totalNID, moduleRecords  # noqa: E402
from Read.ReadNID import readNID  # noqa: E402

BLOCK = """'*', 9, 0, 4, 0, 4, 4, 1, 2
'GRAS_DATA_TITLE',   -1,'TOTAL NID'
'GRAS_DATA_TYPE',   -1,'STAT_DOUBLE'
'GRAS_MODULE_NAME',   -1,'{name}'
'GRAS_MODULE_TYPE',   -1,'NID'
'NID','MeV/g',    1,'NID'
'Error','MeV/g',    1,'Error NID'
'Entries','',    1,'Number of entries'
'Non zero entries','',    1,'Number of non zero entries'
 {nid},   {error},     {entries},          {nz}
'End of Block'
"""

DOSE_BLOCK = """'*', 10, 0, 5, 0, 4, 4, 2, 4
'GRAS_DATA_TITLE',   -1,'TOTAL DOSE FOR INDIVIDUAL VOLUMES'
'GRAS_DATA_TYPE',   -1,'STAT_DOUBLE'
'GRAS_MODULE_NAME',   -1,'doseDetector'
'GRAS_MODULE_TYPE',   -1,'DOSE'
'VOLUME NAMES',   -2,'VT01_gox_0_PV','VT01_die_0_PV'
'Dose','rad',    1,'Dose/energy deposition'
'Error','rad',    1,'Error dose/energy deposition'
'Entries','',    1,'Number of entries'
'Non zero entries','',    1,'Number of non zero entries'
  1.455e-08,  1.4315e-09,       3e+05,         419
 1.8483e-08,  9.3127e-10,       3e+05,         579
'End of Block'
"""


def writeFile(directory, name, blocks):
    path = Path(directory) / name
    path.write_text("".join(blocks) + "'End of File'\n", encoding="utf-8")
    return str(path)


def nidBlock(name="nidActive", nid=4.0e-5, error=4.0e-6, entries=1000, nz=50):
    return BLOCK.format(name=name, nid=nid, error=error, entries=entries, nz=nz)


class TestReadNID(unittest.TestCase):
    def test_reads_several_modules(self):
        with tempfile.TemporaryDirectory() as d:
            f = writeFile(d, "a.csv", [nidBlock("nidActive", 1.0e-5, 1.0e-6, 100, 7),
                                       nidBlock("nidSubstrate", 2.0e-5, 2.0e-6, 100, 9)])
            modules = readNID(f)
        self.assertEqual(sorted(modules), ["nidActive", "nidSubstrate"])
        self.assertAlmostEqual(modules["nidActive"]["nid"], 1.0e-5)
        self.assertEqual(modules["nidSubstrate"]["non-zeros"], 9)
        self.assertEqual(modules["nidActive"]["unit"], "MeV/g")

    def test_ignores_dose_blocks(self):
        """A CRDS run scores TID and DDD into one file; take only the DDD."""
        with tempfile.TemporaryDirectory() as d:
            f = writeFile(d, "a.csv", [DOSE_BLOCK, nidBlock()])
            modules = readNID(f)
        self.assertEqual(list(modules), ["nidActive"])

    def test_tid_only_file_yields_nothing(self):
        with tempfile.TemporaryDirectory() as d:
            f = writeFile(d, "a.csv", [DOSE_BLOCK])
            self.assertEqual(readNID(f), {})

    def test_module_prefix_filters(self):
        with tempfile.TemporaryDirectory() as d:
            f = writeFile(d, "a.csv", [nidBlock("nidActive"), nidBlock("otherThing")])
            self.assertEqual(list(readNID(f, modulePrefix="nid")), ["nidActive"])


class TestTotalNID(unittest.TestCase):
    def pooled(self, files):
        with tempfile.TemporaryDirectory() as d:
            for i, blocks in enumerate(files):
                writeFile(d, "f{}.csv".format(i), blocks)
            return totalNID(d)

    def test_entry_weighted_not_file_weighted(self):
        """GRAS drops snapshots into Res/, so files differ in entry count."""
        results = self.pooled([
            [nidBlock(nid=2.0, error=0.2, entries=100, nz=10)],
            [nidBlock(nid=4.0, error=0.4, entries=300, nz=30)],
        ])
        # Entry weighted: (100*2 + 300*4)/400 = 3.5.  File weighted would be 3.0.
        self.assertAlmostEqual(results["nidActive"]["nid"], 3.5)
        self.assertEqual(results["nidActive"]["entries"], 400)

    def test_pooled_error_formula(self):
        results = self.pooled([
            [nidBlock(nid=2.0, error=0.2, entries=100, nz=10)],
            [nidBlock(nid=4.0, error=0.4, entries=300, nz=30)],
        ])
        expected = ((100 * 0.2) ** 2 + (300 * 0.4) ** 2) ** 0.5 / 400
        self.assertAlmostEqual(results["nidActive"]["error"], expected)

    def test_birge_flags_inconsistent_files(self):
        results = self.pooled([
            [nidBlock(nid=1.0, error=1e-6, entries=100, nz=10)],
            [nidBlock(nid=9.0, error=1e-6, entries=100, nz=10)],
        ])
        self.assertTrue(results["nidActive"]["inconsistent"])

    def test_consistent_files_are_not_flagged(self):
        results = self.pooled([
            [nidBlock(nid=2.0, error=1.0, entries=100, nz=10)],
            [nidBlock(nid=2.1, error=1.0, entries=100, nz=10)],
        ])
        self.assertFalse(results["nidActive"]["inconsistent"])

    def test_file_count_guard(self):
        with tempfile.TemporaryDirectory() as d:
            writeFile(d, "a.csv", [nidBlock()])
            with self.assertRaises(ValueError):
                totalNID(d, expectedFiles=100)

    def test_module_count_guard(self):
        with tempfile.TemporaryDirectory() as d:
            writeFile(d, "a.csv", [nidBlock()])
            with self.assertRaises(ValueError):
                totalNID(d, expectedModules=2)

    def test_module_mismatch_between_files(self):
        with self.assertRaises(ValueError):
            self.pooled([[nidBlock("nidActive")], [nidBlock("nidSubstrate")]])

    def test_tid_only_directory_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            writeFile(d, "a.csv", [DOSE_BLOCK])
            with self.assertRaises(ValueError):
                totalNID(d)

    def test_empty_directory_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):
                totalNID(d)

    def test_particles_for_one_percent_scales_as_inverse_square(self):
        results = self.pooled([[nidBlock(nid=1.0, error=0.1, entries=1000, nz=50)]])
        # 10% now, so 1% needs 100x the primaries.
        self.assertAlmostEqual(results["nidActive"]["relative_error_percent"], 10.0)
        self.assertAlmostEqual(results["nidActive"]["particles_for_one_percent"], 1e5)

    def test_hit_fraction(self):
        results = self.pooled([[nidBlock(entries=1000, nz=50)]])
        self.assertAlmostEqual(results["nidActive"]["hit_fraction"], 0.05)


class TestModuleRecords(unittest.TestCase):
    def test_orders_by_the_given_list(self):
        with tempfile.TemporaryDirectory() as d:
            writeFile(d, "a.csv", [nidBlock("nidActive", 1.0e-5),
                                   nidBlock("nidSubstrate", 2.0e-5)])
            results = totalNID(d)
        records = moduleRecords(results, [
            ("nidSubstrate", "LED_die_0_PV", "substrate"),
            ("nidActive", "LED_active_0_PV", "active"),
        ])
        self.assertEqual([r["label"] for r in records], ["substrate", "active"])

    def test_missing_module_raises(self):
        with tempfile.TemporaryDirectory() as d:
            writeFile(d, "a.csv", [nidBlock("nidActive")])
            results = totalNID(d)
        with self.assertRaises(KeyError):
            moduleRecords(results, [("nidMissing", "v", "l")])


if __name__ == "__main__":
    unittest.main()