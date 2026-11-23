"""Regression test for Read.ReadNID and Dependencies.TotalNID.

What is pinned here is what a fresh displacement-damage aggregation gets wrong:
averaging by file instead of by entries, and pooling errors as a plain
quadrature sum.

Run from the repository root:

    python3 -m unittest discover -s Tests -p "Test*.py"
"""
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.TotalNID import totalNID  # noqa: E402
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


class TestTotalNID(unittest.TestCase):
    def pooled(self, files):
        with tempfile.TemporaryDirectory() as d:
            for i, blocks in enumerate(files):
                writeFile(d, "f{}.csv".format(i), blocks)
            return totalNID(d)

    def test_entry_weighted_mean_and_pooled_error(self):
        """GRAS drops snapshots into Res/, so files differ in entry count."""
        results = self.pooled([
            [nidBlock(nid=2.0, error=0.2, entries=100, nz=10)],
            [nidBlock(nid=4.0, error=0.4, entries=300, nz=30)],
        ])
        # Entry weighted: (100*2 + 300*4)/400 = 3.5.  File weighted would be 3.0.
        self.assertAlmostEqual(results["nidActive"]["nid"], 3.5)
        self.assertEqual(results["nidActive"]["entries"], 400)
        expected = ((100 * 0.2) ** 2 + (300 * 0.4) ** 2) ** 0.5 / 400
        self.assertAlmostEqual(results["nidActive"]["error"], expected)

    def test_birge_flags_only_inconsistent_files(self):
        results = self.pooled([
            [nidBlock(nid=1.0, error=1e-6, entries=100, nz=10)],
            [nidBlock(nid=9.0, error=1e-6, entries=100, nz=10)],
        ])
        self.assertTrue(results["nidActive"]["inconsistent"])
        results = self.pooled([
            [nidBlock(nid=2.0, error=1.0, entries=100, nz=10)],
            [nidBlock(nid=2.1, error=1.0, entries=100, nz=10)],
        ])
        self.assertFalse(results["nidActive"]["inconsistent"])

    def test_module_mismatch_between_files(self):
        with self.assertRaises(ValueError):
            self.pooled([[nidBlock("nidActive")], [nidBlock("nidSubstrate")]])


if __name__ == "__main__":
    unittest.main()
