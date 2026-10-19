"""Read.ReadGPSSources: components of a multi-source mono-energetic GPS macro."""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from Read.ReadGPSSources import readGPSSources

MACRO = """# spectrum
/gps/source/clear
/gps/source/add 3.0
/gps/particle e-
/gps/ene/type Mono
/gps/ene/mono 0.25 MeV
/control/execute ../{BeamMacro}.mac

/gps/source/add 1.0   # second
/gps/particle e-
/gps/ene/type Mono
/gps/energy 500 keV
/gps/source/list
"""


class ReadGPSSources(unittest.TestCase):
    def read(self, text):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "spectrum.mac"
            path.write_text(text, encoding="utf-8")
            return readGPSSources(path)

    def test_energies_and_fractions(self):
        gps = self.read(MACRO)
        np.testing.assert_allclose(gps["energy"], [0.25, 0.5])
        np.testing.assert_allclose(gps["intensity"], [3.0, 1.0])
        np.testing.assert_allclose(gps["fraction"], [0.75, 0.25])
        self.assertEqual(gps["particle"], ["e-", "e-"])

    def test_source_without_energy(self):
        with self.assertRaises(ValueError):
            self.read("/gps/source/add 1\n/gps/particle e-\n")

    def test_non_mono_source(self):
        with self.assertRaises(ValueError):
            self.read("/gps/source/add 1\n/gps/ene/type Arb\n")


if __name__ == "__main__":
    unittest.main()
