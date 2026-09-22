"""Checks for the CRDS macro alias reader."""

import tempfile
import unittest
from pathlib import Path

from Read.ReadCRDSMacro import spectrumName


class TestCRDSMacro(unittest.TestCase):
    def test_spectrum_alias(self):
        with tempfile.TemporaryDirectory() as directory:
            Path(directory, "A.mac").write_text(
                "/control/alias Spectrum Clinical/TrueBeam/HUS-PDD-SSD100\n",
                encoding="utf-8")
            self.assertEqual(spectrumName(directory), "HUS-PDD-SSD100")

    def test_missing_spectrum_alias(self):
        with tempfile.TemporaryDirectory() as directory:
            Path(directory, "A.mac").write_text("/run/beamOn 10\n", encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "no Spectrum alias"):
                spectrumName(directory)
