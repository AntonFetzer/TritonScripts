"""Read.ReadLETHistos and Read.ReadFluenceHistos: values, units and block order of GRAS histogram files."""
import tempfile
import unittest
from pathlib import Path

import numpy as np

from Read.ReadFluenceHistos import readFluenceHistos
from Read.ReadLETHistos import readLETHistos

# Shortened GRAS output, in the block order of real files: a scalar block, two histograms, the 'general' block.
SCALAR = """'*', 9, 0, 4, 0, 4, 4, 1, 3
'GRAS_DATA_TYPE',   -1,'STAT_DOUBLE'
'Entries','',    1,'Number of entries'
  0.0023786,   1.903e-05,  5.2189e+06,  5.2189e+06
'End of Block'
"""

HISTOGRAM = """'*', 32, 1, 24, 0, 6, 6, 2, 2
'GRAS HISTOGRAM 1D'
'OVERFLOW_ENTRIES',    1,     0, 'none'
'lower','{x}',    1,'Bin lower edge'
'upper','{x}',    1,'Bin upper edge'
'mean','{x}',    1,'Bin mean'
'value','{y}',    1,'Bin value'
'error','{y}',    1,'Bin error'
'entries','',    1,'Bin entries'
{rows}
'End of Block'
"""

GENERAL = """'*', 12, 0, 2, 0, 9, 9, 1, 0
'GRAS_MODULE_NAME',   -1,'general'
'NumOfEvt','',    1,'Number of events'
'Gamma','',    1,'Average number of gammas'
      20000,      96.072
'End of Block'
'End of File'
"""

FIRST = """       2330,        4660,        2796,       5e+06,       1e+06,          25
       4660,        9320,           0,           0,           0,           0"""
SECOND = """       4660,       23300,        6990,       3e+06,       2e+06,           7"""


def read(reader, x, y):
    text = SCALAR + HISTOGRAM.format(x=x, y=y, rows=FIRST) + HISTOGRAM.format(x=x, y=y, rows=SECOND) + GENERAL
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "histos.csv"
        path.write_text(text, encoding="ascii")
        return reader(path)


class TestReadHistos(unittest.TestCase):

    def check(self, histogram, expected):
        for key, values in expected.items():
            np.testing.assert_allclose(histogram[key], values, err_msg=key)
        self.assertTrue(np.issubdtype(histogram['entries'].dtype, np.integer))

    def test_let_is_converted_to_silicon_mass_let_per_cm2(self):
        # MeV/cm / (2330 mg/cm3) to MeV cm2 mg-1, counts / (1e6 cm2) to cm-2
        LET, Eff = read(readLETHistos, "MeV/cm", "counts")
        self.check(LET, {'lower': [1.0, 2.0], 'upper': [2.0, 4.0], 'mean': [1.2, 0.0],
                         'value': [5.0, 0.0], 'error': [1.0, 0.0], 'entries': [25, 0]})
        self.check(Eff, {'lower': [2.0], 'upper': [10.0], 'mean': [3.0],
                         'value': [3.0], 'error': [2.0], 'entries': [7]})
        self.assertEqual(LET['events'], 20000.0)
        self.assertEqual(Eff['events'], 20000.0)

    def test_fluence_keeps_its_units_electrons_before_protons(self):
        Electron, Proton = read(readFluenceHistos, "MeV", "counts/cm2")
        self.check(Electron, {'lower': [2330.0, 4660.0], 'upper': [4660.0, 9320.0], 'mean': [2796.0, 0.0],
                              'value': [5e6, 0.0], 'error': [1e6, 0.0], 'entries': [25, 0]})
        self.check(Proton, {'lower': [4660.0], 'upper': [23300.0], 'mean': [6990.0],
                            'value': [3e6], 'error': [2e6], 'entries': [7]})
        self.assertEqual(Electron['events'], 20000.0)
        self.assertEqual(Proton['events'], 20000.0)


if __name__ == "__main__":
    unittest.main()
