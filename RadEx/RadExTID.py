"""Total ionising dose of RadEx GRAS runs, written into each run folder.

For every folder below PATH that contains a Res/ folder, pool the GRAS result
files with Dependencies.TotalDose.totalDose and write TotalDose_<folder>.csv
next to Res/, with the diagnostic plots Dose_, Error_ and NonZeros_<folder>.pdf.
Every tally is written in the order the detector macro defines them (for the
VT01 module, the twelve gate oxides, then the twelve silicon dies). Doses are
in kRad per unit of the normalisation fluence the run was set up with, as GRAS
writes them; any rescaling to a delivered fluence is run-specific
post-processing and does not belong here.

Run from anywhere:

    python3 RadEx/RadExTID.py <run folder, or a folder of run folders>
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from uncertainties import ufloat, ufloat_fromstr

from Dependencies.TotalDose import totalDose
from Plotting.TilePlots import plotDose, plotNonZeros, plotRelativeError

parser = argparse.ArgumentParser(description="Total ionising dose of RadEx GRAS runs.")
parser.add_argument("path", help="run folder, or a folder searched recursively for run folders")
args = parser.parse_args()

for root, dirs, files in sorted(os.walk(args.path)):
    if 'Res' not in dirs:
        continue
    folder_name = os.path.basename(os.path.normpath(root))
    Results = totalDose(os.path.join(root, 'Res'))

    # Round each dose to two significant digits of its error, as the uncertainties
    # package displays it, and read the rounded value back from that string.
    output_file = os.path.join(root, "TotalDose_" + folder_name + ".csv")
    with open(output_file, 'w') as f:
        f.write("Tile, Dose [kRad], Error [kRad], Non-Zero Entries\n")
        for i in range(len(Results['dose'])):
            Rounded = ufloat_fromstr(f"{ufloat(Results['dose'][i], Results['error'][i]):.2u}")
            f.write(f"{i},{Rounded.n:.12g},{Rounded.s:.12g},{Results['non-zeros'][i]}\n")
    print("Wrote", output_file)

    plotDose(Results['dose'], Results['error'], 'Dose per tile ' + folder_name,
             os.path.join(root, "Dose_" + folder_name + ".pdf"))
    plotRelativeError(Results['dose'], Results['error'], os.path.join(root, "Error_" + folder_name + ".pdf"))
    plotNonZeros(Results['non-zeros'], os.path.join(root, "NonZeros_" + folder_name + ".pdf"))
