"""Total ionising dose and displacement damage of CRDS GRAS runs, written into each run folder.

CRDS scores the RadFET TID tallies and the LED NID modules in the same result
files. For every folder below PATH that contains a Res/ folder, pool the TID
with Dependencies.TotalDose.totalDose and the NID with
Dependencies.TotalNID.totalNID, and write TotalDose_<folder>.csv and
TotalNID_<folder>.csv next to Res/, with the TID plots Dose_, Error_ and
NonZeros_<folder>.pdf. Values are per unit of the normalisation fluence the run
was set up with, in the units GRAS writes (kRad, MeV/g); any rescaling to a
delivered fluence is run-specific post-processing and does not belong here.

Run from anywhere:

    python3 CRDS/CRDSDose.py <run folder, or a folder of run folders>
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from uncertainties import ufloat, ufloat_fromstr

from Dependencies.TotalDose import totalDose
from Dependencies.TotalNID import totalNID
from Plotting.TilePlots import plotDose, plotNonZeros, plotRelativeError
from Read.ReadNID import readNID


def rounded(value, error):
    """Round to two significant digits of the error, as the uncertainties package displays it."""
    result = ufloat_fromstr(f"{ufloat(value, error):.2u}")
    return f"{result.n:.12g}", f"{result.s:.12g}"


parser = argparse.ArgumentParser(description="Total TID and NID of CRDS GRAS runs.")
parser.add_argument("path", help="run folder, or a folder searched recursively for run folders")
args = parser.parse_args()

for root, dirs, files in sorted(os.walk(args.path)):
    if 'Res' not in dirs:
        continue
    folder_name = os.path.basename(os.path.normpath(root))
    res = os.path.join(root, 'Res')

    Results = totalDose(res)
    output_file = os.path.join(root, "TotalDose_" + folder_name + ".csv")
    with open(output_file, 'w') as f:
        f.write("Tile, Dose [kRad], Error [kRad], Non-Zero Entries\n")
        for i in range(len(Results['dose'])):
            dose, error = rounded(Results['dose'][i], Results['error'][i])
            f.write(f"{i},{dose},{error},{Results['non-zeros'][i]}\n")
    print("Wrote", output_file)
    plotDose(Results['dose'], Results['error'], 'Dose per tile ' + folder_name,
             os.path.join(root, "Dose_" + folder_name + ".pdf"))
    plotRelativeError(Results['dose'], Results['error'], os.path.join(root, "Error_" + folder_name + ".pdf"))
    plotNonZeros(Results['non-zeros'], os.path.join(root, "NonZeros_" + folder_name + ".pdf"))

    first = sorted(f for f in os.listdir(res) if f.endswith(".csv"))[0]
    if not readNID(os.path.join(res, first)):
        print("No NID modules in", res)
        continue
    NID = totalNID(res)
    output_file = os.path.join(root, "TotalNID_" + folder_name + ".csv")
    with open(output_file, 'w') as f:
        f.write("Module, NID, Error, Unit, Non-Zero Entries, Birge Ratio\n")
        for module in sorted(NID):
            entry = NID[module]
            nid, error = rounded(entry['nid'], entry['error'])
            f.write(f"{module},{nid},{error},{entry['unit']},{entry['non-zeros']},{entry['birge_ratio']:.2f}\n")
    print("Wrote", output_file)
