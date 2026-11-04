import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat, ufloat_fromstr


REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.TotalDose import totalDose
from Dependencies.MergeTotalDose import mergeTotalDose
from Plotting.TilePlots import plotDose, plotNonZeros, plotRelativeError

# First run TotalDose on all species folders of the AllSpecies folder provided as commandline argument
# to get the total dose for each species, then run MergeTotalDose 
# to merge the total dose for all species into a single file, then plot the shielding curves 
# for all species in a single plot with the total as well. 
# Store the merged total dose result in a CSV file with the same format as the individual species total dose files. 
#
# Run as: python3 Plotting/AllSpeciesShieldingCurves.py <AllSpecies folder>
# Doses stay in the normalisation the runs were set up with, as in the individual species files.
# Each species folder also gets its TotalDose_ CSV and its Dose_, Error_ and NonZeros_ plots.
# The proton run is left out of the total and the plot, because the complete cosmic and solar proton runs in
# Shielding_Curves/Carrington replace it. The total and the plot are therefore named Heavy-Ion, not All-Species.

# Species drawn in colour, because they carry most of the dose. Each keeps its colour in the plot of every
# orbit, and the order passes the colour-blindness check for line charts. All other species are grey context.
Highlighted = {
    "Helium": "#eb6834",
    "Carbon": "#1baf7a",
    "Oxygen": "#eda100",
    "Silicon": "#e87ba4",
    "Iron": "#008300",
}
OtherColour = "#c3c2b7"


def _writeTotalDose(Results, OutputFile):
    """Write a pooled dose result as a TotalDose CSV, rounded to two significant digits of the error."""
    with open(OutputFile, 'w') as f:
        f.write("Tile, Dose [kRad], Error [kRad], Non-Zero Entries\n")
        for i in range(len(Results['dose'])):
            Rounded = ufloat_fromstr(f"{ufloat(Results['dose'][i], Results['error'][i]):.2u}")
            f.write(f"{i},{Rounded.n:.12g},{Rounded.s:.12g},{Results['non-zeros'][i]}\n")
    print("Wrote", OutputFile)


if len(sys.argv) != 2:
    sys.exit("Usage: python3 Plotting/AllSpeciesShieldingCurves.py <AllSpecies folder>")

Folder = Path(sys.argv[1]).resolve()

Species = {}
Missing = []
for SpeciesFolder in sorted(d for d in Folder.iterdir() if (d / "Res").is_dir()):
    if any((SpeciesFolder / "Res").glob("*.csv")):
        Name = SpeciesFolder.name
        Results = totalDose(str(SpeciesFolder / "Res"))
        Species[Name.split("-")[-1]] = Results
        _writeTotalDose(Results, SpeciesFolder / ("TotalDose_" + Name + ".csv"))
        plotDose(Results['dose'], Results['error'], 'Dose per tile ' + Name, SpeciesFolder / ("Dose_" + Name + ".pdf"))
        plotRelativeError(Results['dose'], Results['error'], SpeciesFolder / ("Error_" + Name + ".pdf"))
        plotNonZeros(Results['non-zeros'], SpeciesFolder / ("NonZeros_" + Name + ".pdf"))
    else:
        Missing.append(SpeciesFolder.name)

HeavyIons = {s: d for s, d in Species.items() if s != "Proton"}
if not HeavyIons:
    sys.exit("ERROR !!! No heavy-ion species folder with results in: " + str(Folder))
if Missing:
    print("Warning !!! The total leaves out", len(Missing), "species without results:", ", ".join(Missing))

Total = mergeTotalDose(list(HeavyIons.values()))
HeavyName = Folder.name.replace("All-Species", "Heavy-Ion")

NumTiles = len(Total['dose'])
x = np.linspace(0, 10, num=NumTiles, endpoint=True)  # The 100Al10mm geometry: 0 to 10 mm of aluminium

# Store the merged result, rounded as in the individual species files
_writeTotalDose(Total, Folder / ("TotalDose_" + HeavyName + ".csv"))

# Plot every heavy-ion species and the total. Tiles without dose are left out of the logarithmic axis.
plt.figure(figsize=(7.5, 4.5))
Others = [s for s in HeavyIons if s not in Highlighted]
for i, s in enumerate(Others):
    plt.plot(x, np.ma.masked_less_equal(Species[s]['dose'], 0), color=OtherColour, linewidth=0.5,
             label=f"Other species ({len(Others)})" if i == 0 else None)
for s, Colour in Highlighted.items():
    if s in Species:
        plt.errorbar(x, np.ma.masked_less_equal(Species[s]['dose'], 0), Species[s]['error'], color=Colour,
                     linewidth=1.5, elinewidth=0.8, label=s, zorder=3)
plt.errorbar(x, np.ma.masked_less_equal(Total['dose'], 0), Total['error'], color='k', linewidth=2.5, elinewidth=0.8,
             label="Total")

plt.yscale("log")
plt.xlim(x[0], x[-1])
plt.grid(color="#e6e5e0", linewidth=0.5)
plt.gca().set_axisbelow(True)
plt.title(HeavyName)
plt.xlabel("Aluminium Thickness [mm]")
plt.ylabel("Ionising Dose [kRad per normalisation fluence]")

# Legend order: the total, the highlighted species, then the rest
Handles, Labels = plt.gca().get_legend_handles_labels()
Order = sorted(range(len(Labels)), key=lambda i: (Labels[i] != "Total", Labels[i].startswith("Other")))
plt.legend([Handles[i] for i in Order], [Labels[i] for i in Order], loc="center left", bbox_to_anchor=(1.01, 0.5),
           frameon=False)

OutputPlot = Folder / ("ShieldingCurves_" + HeavyName + ".pdf")
plt.savefig(OutputPlot, format='pdf', bbox_inches="tight")
print("Wrote", OutputPlot)
