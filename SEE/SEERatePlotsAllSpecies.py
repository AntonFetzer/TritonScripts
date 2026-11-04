import matplotlib.pyplot as plt
import numpy as np
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Dependencies.SEERates import readRateCSVs, writeRateCSV

BaseDirectory = "/scratch/work/fetzera1/GRAS/1Tile/LET_Histograms/"

# Generate list of all folders with All-Species in their name
FolderList = [BaseDirectory + f_ for f_ in os.listdir(BaseDirectory) if "All-Species" in f_ and os.path.isdir(BaseDirectory + f_)]

# Alternatively process only the folders given as command line arguments
if len(sys.argv) > 1:
    FolderList = [os.path.abspath(p) for p in sys.argv[1:]]

# The model names, such as GEO-CREME96, come from the folder names
Directories = {os.path.basename(p).replace("Carrington-", "").replace("-All-Species", ""): p + "/" for p in FolderList}

CrossectionNames = ["Cypress CY62167GE30-45ZXI", "NanoXplore SEU"]

# Define default errorbar style
default_errorbar_style = {
    'capsize': 4,        # Size of the error bar caps
    'elinewidth': 2,     # Width of the error bar lines
    'capthick': 2,       # Thickness of the error bar caps
}

FigureNumber = 0
for ModelName, Directory in Directories.items():
  for CrossectionName in CrossectionNames:
    FigureNumber += 1

    # Read the CSV file of every run in this set
    Dict = readRateCSVs([Directory + "*/* " + CrossectionName + " SEE-Rate.csv"], CrossectionName)

    # Every species must have the same shielding points, otherwise the total
    # over species would silently miss contributions (e.g. runs that produced
    # no LET histogram and were skipped by the rate-estimate script).
    Shieldings = {tuple(Dict[name]['Shielding']) for name in Dict}
    if len(Shieldings) > 1:
        expected = max(Shieldings, key=len)
        for name in sorted(Dict):
            missing = set(expected) - set(Dict[name]['Shielding'])
            if missing:
                print(f"{name}: missing shielding {sorted(str(m) for m in missing)}")
        raise ValueError(
            f"{ModelName}: species have inconsistent shielding points -- "
            "rerun the missing simulations and the rate-estimate script"
        )

    plt.figure(FigureNumber, figsize=(5, 7))

    # Sort the datasets by their SEE rate at 4mm shielding, so the legend
    # is ordered from highest to lowest contribution
    Shielding4mmIndex = np.where(Dict[list(Dict.keys())[0]]['Shielding'] == '4mm')[0][0]
    sorted_names = sorted(Dict.keys(), key=lambda name: Dict[name]['SEE_Rate'][Shielding4mmIndex], reverse=True)
    Top3Names = sorted_names[:3]
    RestNames = sorted_names[3:]

    # Combine all species into a total rate
    Total = {}
    Total['Shielding'] = Dict[sorted_names[0]]['Shielding']
    Total['SEE_Rate'] = np.sum([Dict[name]['SEE_Rate'] for name in sorted_names], axis=0)
    Total['SEE_Error'] = np.sqrt(np.sum([np.square(Dict[name]['SEE_Error']) for name in sorted_names], axis=0))

    # Combine all species except the top 3 at 4mm shielding
    Total_Ex3 = {}
    Total_Ex3['Shielding'] = Dict[sorted_names[0]]['Shielding']
    Total_Ex3['SEE_Rate'] = np.sum([Dict[name]['SEE_Rate'] for name in RestNames], axis=0)
    Total_Ex3['SEE_Error'] = np.sqrt(np.sum([np.square(Dict[name]['SEE_Error']) for name in RestNames], axis=0))

    # Combine all non-proton species into a total heavy-ion rate and write it
    # next to the data in the structure of the SEE-Rate CSV files, where
    # SEERatePlotsLET.py reads it.
    HeavyNames = [name for name in sorted_names if not name.endswith("Proton")]
    HeavyRate = np.sum([Dict[name]['SEE_Rate'] for name in HeavyNames], axis=0)
    HeavyError = np.sqrt(np.sum([np.square(Dict[name]['SEE_Error']) for name in HeavyNames], axis=0))
    HeavyEntries = np.sum([Dict[name]['Entries_Contributing_To_SEE'] for name in HeavyNames], axis=0)
    SummaryName = (sorted_names[0].rsplit("-", 1)[0]) + "-AllHeavyIons"

    Rows = []
    for i, shielding in enumerate(Total['Shielding']):
        RelativeError = HeavyError[i] / HeavyRate[i] if HeavyRate[i] else float('nan')
        Rows.append((shielding, HeavyRate[i], HeavyError[i], RelativeError, int(HeavyEntries[i])))
    writeRateCSV(Directory + SummaryName + " " + CrossectionName + " SEE-Rate.csv", SummaryName, CrossectionName, Rows)

    plt.errorbar(Total['Shielding'], Total['SEE_Rate'], yerr=Total['SEE_Error'],
                 label="Total", color='black', linewidth=3, **default_errorbar_style)

    for i, name in enumerate(Top3Names):
        plt.errorbar(Dict[name]['Shielding'], Dict[name]['SEE_Rate'],
                     yerr=Dict[name]['SEE_Error'], label=name,
                     color='C' + str(i % 10), **default_errorbar_style)

    plt.errorbar(Total_Ex3['Shielding'], Total_Ex3['SEE_Rate'], yerr=Total_Ex3['SEE_Error'],
                 label="Total excluding top 3", color='C5', linestyle='--', **default_errorbar_style)

    plt.yscale("log")
    plt.grid()
    plt.title(CrossectionName + "\n" + ModelName.replace("-", " ") + " SEU Rate Estimates - All Species")

    plt.xlabel("Aluminium Shielding Thickness")
    plt.ylabel("Single Event Upset Rate [s-1 bit-1]")

    plt.legend(loc='lower left')

    plt.savefig(Directory + CrossectionName + " Rates All Species.pdf", format='pdf', bbox_inches="tight")
    # plt.show()
