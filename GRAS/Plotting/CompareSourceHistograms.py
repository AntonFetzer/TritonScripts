from os import path
import pandas as pd
import matplotlib.pyplot as plt
from GRAS.Read.ReadSourceHistos import readSourceHistos
from glob import iglob

Paths = [
    "/l/triton_work/SourceHistograms/Carrington/CarringtonElectronINTEGRALPowTabelated/Res",
    "/l/triton_work/SourceHistograms/Carrington/CarringtonElectronDiffPowTabelated/Res",
    "/l/triton_work/SourceHistograms/A9-GTO/AE9/Res",
    "/l/triton_work/SourceHistograms/A9-FS1/AE9/Res",
]

# Find the first .csv file in each of the paths
Files = [next(iglob(path.join(Path, '*.csv')), None) for Path in Paths]

# Check if all files were found
if None in Files:
    raise ValueError("Could not locate a .csv file in one or more of the specified directories.")

# List to store the histograms
all_histograms = []

# Read histograms using the provided function
for file in Files:
    histograms = readSourceHistos(file)
    all_histograms.append(histograms)

# Plotting Energy Histograms for Comparison
plt.figure()

for idx, histograms in enumerate(all_histograms):
    energy = histograms['Energy']

    # Derive the label from the folder name in which the file is located
    label = path.basename(path.dirname(path.dirname(Files[idx])))

    #plt.errorbar(energy['lower'], energy['value'], width=energy['upper'] - energy['lower'], yerr=energy['error'], align='edge', alpha=0.5, label=label)
    plt.errorbar(energy['mean'], energy['value'], energy['error'], fmt=' ', capsize=5, elinewidth=1, capthick=1, label=label)

plt.xlabel('Energy [MeV]')
plt.ylabel('Counts per bin')
plt.title('Energy Histogram Comparison')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.grid(which='both')

#plt.show()

# Save the plot to a file
plt.savefig("/l/triton_work/SourceHistograms/Carrington/Energy_Histogram_Comparison.pdf", format='pdf', bbox_inches='tight')
