from os import path
import pandas as pd
import matplotlib.pyplot as plt
from GRAS.Read.ReadSourceHistos import readSourceHistos
from glob import iglob

Paths = [
    "/l/triton_work/SourceHistograms/ISO-GTO/ISO-GTO-Fe/Res/",
    "/l/triton_work/SourceHistograms/ISO-GTO/ISO-GTO-Fe-mission/Res/"
]

# Find the first .csv file in each of the paths
Files = [next(iglob(path.join(Path, '*.csv')), None) for Path in Paths]

# List to store the histograms
all_histograms = []

# Read histograms using the provided function
for file in Files:
    histograms = readSourceHistos(file)
    all_histograms.append(histograms)

# Plotting Energy Histograms for Comparison
plt.figure()

for idx, histograms in enumerate(all_histograms):
    energy_data = histograms['Energy']

    # Derive the label from the folder name in which the file is located
    label = path.basename(path.dirname(path.dirname(Files[idx])))

    plt.bar(energy_data['lower'], energy_data['value'], width=energy_data['upper'] - energy_data['lower'],
            yerr=energy_data['error'], align='edge', alpha=0.5, label=label)

plt.xlabel('Energy [MeV]')
plt.ylabel('Counts per bin')
plt.title('Energy Histogram Comparison')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.show()
