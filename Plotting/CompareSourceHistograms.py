import pandas as pd
import matplotlib.pyplot as plt
from GRAS.Dependencies.TotalSourceHistograms import totalSourceHistos
import os

# Set the grid to be behind the plot
plt.rc('axes', axisbelow=True)

Paths = [
    "/l/triton_work/Source_Histograms/Carrington/CarringtonElectronINTEGRALPowTabelated/Res",
    "/l/triton_work/Source_Histograms/Carrington/CarringtonElectronDiffPowTabelated/Res",
    "/l/triton_work/Source_Histograms/Carrington/CarringtonElectronDiffPow/Res",

    "/l/triton_work/Source_Histograms/Carrington/Van-Allen-Probes-AE9-mission/Res",

    # "/l/triton_work/Source_Histograms/Carrington/Carrington-SEP-Expected-Int-With0/Res",
    # "/l/triton_work/Source_Histograms/Carrington/Carrington-SEP-Plus2Sigma-Int-With0/Res",
    # "/l/triton_work/Source_Histograms/Carrington/Carrington-SEP-Minus2Sigma-Int-With0/Res",

    # "/l/triton_work/Source_Histograms/Carrington/Carrington-SEP-Expected-Diff/Res",
    # "/l/triton_work/Source_Histograms/Carrington/Carrington-SEP-Expected-Int/Res",
]

# List to store the histograms
all_histograms = []

# Read histograms using the provided function
for path in Paths:
    histograms = totalSourceHistos(path)
    all_histograms.append(histograms)

NumberOfHistograms = len(all_histograms)

# Plotting Energy Histograms for Comparison
plt.figure(0)

for i in range(NumberOfHistograms):
    histograms = all_histograms[i]
    energy = histograms['Energy']

    # Derive the label from the second to last folder name in the path
    label = os.path.basename(os.path.dirname(Paths[i]))

    # Plot histogram without error bar
    # plt.bar(energy['lower'], energy['value'], width=energy['upper'] - energy['lower'], align='edge', alpha=0.5, label=label)

    # Plot histogram as error bars
    plt.errorbar(energy['mean'], energy['value'], energy['error'], fmt=' ', capsize=5, elinewidth=1, capthick=1, label=label)

    # Plot histogram with error bars
    # plt.bar(energy['lower'], energy['value'], width=energy['upper'] - energy['lower'], yerr=energy['error'], align='edge')

    # Plot histograms as scatter plot
    # plt.plot(energy['mean'], energy['value'], '.-', label=label)

plt.xlabel('Energy [MeV]')
plt.ylabel('Counts per bin')
plt.title('Energy Histogram Comparison')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.grid(which='both')

# Save the plot to a file
plt.savefig("/l/triton_work/Source_Histograms/Carrington/Energy_Histogram_Comparison.pdf", format='pdf', bbox_inches='tight')

# plt.show()