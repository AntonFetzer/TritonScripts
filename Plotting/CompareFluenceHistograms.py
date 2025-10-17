import matplotlib.pyplot as plt
import numpy as np
from Dependencies.TotalFluenceHistograms import totalFluenceHistos
import os

# Set the grid to be behind the plot
#plt.rc('axes', axisbelow=True)
Path = "/l/triton_work/Fluence_Histograms/Carrington/"

Folders = [
    
    # "CarringtonElectronDiffPowTabelated/Res",
    # "CarringtonElectronINTEGRALPowTabelated/Res",

    # "Carrington-SEP-Expected-Diff/Res",
    # "Carrington-SEP-Expected-Int/Res",
    "Carrington-SEP-Expected-Int-With0/Res",
    # "Carrington-SEP-Plus2Sigma-Int-With0/Res",
    # "Carrington-SEP-Minus2Sigma-Int-With0/Res",

    # "GEO-cosmic-proton/Res",
    # "GEO-electron/Res",
    # "GEO-solar-proton/Res",
    "GEO-trapped-proton/Res",

    # "MEO-cosmic-proton/Res",
    # "MEO-electron/Res",
    # "MEO-solar-proton/Res",
    "MEO-trapped-proton/Res",

    # "VAP-cosmic-proton/Res",
    # "VAP-electron/Res",
    # "VAP-electron_Old/Res",
    # "VAP-solar-proton/Res",
    "VAP-trapped-proton/Res",
    "VAP-trapped-proton_Old/Res",

    # "SolarProton-5minPeakFlux/Res",

    # "LEO-cosmic-proton/Res",
    # "LEO-electron/Res",
    # "LEO-solar-proton-NoZero/Res",
    "LEO-trapped-proton/Res",
]

# List to store the histograms
all_histograms = []

# Read histograms using the provided function
for folder in Folders:
    histograms = totalFluenceHistos(Path + folder)
    all_histograms.append(histograms)

NumberOfHistograms = len(all_histograms)


# Plotting Electron Histograms
plt.figure(0)

for i in range(NumberOfHistograms):
    histograms = all_histograms[i]
    Elec = histograms['Electrons']

    if sum(Elec['value']) == 0:
        continue

    # Derive the label from the second to last folder name in the path
    label = Folders[i].split('/')[0]

    if "GEO" in label or "VAP" in label or "LEO" in label or "MEO" in label:
        # Convert 11 year mission fluence to flux per second
        Elec['value'] /= (4015 * 24 * 60 * 60)
        Elec['error'] /= (4015 * 24 * 60 * 60)

    if "Old" in label:
        Elec['value'] *= (4015/30)
        Elec['error'] *= (4015/30)

    # Replace 0 with NaN using numpy
    Elec['value'] = np.where(Elec['value'] == 0, np.nan, Elec['value'])

    # Plot histogram without error bar
    # plt.bar(Elec['lower'], Elec['value'], width=Elec['upper'] - Elec['lower'], align='edge', alpha=0.5, label=label)

    # Plot histogram as error bars
    # plt.errorbar(Elec['mean'], Elec['value'], Elec['error'], fmt=' ', capsize=2, elinewidth=1, capthick=1, label=label)

    # Plot histogram with error bars
    # plt.bar(Elec['lower'], Elec['value'], width=Elec['upper'] - Elec['lower'], yerr=Elec['error'], align='edge')

    # Plot histograms as line plot
    # plt.plot(Elec['mean'], Elec['value'], label=label)

    # Plot histograms as line plot with error bars
    plt.errorbar(Elec['mean'], Elec['value'], Elec['error'], capsize=2, elinewidth=1, capthick=1, label=label)

plt.xlabel('Energy [MeV]')
plt.ylabel('Fluence [counts/cm2]')
plt.title('Electron Fluence Spectrum')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.grid(which='both')

# Save the plot to a file
plt.savefig(Path + "Electrons_Histogram_Comparison.pdf", format='pdf', bbox_inches='tight')


# Plotting Proton Histograms
plt.figure(1)

for i in range(NumberOfHistograms):
    histograms = all_histograms[i]
    Prot = histograms['Protons']

    if sum(Prot['value']) == 0:
        continue

    # Derive the label from the Folder name
    label = Folders[i].split('/')[0]

    if "GEO" in label or "VAP" in label or "LEO" in label or "MEO" in label:
        # Convert 11 year mission fluence to flux per second
        Prot['value'] /= 4015 * 24 * 60 * 60  
        Prot['error'] /= 4015 * 24 * 60 * 60

    if "Old" in label:
        Prot['value'] *= (4015/30)
        Prot['error'] *= (4015/30)

    # Replace 0 with NaN using numpy
    Prot['value'] = np.where(Prot['value'] == 0, np.nan, Prot['value'])

    # Plot histogram without error bar
    # plt.bar(Prot['lower'], Prot['value'], width=Prot['upper'] - Prot['lower'], align='edge', alpha=0.5, label=label)

    # Plot histogram as error bars
    plt.errorbar(Prot['mean'], Prot['value'], Prot['error'], fmt='', capsize=2, elinewidth=1, capthick=1, label=label)

    # Plot histogram with error bars
    # plt.bar(Prot['lower'], Prot['value'], width=Prot['upper'] - Prot['lower'], yerr=Prot['error'], align='edge')

    # Plot histograms as scatter plot
    # plt.plot(Prot['mean'], Prot['value'], '.-', label=label)

plt.xlabel('Energy [MeV]')
plt.ylabel('Fluence [counts/cm2]')
plt.title('Proton Fluence Spectrum')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.grid(which='both')

# Save the plot to a file
plt.savefig(Path + "Protons_Histogram_Comparison.pdf", format='pdf', bbox_inches='tight')
# plt.show()