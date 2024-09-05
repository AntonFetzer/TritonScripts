import matplotlib.pyplot as plt
import numpy as np
from GRAS.Dependencies.TotalFluenceHistograms import totalFluenceHistos
import os

# Set the grid to be behind the plot
#plt.rc('axes', axisbelow=True)

Paths = [
    "/l/triton_work/Fluence_Histograms/Carrington/CarringtonElectronINTEGRALPowTabelated/Res",
    "/l/triton_work/Fluence_Histograms/Carrington/CarringtonElectronDiffPowTabelated/Res",
    "/l/triton_work/Fluence_Histograms/Carrington/CarringtonElectronDiffPow/Res",

    "/l/triton_work/Fluence_Histograms/Carrington/Van-Allen-Belt-Probes-AE9-mission/Res",

    "/l/triton_work/Fluence_Histograms/Carrington/Carrington-SEP-Expected-Int-With0/Res",
    "/l/triton_work/Fluence_Histograms/Carrington/Carrington-SEP-Plus2Sigma-Int-With0/Res",
    "/l/triton_work/Fluence_Histograms/Carrington/Carrington-SEP-Minus2Sigma-Int-With0/Res",

    # "/l/triton_work/Fluence_Histograms/Carrington/Carrington-SEP-Expected-Diff/Res",
    # "/l/triton_work/Fluence_Histograms/Carrington/Carrington-SEP-Expected-Int/Res",

    "/l/triton_work/Fluence_Histograms/Carrington/Van-Allen-Belt-Probes-AP9-mission/Res",

    "/l/triton_work/Fluence_Histograms/Carrington/GEO-SolarProton-5minPeakFlux/Res",
]

# List to store the histograms
all_histograms = []

# Read histograms using the provided function
for path in Paths:
    histograms = totalFluenceHistos(path)
    all_histograms.append(histograms)

NumberOfHistograms = len(all_histograms)


# Plotting Electron Histograms
plt.figure(0)

for i in range(NumberOfHistograms):
    histograms = all_histograms[i]
    Elec = histograms['Electrons']


    if Elec['value'].sum() != 0:
        # Derive the label from the second to last folder name in the path
        label = os.path.basename(os.path.dirname(Paths[i]))

        if label == 'Van-Allen-Belt-Probes-AE9-mission' or label == 'Van-Allen-Belt-Probes-AP9-mission':
            # Convert from monthly fluence to flux per second
            Elec['value'] /= 30 * 24 * 60 * 60  
            Elec['error'] /= 30 * 24 * 60 * 60

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
plt.savefig("/l/triton_work/Fluence_Histograms/Carrington/Electrons_Histogram_Comparison.pdf", format='pdf', bbox_inches='tight')


# Plotting Proton Histograms
plt.figure(1)

for i in range(NumberOfHistograms):
    histograms = all_histograms[i]
    Prot = histograms['Protons']

    # Derive the label from the second to last folder name in the path
    label = os.path.basename(os.path.dirname(Paths[i]))

    if label == 'Van-Allen-Belt-Probes-AE9-mission' or label == 'Van-Allen-Belt-Probes-AP9-mission':
        # Convert from monthly fluence to flux per second
        Prot['value'] /= 30 * 24 * 60 * 60  
        Prot['error'] /= 30 * 24 * 60 * 60

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
plt.savefig("/l/triton_work/Fluence_Histograms/Carrington/Protons_Histogram_Comparison.pdf", format='pdf', bbox_inches='tight')
# plt.show()