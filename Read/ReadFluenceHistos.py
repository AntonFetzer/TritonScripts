import numpy as np
import matplotlib.pyplot as plt
import csv
import os


def readFluenceHistos(file):
    """
    Read fluence histograms from a .csv file and return two histogram dicts

    Args:
        file: path to the .csv file containing two histogram blocks

    Returns:
        Tuple (ElectronHistDict, ProtonHistDict) where each is a dict with keys:
            'lower','upper','mean','value','error','entries' (numpy arrays)
    """
    keys = ['lower', 'upper', 'mean', 'value', 'error', 'entries']

    # Initialize dictionaries for Electron and Proton histograms
    ElectronHist = {key: [] for key in keys}
    ProtonHist = {key: [] for key in keys}

    # Read flags follow the same pattern as ReadDoseHistos:
    # 0 = looking for first "'Bin entries'" (start electrons)
    # 1 = reading electrons until "'End of Block'"
    # 2 = looking for second "'Bin entries'" (start protons)
    # 3 = reading protons until "'End of Block'"
    ReadFlag = 0

    with open(file, 'r') as f:
        reader = csv.reader(f)
        for line in reader:
            if ReadFlag == 0:
                if "'Bin entries'" in line:
                    ReadFlag = 1
            elif ReadFlag == 1:
                if "'End of Block'" in line:
                    ReadFlag = 2
                else:
                    for i, key in enumerate(keys):
                        ElectronHist[key].append(float(line[i]))
            elif ReadFlag == 2:
                if "'Bin entries'" in line:
                    ReadFlag = 3
            elif ReadFlag == 3:
                if "'End of Block'" in line:
                    ReadFlag = 4
                else:
                    for i, key in enumerate(keys):
                        ProtonHist[key].append(float(line[i]))

    # Check if any data was read
    if not ElectronHist['lower'] and not ProtonHist['lower']:
        raise ValueError("No histogram data found in file: {}".format(file))
    
    # Convert entries to integers
    ElectronHist['entries'] = [int(x) for x in ElectronHist['entries']]
    ProtonHist['entries'] = [int(x) for x in ProtonHist['entries']]

    # Convert to numpy arrays
    for key in ElectronHist:
        ElectronHist[key] = np.array(ElectronHist[key])
        ProtonHist[key] = np.array(ProtonHist[key])

    return ElectronHist, ProtonHist


if __name__ == "__main__":
    Path = "/l/triton_work/Fluence_Histograms/CarringtonShielded/Carrington-SEP-Expected-Int/1mm/Res/"
    # Find the first .csv file in the specified directory
    FullPath = next((os.path.join(Path, f) for f in os.listdir(Path) if f.endswith('.csv')), None)
    if FullPath is None:
        raise FileNotFoundError("No .csv file found in {}".format(Path))

    print("Reading in file:", FullPath)
    PlotDir = os.path.join(os.path.dirname(Path), '../')

    ElectronFluence, ProtonFluence = readFluenceHistos(FullPath)

    # Check that Electron and Proton histograms have data
    if ElectronFluence['lower'].size == 0:
        raise ValueError("No Electron histogram data found in file: {}".format(FullPath))
    if ProtonFluence['lower'].size == 0:
        raise ValueError("No Proton histogram data found in file: {}".format(FullPath))

    # Check that entries are integers
    if not np.issubdtype(ElectronFluence['entries'].dtype, np.integer):
        raise TypeError("Electron histogram 'entries' are not integers.")
    if not np.issubdtype(ProtonFluence['entries'].dtype, np.integer):
        raise TypeError("Proton histogram 'entries' are not integers.")

    plt.rc('axes', axisbelow=True)
    plt.figure(0)
    plt.bar(ElectronFluence['lower'], ElectronFluence['value'], width=ElectronFluence['upper'] - ElectronFluence['lower'], yerr=ElectronFluence['error'], align='edge')
    plt.xlabel('Energy [MeV]')
    plt.ylabel('Fluence [counts/cm2]')
    plt.title('Electron Fluence Spectrum')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(which='both')

    print("Saving Electron Histogram:", os.path.join(PlotDir, 'Electron Fluence Spectrum.pdf'))
    plt.savefig(os.path.join(PlotDir, 'Electron Fluence Spectrum.pdf'), format='pdf', bbox_inches='tight')


    plt.figure(1)
    plt.bar(ProtonFluence['lower'], ProtonFluence['value'], width=ProtonFluence['upper'] - ProtonFluence['lower'], yerr=ProtonFluence['error'], align='edge')
    plt.xlabel('Energy [MeV]')
    plt.ylabel('Fluence [counts/cm2]')
    plt.title('Proton Fluence Spectrum')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(which='both')

    print("Saving Proton Histogram:", os.path.join(PlotDir, 'Proton Fluence Spectrum.pdf'))
    plt.savefig(os.path.join(PlotDir, 'Proton Fluence Spectrum.pdf'), format='pdf', bbox_inches='tight')

    # plt.show()

