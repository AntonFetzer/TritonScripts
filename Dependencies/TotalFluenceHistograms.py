import matplotlib.pyplot as plt
from Read.ReadFluenceHistos import readFluenceHistos
from Dependencies.MergeHistograms import mergeHistograms
import sys
import os

def totalFluenceHistos(path):
    """
    Accumulates and processes GRAS Fluence Histograms from multiple paralel GRAS Result files.

    Args:
        path (str): The path to the Res folder containing the csv files.

    Returns:
        dict: A dictionary containing the accumulated Fluence histos.

    Raises:
        SystemExit: If no files are found in the specified path.
        SystemExit: If the histogram bins in a file do not match the bins in the first file.
    """
    print("\nReading in all csv files in folder:", path)

    # Get list of all csv files in Path
    Files = [f for f in os.listdir(path) if "Fluence" in f and f.endswith(".csv")]
    
    if not Files:
        print("ERROR !!! No files found")
        # Wait for user input before exiting
        input("Press Enter to continue")
        return None
    
    NumFiles = len(Files)
    # print("Number of Files:", NumFiles)

    # Initialise the histogram lists
    ElectronHistos = []
    ProtonHistos = []

    for File in Files:
        ElectronHist, ProtonHist = readFluenceHistos(os.path.join(path, File))
        ElectronHistos.append(ElectronHist)
        ProtonHistos.append(ProtonHist)

    if NumFiles == 1:
        print("Only one file found, no merging required.")
        return ElectronHistos[0], ProtonHistos[0]

    # Merge the histograms
    TotalElectronHistos = mergeHistograms(ElectronHistos)
    TotalProtonHistos = mergeHistograms(ProtonHistos)

    return TotalElectronHistos, TotalProtonHistos


if __name__ == "__main__":
    Path = "/l/triton_work/Fluence_Histograms/CarringtonShielded/Carrington-SEP-Expected-Int/1mm/Res/"
    ElectronFluence, ProtonFluence = totalFluenceHistos(Path)

    # Construct the path for the 'Plot' directory
    PlotDir = os.path.join(os.path.dirname(Path), '../')
    
    # Set the grid to be behind the plot
    plt.rc('axes', axisbelow=True)

    if ElectronFluence['value'].sum() != 0:

        plt.figure(0)
        plt.bar(ElectronFluence['lower'], ElectronFluence['value'], width=ElectronFluence['upper'] - ElectronFluence['lower'], yerr=ElectronFluence['error'], align='edge')
        plt.xlabel('Energy [MeV]')
        plt.ylabel('Fluence [counts/cm2]')
        plt.title('Electron Fluence Spectrum')
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(which='both')

        print("Saving Electron Histogram: ", os.path.join(PlotDir, 'Electron Fluence Spectrum.pdf'))
        plt.savefig(os.path.join(PlotDir, 'Electron Fluence Spectrum.pdf'), format='pdf', bbox_inches='tight')

    if ProtonFluence['value'].sum() != 0:

        plt.figure(1)
        plt.bar(ProtonFluence['lower'], ProtonFluence['value'], width=ProtonFluence['upper'] - ProtonFluence['lower'], yerr=ProtonFluence['error'], align='edge')
        plt.xlabel('Energy [MeV]')
        plt.ylabel('Fluence [counts/cm2]')
        plt.title('Proton Fluence Spectrum')
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(which='both')

        print("Saving Proton Histogram: ", os.path.join(PlotDir, 'Proton Fluence Spectrum.pdf'))
        plt.savefig(os.path.join(PlotDir, 'Proton Fluence Spectrum.pdf'), format='pdf', bbox_inches='tight')

    # plt.show()

