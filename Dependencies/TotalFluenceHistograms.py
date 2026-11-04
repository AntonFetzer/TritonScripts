import matplotlib.pyplot as plt
from Read.ReadFluenceHistos import readFluenceHistos
from Dependencies.MergeHistograms import mergeHistograms, completeResultFiles
import os

def totalFluenceHistos(path):
    """
    Accumulates and processes GRAS Fluence Histograms from multiple paralel GRAS Result files.

    Args:
        path (str): The path to the Res folder containing the csv files.

    Returns:
        tuple: The accumulated electron and proton histograms, or (None, None) if the
            folder has no complete csv file.

    Raises:
        SystemExit: If the histogram bins in a file do not match the bins in the first file.
    """
    print("\nReading in all csv files in folder:", path)

    # Get list of all complete csv files in Path
    Files = completeResultFiles(path, "Fluence")

    if not Files:
        print(f"ERROR !!! No complete csv files in {path}")
        return None, None

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

