import numpy as np
import matplotlib.pyplot as plt
from Read.ReadSourceHistos import readSourceHistos
from Dependencies.MergeHistograms import mergeHistograms
import sys
import os

def totalSourceHistos(path):
    """
    Accumulates and processes GRAS Source Histograms from multiple paralel GRAS Result files.

    Args:
        path (str): The path to the Res folder containing the csv files.

    Returns:
        dict: A dictionary containing the accumulated Source histos.

    Raises:
        SystemExit: If no files are found in the specified path.
        SystemExit: If the histogram bins in a file do not match the bins in the first file.
    """
    print("\nReading in all csv files in folder:", path)

    # Get list of all csv files in Path
    Files = [f for f in os.listdir(path) if "Source" in f and f.endswith(".csv")]

    NumFiles = len(Files)
    print("Number of Files:", NumFiles)

    if not Files:
        sys.exit("ERROR !!! No files found")


    # List of known histograms in the expected order
    histo_names = ['Energy', 'Momentum', 'Phi', 'Theta', 'Weights']

    # Initialise the total histograms dict
    # Each histogram type initialised with an empty list, to store the actual histograms
    TotalSourceHistos = {name: [] for name in histo_names}

    for File in Files:
        SourceHistos = readSourceHistos(os.path.join(path, File))
        
        for name in histo_names:
            TotalSourceHistos[name].append(SourceHistos[name])

    # Merge the histograms
    for name in histo_names:
        TotalSourceHistos[name] = mergeHistograms(TotalSourceHistos[name])

    return TotalSourceHistos


if __name__ == "__main__":
    Path = "/l/triton_work/SourceHistograms/Carrington/CarringtonElectronDiffPow/Res/"
    TotalSourceHists = totalSourceHistos(Path)

    # Construct the path for the 'Plot' directory
    PlotDir = os.path.join(os.path.dirname(Path), '../Plot')
    
    # Set the grid to be behind the plot
    plt.rc('axes', axisbelow=True)

    # Plot the total energy histogram
    TotalEnergyHist = TotalSourceHists['Energy']
    plt.figure(0)
    plt.bar(TotalEnergyHist['lower'], TotalEnergyHist['value'], width=TotalEnergyHist['upper'] - TotalEnergyHist['lower'], yerr=TotalEnergyHist['error'], align='edge')
    plt.xlabel('Energy [MeV]')
    plt.ylabel('Counts per bin')
    plt.title('Total Energy Histogram')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(which='both')

    print("Saving Energy Histogram: ", os.path.join(PlotDir, 'Energy Source Spectrum.pdf'))
    plt.savefig(os.path.join(PlotDir, 'Energy Source Spectrum.pdf'), format='pdf', bbox_inches='tight')
    
    # Plot the total momentum histogram
    TotalMomentumHist = TotalSourceHists['Momentum']
    plt.figure(1)
    plt.bar(TotalMomentumHist['lower'], TotalMomentumHist['value'], width=TotalMomentumHist['upper'] - TotalMomentumHist['lower'], yerr=TotalMomentumHist['error'], align='edge')
    plt.xlabel('Momentum [MeV]')
    plt.ylabel('Counts per bin')
    plt.title('Total Momentum Histogram')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(which='both')

    print("Saving Momentum Histogram: ", os.path.join(PlotDir, 'Momentum Source Spectrum.pdf'))
    plt.savefig(os.path.join(PlotDir, 'Momentum Source Spectrum.pdf'), format='pdf', bbox_inches='tight')
    
    # Plot the total phi histogram
    TotalPhiHist = TotalSourceHists['Phi']
    plt.figure(2)
    plt.bar(TotalPhiHist['lower'], TotalPhiHist['value'], width=TotalPhiHist['upper'] - TotalPhiHist['lower'], yerr=TotalPhiHist['error'], align='edge')
    plt.xlabel('Phi [degrees]')
    plt.ylabel('Counts per bin')
    plt.title('Total Phi Histogram')
    plt.grid(which='both')

    print("Saving Phi Histogram: ", os.path.join(PlotDir, 'Phi Source Distribution.pdf'))
    plt.savefig(os.path.join(PlotDir, 'Phi Source Distribution.pdf'), format='pdf', bbox_inches='tight')

    # Plot the total theta histogram
    TotalThetaHist = TotalSourceHists['Theta']
    plt.figure(3)
    plt.bar(TotalThetaHist['lower'], TotalThetaHist['value'], width=TotalThetaHist['upper'] - TotalThetaHist['lower'], yerr=TotalThetaHist['error'], align='edge')
    plt.xlabel('Theta [degrees]')
    plt.ylabel('Counts per bin')
    plt.title('Total Theta Histogram')
    plt.grid(which='both')

    print("Saving Theta Histogram: ", os.path.join(PlotDir, 'Theta Source Distribution.pdf'))
    plt.savefig(os.path.join(PlotDir, 'Theta Source Distribution.pdf'), format='pdf', bbox_inches='tight')

    # Plot the total weights histogram
    TotalWeightsHist = TotalSourceHists['Weights']
    plt.figure(4)
    plt.bar(TotalWeightsHist['lower'], TotalWeightsHist['value'], width=TotalWeightsHist['upper'] - TotalWeightsHist['lower'], yerr=TotalWeightsHist['error'], align='edge')
    plt.xlabel('Weights')
    plt.ylabel('Counts')
    plt.title('Total Weights Histogram')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(which='both')

    print("Saving Weights Histogram: ", os.path.join(PlotDir, 'Weights Source Distribution.pdf'))
    plt.savefig(os.path.join(PlotDir, 'Weights Source Distribution.pdf'), format='pdf', bbox_inches='tight')

    # plt.show()