import numpy as np
import matplotlib.pyplot as plt
from GRAS.Read.ReadLETHistos import readLETHistos
from GRAS.Dependencies.MergeHistograms import mergeHistograms
import sys
import os

def totalLETHistos(path):
    """
    Accumulates and processes GRAS LEThistos from multiple files.

    Args:
        path (str): The path to the folder containing the csv files.

    Returns:
        dict: A dictionary containing the accumulated LEThistos.

    Raises:
        SystemExit: If no files are found in the specified path.
        SystemExit: If the histogram bins in a file do not match the bins in the first file.

    """
    # print("\nReading in all csv files in folder:", path)

    # Get list of all csv files in Path
    Files = [f for f in os.listdir(path) if "LET" in f and f.endswith(".csv")]

    NumFiles = len(Files)
    # print("Number of Files:", NumFiles)

    if not Files:
        sys.exit("ERROR !!! No files found")

    LETList = []
    EffList = []

    # Read the files and accumulate histogram dicts in the list of histograms
    for File in Files:
        LETHistDict, EffHistDict = readLETHistos(os.path.join(path, File))
        LETList.append(LETHistDict)
        EffList.append(EffHistDict)

    # Merge the histograms
    TotalLETHist = mergeHistograms(LETList)
    TotalEffHist = mergeHistograms(EffList)
    
    return TotalLETHist, TotalEffHist


if __name__ == "__main__":

    # Only works if all input files have the same number of particle!!!!!
    path = "/l/triton_work/LET_Histograms/Carrington/ISS-AP9-mission/1mm/Res/"

    LET, Eff = totalLETHistos(path)

    ### LET by Entries ###############
    # Calculate sums and total LET by entries for LETHist
    NumberEntriesLETHist = np.sum(LET['entries'])
    TotalLETbyEntries = np.sum(LET['mean'] * LET['entries'])

    # Plotting LET Histogram by entries
    plt.figure(0)
    plt.bar(LET['lower'], LET['entries'], width=LET['upper'] - LET['lower'], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("LET Histogram " + f"{NumberEntriesLETHist:.2e}" + " entries\nTotal LET by Entries " + f"{TotalLETbyEntries:.2e}" + " [MeV cm2 mg-1]")
    plt.xlabel("LET [MeV cm2 mg-1]")
    plt.ylabel("Number of entries per LET bin")

    plt.savefig(path + "../LET-Entries.pdf", format='pdf', bbox_inches="tight")

    # Calculate total LET by values for LETHist
    TotalLETbyValues = np.sum(LET['mean'] * LET['value'])

    # Plotting LET Histogram by values and error bars
    plt.figure(1)
    plt.bar(LET['lower'], LET['value'], width=LET['upper'] - LET['lower'], align='edge', alpha=0.3)
    plt.errorbar(LET['mean'], LET['value'], LET['error'], fmt=' ', capsize=5, elinewidth=1, capthick=1, label="LET Histogram")
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("LET Histogram " + f"{NumberEntriesLETHist:.2e}" + " entries\nTotal LET by Values " + f"{TotalLETbyValues:.2e}" + " [MeV cm2 mg-1]")
    plt.xlabel("LET [MeV cm2 mg-1]")
    plt.ylabel("Rate per LET bin [cm-2 s-1]")

    plt.savefig(path + "../LET-Values.pdf", format='pdf', bbox_inches="tight")

    # plt.show()
