import numpy as np
import matplotlib.pyplot as plt
from Read.ReadLETHistos import readLETHistos
from Dependencies.MergeHistograms import mergeHistograms
import sys
import os



def totalLETHistos(path):
    """
    Merges GRAS LET histograms from the csv files of identical runs.

    Args:
        path (str): The path to the folder containing the csv files.

    Returns:
        dict: A dictionary containing the accumulated LET histograms.

    Raises:
        SystemExit: If no files are found in the specified path.
        SystemExit: If the histogram bins in a file do not match the bins in the first file.

    """
    print("\nReading in all csv files in folder:", path)

    # Get list of all csv files in Path
    Files = [f for f in os.listdir(path) if "LET" in f and f.endswith(".csv")]

    if not Files:
        print("ERROR !!! No files found")
        # Wait for user input before exiting
        input("Press Enter to continue")
        return None, None

    LETList = []
    EffList = []

    # Read the files and append histogram dicts to the list of histograms
    for File in Files:
        LETHistDict, EffHistDict = readLETHistos(os.path.join(path, File))
        LETList.append(LETHistDict)
        EffList.append(EffHistDict)

    # Merge the histograms
    TotalLETHist = mergeHistograms(LETList)
    TotalEffHist = mergeHistograms(EffList)
    
    return TotalLETHist, TotalEffHist


if __name__ == "__main__":

    path = "/l/triton_work/LET_Histograms/Carrington/"

    # Find all subdirectories in the given path that contain a "Res" subfolder
    # and calculate the total LET histos for each of them
    for root, dirs, files in os.walk(path):
        if 'Res' in dirs:
            folder = os.path.join(root, 'Res/')

            LET, Eff = totalLETHistos(folder)

            # Check if the LET histogram is empty
            if LET['entries'].sum() == 0:
                print("LET histogram is empty, skipping folder:", folder)
                continue

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

            plt.savefig(folder + "../LET-Entries.pdf", format='pdf', bbox_inches="tight")

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

            plt.savefig(folder + "../LET-Values.pdf", format='pdf', bbox_inches="tight")
            plt.close("all")
            #plt.show()
