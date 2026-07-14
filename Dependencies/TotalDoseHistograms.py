import numpy as np
import matplotlib.pyplot as plt
from Read.ReadDoseHistos import readDoseHistos
from Dependencies.MergeHistograms import mergeHistograms
import os


def totalDoseHistograms(path, particle=""):
    """
    Merges the GRAS dose and primary-energy histograms from the csv files of
    identical runs (jobs differing only by seed), weighted by each job's
    number of primary events.

    Replaces the retired TotalDoseHistos.totalGRASHistos, which summed the
    per-file standard errors in quadrature and then divided by
    sqrt(total entries) again (shrinking the error by roughly sqrt(N/K)) and
    mis-accumulated 'value' when merging three or more files.

    Args:
        path (str): The path to the folder containing the csv files.
        particle (str): Optional filename filter, e.g. "Electrons".

    Returns:
        tuple: (DoseHist, PrimaryHist) merged histogram dicts, or (None, None)
            if no usable files are found.
    """
    print("\nReading in all", particle or "csv", "files in folder:", path)

    Files = [f for f in os.listdir(path) if particle in f and f.endswith(".csv")]

    if not Files:
        print(f"ERROR !!! No csv files found in {path}")
        return None, None

    # A job killed while writing its checkpoint leaves an empty or truncated
    # csv. Complete GRAS csv files always end with an 'End of File' line, so
    # anything else is corrupt and gets deleted before reading.
    Complete = []
    for File in Files:
        FilePath = os.path.join(path, File)
        with open(FilePath, 'rb') as f:
            Tail = f.read()[-64:]
        if b"'End of File'" in Tail:
            Complete.append(File)
        else:
            print("WARNING: deleting incomplete csv", FilePath)
            os.remove(FilePath)
    Files = Complete

    if not Files:
        print(f"ERROR !!! No complete csv files in {path}")
        return None, None

    DoseList = []
    PrimaryList = []

    for File in Files:
        DoseHistDict, PrimaryHistDict = readDoseHistos(os.path.join(path, File))
        if len(DoseHistDict['lower']) == 0:
            continue
        DoseList.append(DoseHistDict)
        PrimaryList.append(PrimaryHistDict)

    if not DoseList:
        print(f"ERROR !!! None of the {len(Files)} csv files in {path} contains histogram data")
        return None, None

    TotalDoseHist = mergeHistograms(DoseList)
    TotalPrimaryHist = mergeHistograms(PrimaryList)

    return TotalDoseHist, TotalPrimaryHist


if __name__ == "__main__":
    path = "/l/triton_work/Histograms/AE9500keV/Res/"
    DoseHist, PrimaryHist = totalDoseHistograms(path)

    ####    Dose hist Entries #######
    NumberEntries = np.sum(DoseHist['entries'])

    plt.figure(0)
    plt.bar(DoseHist['lower'], DoseHist['entries'], width=DoseHist['upper'] - DoseHist['lower'], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("Dose Entries Histogram\n" + f"{NumberEntries:.2e}" + " entries")
    plt.xlabel("Dose [krad per Month]")
    plt.ylabel("Number of entries per dose bin")

    ####    Primary hist Values #######
    TotalDose = np.sum(PrimaryHist['value'])

    plt.figure(1)
    plt.bar(PrimaryHist['lower'], PrimaryHist['value'], width=PrimaryHist['upper'] - PrimaryHist['lower'], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("Dose deposited VS primary kinetic energy\n" + f"{TotalDose:.2e}" + " krad total dose")
    plt.xlabel("Kinetic energy [MeV]")
    plt.ylabel("Dose [krad per Month]")

    plt.show()
