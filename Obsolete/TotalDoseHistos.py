import numpy as np
import matplotlib.pyplot as plt
from GRAS.Read.ReadDoseHistos import readDoseHistos
import sys
import os

def totalGRASHistos(path, particle: str):
    print("")
    print("Reading in all", particle, "files in folder:", path)

    # Get list of all csv files in Path
    Files = [f for f in os.listdir(path) if particle in f]

    NumFiles = len(Files)
    print("Number of files:", NumFiles)

    if not Files:
        sys.exit("ERROR !!! No files found")

    # Initialize Data dictionary with the first file
    DoseData, PrimaryData = readDoseHistos(path + Files[0])

    TotalDoseEntries = DoseData['entries'].copy()
    TotalPrimaryEntries = PrimaryData['entries'].copy()

    # Start from the second file
    for File in Files[1:]:
        DoseFileData, PrimaryFileData = readDoseHistos(path + File)

        # Check that the bins for DoseHist align
        if not np.allclose(DoseData['lower'], DoseFileData['lower']) or \
           not np.allclose(DoseData['upper'], DoseFileData['upper']):
            print("CANNOT ADD HISTOGRAMS Lower or Upper missmatch for DoseHist")
            return None, None

        # Check that the bins for PrimaryHist align
        if not np.allclose(PrimaryData['lower'], PrimaryFileData['lower']) or \
           not np.allclose(PrimaryData['upper'], PrimaryFileData['upper']):
            print("CANNOT ADD HISTOGRAMS Lower or Upper missmatch for PrimaryHist")
            return None, None

        # Update mean with weighted average
        DoseData['mean'] = (DoseData['mean'] * DoseData['entries'] + DoseFileData['mean'] * DoseFileData['entries']) / \
                           (DoseData['entries'] + DoseFileData['entries'])
        PrimaryData['mean'] = (PrimaryData['mean'] * PrimaryData['entries'] + PrimaryFileData['mean'] * PrimaryFileData['entries']) / \
                              (PrimaryData['entries'] + PrimaryFileData['entries'])

        # Update error using error propagation rule
        DoseData['error'] = np.sqrt(DoseData['error']**2 + DoseFileData['error']**2)
        PrimaryData['error'] = np.sqrt(PrimaryData['error']**2 + PrimaryFileData['error']**2)

        # Weighted addition for 'value' data
        DoseData['value'] = DoseData['value'] * DoseData['entries'] + DoseFileData['value'] * DoseFileData['entries']
        PrimaryData['value'] = PrimaryData['value'] * PrimaryData['entries'] + PrimaryFileData['value'] * PrimaryFileData['entries']

        # Add the entries for DoseHist and PrimaryHist
        DoseData['entries'] += DoseFileData['entries']
        PrimaryData['entries'] += PrimaryFileData['entries']

        # Update total entries
        TotalDoseEntries += DoseFileData['entries']
        TotalPrimaryEntries += PrimaryFileData['entries']

    # At the end, perform division for 'value' and 'error' using total entries
    DoseData['value'] /= TotalDoseEntries
    PrimaryData['value'] /= TotalPrimaryEntries
    DoseData['error'] /= np.sqrt(TotalDoseEntries)
    PrimaryData['error'] /= np.sqrt(TotalPrimaryEntries)

    return DoseData, PrimaryData





if __name__ == "__main__":

    # Only works if all input files have the same number of particles!!!!!
    path = "/l/triton_work/Histograms/AE9500keV/Res/"
    DoseHist, PrimaryHist = totalGRASHistos(path, "")

    ####    Dose hist Entries #######
    NumberEntries = np.sum(DoseHist['entries'])
    DoseEntries = np.sum(DoseHist['mean'] * DoseHist['entries'])

    plt.figure(0)
    plt.bar(DoseHist['lower'], DoseHist['entries'], width=DoseHist['upper'] - DoseHist['lower'], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title(
        "Dose Entries Histogram\n" + f"{NumberEntries:.2}" + " entries " + f"{DoseEntries:.2}" + " krad total dose ?!?")
    plt.xlabel("Dose [krad per Month]")
    plt.ylabel("Number of entries per dose bin")

    ####    Dose hist Values #######
    SumValues = np.sum(DoseHist['value'])
    DoseValues = np.sum(DoseHist['mean'] * DoseHist['value'])

    plt.figure(1)
    plt.bar(DoseHist['lower'], DoseHist['value'], width=DoseHist['upper'] - DoseHist['lower'], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title(
        "Dose Values Histogram\n" + f"{SumValues:.2}" + " SumValues " + f"{DoseValues:.2}" + " krad total dose ?!?")
    plt.xlabel("Dose [krad per Month] ?")
    plt.ylabel("Number of entries per dose bin")

    ####    Primary hist Entries #######
    NumberEntries = np.sum(PrimaryHist['entries'])

    plt.figure(2)
    plt.bar(PrimaryHist['lower'], PrimaryHist['entries'], width=PrimaryHist['upper'] - PrimaryHist['lower'],
            align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("Entries VS primary kinetic energy\n" + f"{NumberEntries:.2}" + " total Entries")
    plt.xlabel("Kinetic energy [MeV]")
    plt.ylabel("Number of entries")

    ####    Primary hist Values #######
    TotalDose = np.sum(PrimaryHist['value'])

    plt.figure(3)
    plt.bar(PrimaryHist['lower'], PrimaryHist['value'], width=PrimaryHist['upper'] - PrimaryHist['lower'], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title(
        "Dose deposited VS primary kinetic energy\n" + f"{NumberEntries:.2}" + " total Entries " + f"{TotalDose:.2}" + " krad total dose")
    plt.xlabel("Kinetic energy [MeV]")
    plt.ylabel("Dose [krad per Month]")

    plt.show()
