import numpy as np
import matplotlib.pyplot as plt
from GRAS.Read.ReadGRASHistos import readGRASHistos
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

    RawData = []

    for File in Files:
        RawData.append(readGRASHistos(path + File))

    #print("RawDataShape", np.shape(RawData))
    # ( File# , Dose or Primary , Bin#, Var# )

    Data = RawData[0]

                    #   Dose            Primary
    lowerID = 0     #   rad/s           MeV
    upperID = 1     #   rad/s           MeV
    meanID = 2      #   rad/s           MeV
    valueID = 3     #   counts/cm2      rad/s
    errorID = 4     #   counts/cm2      rad/s
    entriesID = 5   #   Num             Num

    for f, file in enumerate(RawData):
        if f == 0:
            continue
        for h, hist in enumerate(file):
            for b, Bin in enumerate(hist):
                if Data[h][b][lowerID] != Bin[lowerID]:
                    print("CANNOT ADD HISTOGRAMS Lower missmatch")
                    print("Data", Data[h][b][lowerID], "Bin", Bin[lowerID])
                    Data = 0
                    break
                if Data[h][b][upperID] != Bin[upperID]:
                    print("CANNOT ADD HISTOGRAMS Upper missmatch")
                    Data = 0
                    break

                Data[h][b][valueID] += Bin[valueID]
                Data[h][b][errorID] += Bin[errorID]
                Data[h][b][entriesID] += Bin[entriesID]

    DoseHist, PrimaryHist = Data

    PrimaryHist[:, 3] = PrimaryHist[:, 3] / NumFiles
    PrimaryHist[:, 4] = PrimaryHist[:, 4] / NumFiles

    DoseHist[:, 3] = DoseHist[:, 3] / NumFiles
    DoseHist[:, 4] = DoseHist[:, 4] / NumFiles

    return Data



if __name__ == "__main__":

    # Only works if all input files have the same number of particles!!!!!
    path = "/home/anton/Desktop/triton_work/CARRINGTON/HistogramSEP/Res/"
    DoseHist, PrimaryHist = totalGRASHistos(path, "Prot")

    print("DoseHist Shape", np.shape(DoseHist))

    lowerID = 0
    upperID = 1
    meanID = 2
    valueID = 3
    errorID = 4
    entriesID = 5

    ####    Dose hist Entries #######
    NumberEntries = sum(DoseHist[:, entriesID])
    DoseEntries = sum(DoseHist[:, meanID] * DoseHist[:, entriesID])

    plt.figure(0)
    plt.bar(DoseHist[:, lowerID], DoseHist[:, entriesID], width=DoseHist[:, upperID] - DoseHist[:, lowerID], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("Dose Entries Histogram\n" + f"{NumberEntries:.2}" + " entries " + f"{DoseEntries:.2}" + " krad total dose ?!?")
    plt.xlabel("Dose [krad per Month]")
    plt.ylabel("Number of entries per dose bin")

    ####    Dose hist Values #######
    SumValues = sum(DoseHist[:, valueID])
    DoseValues = sum(DoseHist[:, meanID] * DoseHist[:, valueID])

    plt.figure(1)
    plt.bar(DoseHist[:, lowerID], DoseHist[:, valueID], width=DoseHist[:, upperID] - DoseHist[:, lowerID], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("Dose Values Histogram\n" + f"{SumValues:.2}" + " SumValues " + f"{DoseValues:.2}" + " krad total dose ?!?")
    plt.xlabel("Dose [krad per Month] ?")
    plt.ylabel("Number of entries per dose bin")

    ####    Primary hist Entries #######
    NumberEntries = sum(PrimaryHist[:, entriesID])

    plt.figure(2)
    plt.bar(PrimaryHist[:, lowerID], PrimaryHist[:, entriesID], width=PrimaryHist[:, upperID] - PrimaryHist[:, lowerID], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("Entries VS primary kinetic energy\n" + f"{NumberEntries:.2}" + " total Entries")
    plt.xlabel("Kinetic energy [MeV]")
    plt.ylabel("Number of entries")
    plt.savefig(path + "../PrimaryHistEntries.eps", format='eps', bbox_inches="tight")

    ####    Primary hist Values #######
    TotalDose = sum(PrimaryHist[:, valueID])

    plt.figure(3)
    plt.bar(PrimaryHist[:, lowerID], PrimaryHist[:, valueID], width=PrimaryHist[:, upperID] - PrimaryHist[:, lowerID], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("Dose deposited VS primary kinetic energy\n" + f"{NumberEntries:.2}" + " total Entries " + f"{TotalDose:.2}" + " krad total dose")
    plt.xlabel("Kinetic energy [MeV]")
    plt.ylabel("Dose [krad per Month]")

    plt.show()
