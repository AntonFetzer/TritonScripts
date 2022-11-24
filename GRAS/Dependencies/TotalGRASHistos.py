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
    path1 = "/home/anton/Desktop/triton_work/CARRINGTON/HistogramTestIntegralSpectrum/Res/"
    path2 = "/home/anton/Desktop/triton_work/CARRINGTON/HistogramTestDifferentialPowSpec/Res/"
    DoseHist, PrimaryHist = totalGRASHistos(path1, "Elec")

    print("DoseHist Shape", np.shape(DoseHist))

    lowerID = 0
    upperID = 1
    meanID = 2
    valueID = 3
    errorID = 4
    entriesID = 5

    '''
    NumberEntries = sum(DoseHist[:, entriesID])
    TotalDose = sum(DoseHist[:, meanID] * DoseHist[:, entriesID])

    plt.figure(1)
    plt.bar(DoseHist[:, lowerID], DoseHist[:, entriesID], width=DoseHist[:, upperID] - DoseHist[:, lowerID], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("Dose Depositions Histogram\n" + f"{NumberEntries:.2}" + " entries " + f"{TotalDose:.2}" + " krad total dose ?!?")
    plt.xlabel("Dose [krad per Month]")
    plt.ylabel("Number of entries per dose bin")
    plt.savefig(path + "../DoseDepositionsHistogram.eps", format='eps', bbox_inches="tight")
    '''

    NumberEntries = sum(PrimaryHist[:, entriesID])
    TotalDose = sum(PrimaryHist[:, valueID])

    plt.figure(2)
    plt.bar(PrimaryHist[:, lowerID], PrimaryHist[:, valueID], width=PrimaryHist[:, upperID] - PrimaryHist[:, lowerID], align='edge', alpha=0.5)
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("Dose deposited VS primary kinetic energy\n" + f"{NumberEntries:.2}" + " simulated particles " + f"{TotalDose:.2}" + " krad total dose")
    plt.xlabel("Kinetic energy [MeV]")
    plt.ylabel("Dose [krad per Month]")

    DoseHist, PrimaryHist = totalGRASHistos(path2, "Elec")
    plt.bar(PrimaryHist[:, lowerID], PrimaryHist[:, valueID], width=PrimaryHist[:, upperID] - PrimaryHist[:, lowerID], align='edge', alpha=0.5)

    plt.savefig(path1 + "../../Comp/DoseVSPrimaryComparison.pdf", format='pdf', bbox_inches="tight")

    DoseHist, PrimaryHist = totalGRASHistos(path1, "Elec")

    NumberEntries = sum(PrimaryHist[:, entriesID])

    plt.figure(3)
    plt.bar(PrimaryHist[:, lowerID], PrimaryHist[:, entriesID], width=PrimaryHist[:, upperID] - PrimaryHist[:, lowerID], align='edge', alpha=0.5)
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("Particle count VS primary kinetic energy\n" + f"{NumberEntries:.2}" + " simulated particles")
    plt.xlabel("Kinetic energy [MeV]")
    plt.ylabel("Number of entries")

    DoseHist, PrimaryHist = totalGRASHistos(path2, "Elec")
    plt.bar(PrimaryHist[:, lowerID], PrimaryHist[:, entriesID], width=PrimaryHist[:, upperID] - PrimaryHist[:, lowerID], align='edge', alpha=0.5)

    plt.savefig(path1 + "../../Comp/PrimaryHistogramComparison.pdf", format='pdf', bbox_inches="tight")
