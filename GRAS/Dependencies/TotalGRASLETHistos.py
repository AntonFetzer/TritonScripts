import numpy as np
import matplotlib.pyplot as plt
from GRAS.Read.ReadGRASLETHistos import readGRASLETHistos
import sys
import os

def totalGRASLETHistos(path, particle: str):

    print("")
    print("Reading in all", particle, "files in folder:", path)

    # Get list of all csv files in Path
    Files = [f for f in os.listdir(path) if particle in f]

    NumFiles = len(Files)

    print("Number of Files:", NumFiles)

    if not Files:
        sys.exit("ERROR !!! No files found")

    RawData = []

    for File in Files:
        RawData.append(readGRASLETHistos(path + File))

    #print("RawDataShape", np.shape(RawData))
    # ( File# , Dose or Primary , Bin#, Var# )

    Data = np.array(RawData[0])
    #print(np.shape(Data))
                    # LET            Eff
    lowerID = 0     # MeV/cm          MeV/cm
    upperID = 1     # MeV/cm          MeV/cm
    meanID = 2      # MeV/cm          MeV/cm
    valueID = 3     # counts          counts
    errorID = 4     # counts          counts
    entriesID = 5   # Num             Num

    MeanCountList = np.zeros(len(Data[0, :, 0]))

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

                if Bin[meanID] != 0 and h==0: # This messes up the Eff histograms, but I don't care
                    MeanCountList[b] += 1
                    Data[h][b][meanID] += Bin[meanID]

                Data[h][b][valueID] += Bin[valueID]
                Data[h][b][errorID] += Bin[errorID]**2
                Data[h][b][entriesID] += Bin[entriesID]


    for b, MeanCount in enumerate(MeanCountList):
        if MeanCount != 0:
            Data[:, b, meanID] = Data[:, b, meanID] / MeanCount
    Data[:, :, valueID] = Data[:, :, valueID] / NumFiles
    Data[:, :, errorID] = np.sqrt(Data[:, :, errorID]) / NumFiles

    return Data


if __name__ == "__main__":

    # Only works if all input files have the same number of particle!!!!!
    path = "/home/anton/Desktop/triton_work/LET/LETMono10MeV/1mm/Res/"

    LETHist, EffHist = totalGRASLETHistos(path, "Prot")

    print("LETHist Shape", np.shape(LETHist))

    lowerID = 0
    upperID = 1
    meanID = 2
    valueID = 3
    errorID = 4
    entriesID = 5

    ### LET by Entries ###############
    NumberEntriesLETHist = sum(LETHist[:, entriesID])
    TotalLETbyEntries = sum(LETHist[:, meanID] * LETHist[:, entriesID])

    plt.figure(0)
    plt.bar(LETHist[:, lowerID], LETHist[:, entriesID], width=LETHist[:, upperID] - LETHist[:, lowerID], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("LET Histogram " + f"{NumberEntriesLETHist:.2}" + " entries\nTotal LET by Entries " + f"{TotalLETbyEntries:.2}" + " MeV/cm")
    plt.xlabel("LET [MeV/cm]")
    plt.ylabel("Number of entries per LET bin")

    ### LET by Values ###############
    TotalLETbyValues = sum(LETHist[:, meanID] * LETHist[:, valueID])

    plt.figure(1)
    plt.bar(LETHist[:, lowerID], LETHist[:, valueID], width=LETHist[:, upperID] - LETHist[:, lowerID], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("LET Histogram " + f"{NumberEntriesLETHist:.2}" + " entries\nTotal LET by Values " + f"{TotalLETbyValues:.2}" + " MeV/cm")
    plt.xlabel("LET [MeV/cm]")
    plt.ylabel("Rate per LET bin [cm-2 s-1]")

    ### Eff by Entries ###############
    NumberEntriesEffHist = sum(EffHist[:, entriesID])
    TotalEffbyEntries = sum(EffHist[:, meanID] * EffHist[:, entriesID])

    plt.figure(2)
    plt.bar(EffHist[:, lowerID], EffHist[:, entriesID], width=EffHist[:, upperID] - EffHist[:, lowerID], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("EffLET Histogram " + f"{NumberEntriesEffHist:.2}" + " entries\nTotal EffLET by Entries " + f"{TotalLETbyEntries:.2}" + " MeV/cm")
    plt.xlabel("EffLET [MeV/cm]")
    plt.ylabel("Number of entries per EffLET bin")

    TotalEffLETbyValues = sum(EffHist[:, meanID] * EffHist[:, valueID])

    plt.figure(3)
    plt.bar(EffHist[:, lowerID], EffHist[:, valueID], width=EffHist[:, upperID] - EffHist[:, lowerID], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("EffLET Histogram " + f"{NumberEntriesEffHist:.2}" + " entries\nTotal EffLET by Values " + f"{TotalEffLETbyValues:.2}" + " MeV/cm")
    plt.xlabel("EffLET [MeV/cm]")
    plt.ylabel("Rate per LET bin [cm-2 s-1]")

    plt.show()
    #plt.savefig(path1 + "../../Comp/EffHistogramComparison.pdf", format='pdf', bbox_inches="tight")
