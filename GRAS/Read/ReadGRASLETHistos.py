import numpy as np
import matplotlib.pyplot as plt
import csv


def readGRASLETHistos(file):

    LETHist = []
    EffHist = []

    ReadFlag = 0
    #print("Reading in File: " + file)
    with open(file, 'r') as f:
        reader = csv.reader(f)
        for line in reader:
            if ReadFlag == 0:
                if "'Bin entries'" in line:
                    ReadFlag = 1
            elif ReadFlag == 1:
                if "'End of Block'" in line:
                    ReadFlag = 2
                else:
                    LETHist.append([float(x) for x in line])
            elif ReadFlag == 2:
                if "'Bin entries'" in line:
                    ReadFlag = 3
            elif ReadFlag == 3:
                if "'End of Block'" in line:
                    ReadFlag = 4
                else:
                    EffHist.append([float(x) for x in line])

    LETHist = np.asarray(LETHist)
    EffHist = np.asarray(EffHist)

    # The LET spectra are in 'counts' per MeV/cm bin ???
    # If the input spcectrum is in [cm-2 s-1] then the counts are in [s-1] per interface area between the shield and detector
    # This area is 1e6 cm2 in the 1Tile.gdml file written in December 2022
    # To get the values in [cm-2 s-1] they have to be divided by 1e6
                        # LET            Eff
    valueID = 3     #   counts          counts
    errorID = 4     #   counts          counts

    AreaNormFactor = 1e6

    LETHist[:, valueID] = LETHist[:, valueID] / AreaNormFactor
    LETHist[:, errorID] = LETHist[:, errorID] / AreaNormFactor
    EffHist[:, valueID] = EffHist[:, valueID] / AreaNormFactor
    EffHist[:, errorID] = EffHist[:, errorID] / AreaNormFactor

    return LETHist, EffHist



if __name__ == "__main__":
    file = "/home/anton/Desktop/triton_work/LET/0mm/Res/Protons1Tile_700548_108379.csv"

    LETHist, EffHist = readGRASLETHistos(file)

    print("LETHist Shape", np.shape(LETHist))

                    # LET            Eff
    lowerID = 0     #   MeV/cm          MeV/cm
    upperID = 1     #   MeV/cm          MeV/cm
    meanID = 2      #   MeV/cm          MeV/cm
    valueID = 3     #   counts          counts
    errorID = 4     #   counts          counts
    entriesID = 5   #   Num             Num

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

    TotalLETbyValues = sum(LETHist[:, meanID] * LETHist[:, valueID])

    plt.figure(1)
    plt.bar(LETHist[:, lowerID], LETHist[:, valueID], width=LETHist[:, upperID] - LETHist[:, lowerID], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("LET Histogram " + f"{NumberEntriesLETHist:.2}" + " entries\nTotal LET by Values " + f"{TotalLETbyValues:.2}" + " MeV/cm")
    plt.xlabel("LET [MeV/cm]")
    plt.ylabel("Rate per LET bin [cm-2 s-1]")

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

    NumberEntries = sum(EffHist[:, entriesID])

    plt.show()
