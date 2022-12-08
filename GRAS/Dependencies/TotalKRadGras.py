import os
import numpy as np
from GRAS.Read.ReadGRASCSV import readGrasCsv
import matplotlib.pyplot as plt
import sys


def totalkRadGras(path, particle: str):
    print("")
    print("Reading in all", particle, "files in folder:", path)

    # Get list of all csv files in Path
    Files = [f for f in os.listdir(path) if particle in f]

    if not Files:
        sys.exit("ERROR !!! No files found")

    RawData = []

    for File in Files:
        RawData.append(readGrasCsv(path + File))

    # RawData[File][Variable][Tile]

    Dose = np.zeros(np.shape(RawData[0][0]), dtype=np.float64)
    Entries = np.zeros(np.shape(RawData[0][2]), dtype=np.float64)
    NonZeroEntries = np.zeros(np.shape(RawData[0][3]), dtype=np.float64)

    for Temp in RawData:
        #print(Temp[0] * Temp[2])
        Dose += Temp[0] * Temp[2]  # Weighted average of all the dose results
        Entries += Temp[2]  # Total number of Entries
        NonZeroEntries += Temp[3]  # Total number of Non Zero Entries

    Dose = Dose / Entries  # Divide by total number of entries to get the weighted average of the dose results

    ErrTerm = np.zeros(np.shape(RawData[0][1]), dtype=np.float64)

    for Temp in RawData:
        ErrTerm += Temp[2] ** 2 * (Temp[1] ** 2 + (Temp[0] - Dose) ** 2)
        #ErrTerm += Temp[2] ** 2 * (Temp[1] ** 2)

    Error = np.sqrt(ErrTerm) / Entries

    # plt.plot(100 * Error / Dose, '.')
    # plt.show()

    NumTiles = len(Dose)

    MinZE = min(NonZeroEntries)
    LowestTile = np.argmin(NonZeroEntries)
    PartPerTile = Entries[0] / NumTiles

    print("Total number of particles:", f"{Entries[0]:.2}")
    print("Particles per tile:", f"{PartPerTile:.2}")
    print("Minimum number of non Zero Entries:", f"{MinZE:.2}", "at Tile number", str(LowestTile))

    if MinZE < 100:
        print("ERROR !!! Tile", LowestTile, "has only", MinZE, "Non-Zero entries !!!")

    RelativeError = Error / Dose
    MaxRelativeError = max(RelativeError)
    # print(RelativeError)
    MaxRelativeErrorTile = np.argmax(RelativeError)
    print("Maximum relative error =", str(round(MaxRelativeError * 100, 2)), "% at Tile number", str(MaxRelativeErrorTile))

    PartNumRequired = (MaxRelativeError/0.01)**2 * Entries[0]
    FullRunns = PartNumRequired/2e9

    print("The number of particles required to achieve 1% err is:", f"{PartNumRequired:.2}", "or", int(FullRunns), "2e9 Runns")

    if MaxRelativeError > 1:
        print("ERROR !!! Tile " + str(MaxRelativeErrorTile) + " has " + str(MaxRelativeError * 100) + " % relative error!!!")

    Data = np.asarray([Dose, Error, Entries, NonZeroEntries])

    return Data


if __name__ == "__main__":
    Path = "/home/anton/Desktop/triton_work/3MatTriangles/Al-FlatTest/Res/"

    AorB = "B"

    #Electrons = totalkRadGras(Path, "Elec")
    #Electrons = totalkRadGras(Path, "Electrons2Mat" + AorB)
    Protons = totalkRadGras(Path, "Prot")
    #Protons = totalkRadGras(Path, "Protons2Mat" + AorB)

    #Total = Electrons + Protons
    #Total[1] = np.sqrt(Electrons[1] * Electrons[1] + Protons[1] * Protons[1])

#    for i, x in enumerate(Protons[0]):
#        print(i, Electrons[0][i], Electrons[1][i], Protons[0][i], Protons[1][i], Total[0][i], Total[1][i])
        # print(i, Electrons[0][i], 100 * Electrons[1][i] / Electrons[0][i], Protons[0][i], 100 * Protons[1][i] / Protons[0][i])

    # x = np.linspace(0, 99, num=100, dtype=np.float64)

    plt.figure(1)
    #plt.plot(100 * Electrons[1] / Electrons[0], '.', label="Electrons")
    plt.plot(100 * Protons[1] / Protons[0], '.', label="Protons")
    plt.legend()
    plt.title("Relative Error in %")

    plt.figure(2)
    #plt.plot(Electrons[3], '.', label="Electrons")
    plt.plot(Protons[3], '.', label="Protons")
    plt.legend()
    plt.yscale("log")
    plt.title("Number of non-Zero Entries")

    plt.figure(3)
    #plt.plot(Electrons[0], '.', label="Electrons")
    plt.plot(Protons[0], '.', label="Protons")
    #plt.plot(Protons[0] + Electrons[0], '.', label="Total")
    plt.legend()
    plt.grid()
    plt.yscale("log")
    plt.title("Dose")

    plt.show()


#print(np.geomspace(0.2, 3.2, num=9, endpoint=True))
