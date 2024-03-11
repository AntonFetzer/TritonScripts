import os
import numpy as np
from GRAS.Read.ReadGRASdose import readGRASdose
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
        RawData.append(readGRASdose(path + File))

    # RawData[File][Variable][Tile]

    # Ionising dose in rad
    Dose = np.zeros(np.shape(RawData[0][0]), dtype=np.float64)

    # Number of entries. This number is the same for each tile, because each detected particle causes an entry at all tiles. Most of these entries are therefore empty.
    Entries = np.zeros(np.shape(RawData[0][2]), dtype=np.float64)

    # Number of non-zero entries for each tile. This number is different for each tile.
    NonZeroEntries = np.zeros(np.shape(RawData[0][3]), dtype=np.float64)

    for Temp in RawData:
        #print(Temp[0] * Temp[2])
        Dose += Temp[0] * Temp[2]  # Weighted average of all the dose results. Multiply dose by number of entries to
        Entries += Temp[2]  # Summed up number of entries from all files
        NonZeroEntries += Temp[3]  # Summed up number of Non Zero Entries from all files

    Dose = Dose / Entries  # Divide by summed up number of entries to get the weighted average of the dose results


    # Error estimation based on the standard deviaiton of the poulation and distance from the mean.
    ErrTerm = np.zeros(np.shape(RawData[0][1]), dtype=np.float64)

    for Temp in RawData:
        ErrTerm += Temp[2] ** 2 * (Temp[1] ** 2 + (Temp[0] - Dose) ** 2)
        #ErrTerm += Temp[2] ** 2 * (Temp[1] ** 2)

    Error = np.sqrt(ErrTerm) / Entries


    # plt.plot(100 * Error / Dose, '.')
    # plt.show()

    # Entries is the same for all tiles, because each tile gets an entry if a particle is detected in any tile.
    TotalEntries = Entries[0]

    NumTiles = len(Dose)

    MinZE = min(NonZeroEntries)
    LowestTile = np.argmin(NonZeroEntries)
    PartPerTile = Entries / NumTiles


    print("Total number of particles:", f"{TotalEntries:.2}")
    print("Minimum Particles per tile:", f"{min(PartPerTile):.2}")
    print("Minimum number of non Zero Entries:", f"{MinZE:.2}", "at Tile number", str(LowestTile))

    if MinZE < 100:
        print("ERROR !!! Tile", LowestTile, "has only", MinZE, "Non-Zero entries !!!")

    RelativeError = Error / Dose
    MaxRelativeError = max(RelativeError)
    # print(RelativeError)
    MaxRelativeErrorTile = np.argmax(RelativeError)
    print("Maximum relative error =", str(round(MaxRelativeError * 100, 2)), "% at Tile number", str(MaxRelativeErrorTile))

    PartNumRequired = (MaxRelativeError/0.01)**2 * TotalEntries  # Number of Particles required to get less than 1% relative error
    FullRunns = PartNumRequired/2e9

    print("The number of particles required to achieve 1% err is:", f"{PartNumRequired:.2}", "or", round((MaxRelativeError/0.01)**2, 3), " times the number of particles in the run")

    if MaxRelativeError > 1:
        print("ERROR !!! Tile " + str(MaxRelativeErrorTile) + " has " + str(MaxRelativeError * 100) + " % relative error!!!")

    Data = np.asarray([Dose, Error, Entries, NonZeroEntries])

    return Data


if __name__ == "__main__":
    Path = "/l/triton_work/ShieldingCurves/Carrington/CarringtonElectronDiffPowTabelated-10mm/Res/"

    TID = totalkRadGras(Path, "")

    # Relative Error in %
    plt.figure(1)
    plt.plot(100 * TID[1] / TID[0], '.')  # Relative Error in %
    plt.grid(which="both")
    plt.title("Relative Error in %")

    plt.savefig("/l/triton_work/ShieldingCurves/Carrington/CarringtonElectronDiffPowTabelated-10mm/Plot/RelativeError.pdf", format='pdf', bbox_inches="tight")

    # Number of non-zero entries
    plt.figure(2)
    plt.plot(TID[3], '.')
    plt.grid(which="both")
    plt.yscale("log")
    plt.title("Number of non-Zero Entries")

    plt.savefig("/l/triton_work/ShieldingCurves/Carrington/CarringtonElectronDiffPowTabelated-10mm/Plot/NonZeroEntries.pdf", format='pdf', bbox_inches="tight")

    # Total Dose
    x = np.linspace(0, len(TID[0])-1, num=len(TID[0]), dtype=int, endpoint=True)
    plt.figure(3)
    plt.errorbar(x, TID[0], yerr=TID[1], fmt=' ', capsize=5)
    plt.grid(which="both")
    plt.yscale("log")
    plt.title("Dose")

    plt.savefig("/l/triton_work/ShieldingCurves/Carrington/CarringtonElectronDiffPowTabelated-10mm/Plot/Dose.pdf", format='pdf', bbox_inches="tight")

    #plt.show()



