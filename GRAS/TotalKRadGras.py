import os
import numpy as np
from ReadGRASCSV import readGrasCsv
import matplotlib.pyplot as plt


def totalkRadGras(path, particle: str):
    print("Reading in all", particle, "files in folder:", path)

    # Get list of all csv files in Path
    Files = [f for f in os.listdir(path) if particle in f]

    if not Files:
        print("No files found")

    Dose = 0
    Error = 0
    Entries = 0
    NonZeroEntries = 0

    for i, File in enumerate(Files):
        Temp = readGrasCsv(path + File)
        Dose += Temp[0]
        Error += Temp[1]*Temp[1]
        Entries += Temp[2]
        NonZeroEntries += Temp[3]

    NumTiles = len(Dose)

    MinZE = min(NonZeroEntries)
    LowestTile = np.argmin(NonZeroEntries)
    PartPerTile = Entries[0]/NumTiles

    print("Particles per tile: ", str(PartPerTile), " Minimum number of non Zero Entries: ", str(MinZE))
    if MinZE < 10:
        print("WARNING !!! Tile ", str(LowestTile), " has only ", str(MinZE), " Non-Zero entries !!! ")
        return 0

    Data = np.asarray([Dose, np.sqrt(Error), Entries, NonZeroEntries])

    # Dose data is in rad/s --> multiply with number of seconds in a month to get to dose per months.
    # Dose is given per generated particle --> need to divide by the number of files
    # Dose is given in rad --> divide by 1000 to get to krad.
    Data[0] = Data[0] * 30 * 24 * 60 * 60 / len(Files) / 1000
    Data[1] = Data[1] * 30 * 24 * 60 * 60 / len(Files) / 1000

    return Data


if __name__ == "__main__":
    Path = "/home/anton/Desktop/triton_work/Permutations/2Layer/Res/"

    Electrons = totalkRadGras(Path, "Elec")
    Protons = totalkRadGras(Path, "Prot")
    '''
    for i, x in enumerate(Protons[0]):
        print(i, Electrons[0][i], 100 * Electrons[1][i] / Electrons[0][i], Protons[0][i], 100 * Protons[1][i] / Protons[0][i])

    x = np.linspace(0, 2.5, num=101, dtype=float)

    plt.plot(Electrons[0], ".")
    plt.plot(Protons[0], ".")

    plt.show()
    '''