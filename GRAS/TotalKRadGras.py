import os
import numpy as np
from ReadGRASCSV import readGrasCsv
from MeVtokRad2DGras import MeVtokRad_2D
import matplotlib.pyplot as plt


def totalkRadGras(path, particle: str):
    print("Reading in all", particle, "files in folder:", path)

    # Get list of all csv files in Path
    Files = [f for f in os.listdir(path) if particle in f]

    NumTiles = np.shape(readGrasCsv(path + Files[0]))[1]

    Data = np.zeros((2, NumTiles), dtype=float)

    for i, File in enumerate(Files):
        Data += readGrasCsv(path + File)

    Data = Data / len(Files)

    if "lec" in particle:
        NORM_FACTOR_SPECTRUM = 7.891281E+14
    elif "rot" in particle:
        NORM_FACTOR_SPECTRUM = 3.389664E+11
    else:
        NORM_FACTOR_SPECTRUM = "NaN"
        print("!!!! ERROR !!!! WRONG PARTICLE TYPE SPECIFIED")

    Data = MeVtokRad_2D(Data, NORM_FACTOR_SPECTRUM)

    return Data


if __name__ == "__main__":
    Path = "/home/anton/Desktop/triton_work/GRAS-2Mat/Test/Results/"

    Electrons = totalkRadGras(Path, "Elec")*99

    Protons = totalkRadGras(Path, "Prot")*99

    plt.plot(Electrons[0])
    plt.plot(Protons[0])
    plt.show()

    # print("The total electron dose is: ", Electrons[0], " kRad")
    # print("with ", 100 * Electrons[1] / Electrons[0], " % error")

    # print("The total proton energy dose is: ", Results["Prot"], " kRad")
    # print("with ", 100 * Results["PrErr"] / Results["Prot"], " % error")
