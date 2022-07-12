import os
from ReadGRASCSV1Tile import readGrasCsv1Tile
from MeVtokRad2DGras import MeVtokRad_2D


def totalkRadGras1Tile(path):
    print("Reading in all Electrons and Protons files in folder:", path)

    MeV = {}  # Initialise empty dict !

    # Get list of all csv files in Path
    ElecFiles = [f for f in os.listdir(path) if "Electrons" in f]
    ProtFiles = [f for f in os.listdir(path) if "Protons" in f]

    MeV["Elec"] = 0
    MeV["ElErr"] = 0

    MeV["Prot"] = 0
    MeV["PrErr"] = 0

    for File in ElecFiles:
        Res = readGrasCsv1Tile(path + File)
        MeV["Elec"] += Res["Dose"]
        MeV["ElErr"] += Res["Error"]

    MeV["Elec"] = MeV["Elec"] / len(ElecFiles)
    MeV["ElErr"] = MeV["ElErr"] / len(ElecFiles)

    for File in ProtFiles:
        Res = readGrasCsv1Tile(path + File)
        MeV["Prot"] += Res["Dose"]
        MeV["PrErr"] += Res["Error"]

    MeV["Prot"] = MeV["Prot"] / len(ProtFiles)
    MeV["PrErr"] = MeV["PrErr"] / len(ProtFiles)

    kRad = {}
    NORM_FACTOR_SPECTRUM_Elec = 7.891281E+14
    NORM_FACTOR_SPECTRUM_Prot = 3.389664E+11

    kRad["Elec"] = MeVtokRad_2D(MeV["Elec"], NORM_FACTOR_SPECTRUM_Elec)
    kRad["ElErr"] = MeVtokRad_2D(MeV["ElErr"], NORM_FACTOR_SPECTRUM_Elec)

    kRad["Prot"] = MeVtokRad_2D(MeV["Prot"], NORM_FACTOR_SPECTRUM_Prot)
    kRad["PrErr"] = MeVtokRad_2D(MeV["PrErr"], NORM_FACTOR_SPECTRUM_Prot)

    return kRad


if __name__ == "__main__":
    Path = "/home/anton/Desktop/triton_work/GRAS-1Mat/ParallelTest/Results/"

    Results = totalkRadGras1Tile(Path)

    print("The total electron dose is: ", Results["Elec"], " kRad")
    print("with ", 100 * Results["ElErr"] / Results["Elec"], " % error")

    print("The total proton energy dose is: ", Results["Prot"], " kRad")
    print("with ", 100 * Results["PrErr"] / Results["Prot"], " % error")
