import os
from ReadGRASCSV1Tile import readGrasCsv1Tile


def totalMeVGras(path):
    print("Reading in all Electrons and Protons files in folder:", path)

    Total = {}  # Initialise empty dict !

    # Get list of all csv files in Path
    ElecFiles = [f for f in os.listdir(path + "/Results/") if "Electrons" in f]
    ProtFiles = [f for f in os.listdir(path + "/Results/") if "Protons" in f]

    Total["ElecMeV"] = 0
    Total["ElecError"] = 0

    Total["ProtMeV"] = 0
    Total["ProtError"] = 0

    for File in ElecFiles:
        Res = readGrasCsv1Tile(path + "/Results/" + File)
        Total["ElecMeV"] += Res["Dose"]
        Total["ElecError"] += Res["Error"]

    Total["ElecMeV"] = Total["ElecMeV"] / len(ElecFiles)
    Total["ElecError"] = Total["ElecError"] / len(ElecFiles)

    for File in ProtFiles:
        Res = readGrasCsv1Tile(path + "/Results/" + File)
        Total["ProtMeV"] += Res["Dose"]
        Total["ProtError"] += Res["Error"]

    Total["ProtMeV"] = Total["ProtMeV"] / len(ProtFiles)
    Total["ProtError"] = Total["ProtError"] / len(ProtFiles)

    return Total

if __name__ == "__main__":
    Path = "/home/anton/Desktop/triton_work/GRAS-2Mat/Test"

    Results = totalMeVGras(Path)

    print("The total electron energy deposition per number of events is: ", Results["ElecMeV"], " MeV")
    print("with ", 100 * Results["ElecError"] / Results["ElecMeV"], " % error")

    print("The total proton energy deposition per number of events is: ", Results["ProtMeV"], " MeV")
    print("with ", 100 * Results["ProtError"] / Results["ProtMeV"], " % error")
