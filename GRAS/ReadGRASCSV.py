import matplotlib.pyplot as plt
import numpy as np
import os
import csv


def readGrasCsv(file):
    data = {}  # Initialise empty dict !

    ReadFlag = 0
    # print("Reading in File: " + file)
    with open(file, 'r') as f:
        reader = csv.reader(f)
        for line in reader:
            # print(line)
            if ReadFlag == 0:
                if "'Non zero entries'" in line:
                    ReadFlag = 1
            elif ReadFlag == 1:
                if "'End of Block'" in line:
                    ReadFlag = 2
                else:
                    data["Dose"] = float(line[0])
                    data["Error"] = float(line[1])
    return data


if __name__ == "__main__":
    Path = "/home/anton/Desktop/triton_work/GRAS-2Mat/Test"

    # Get list of all csv files in Path
    ElecFiles = [f for f in os.listdir(Path + "/Results/") if "Electrons" in f]
    ProtFiles = [f for f in os.listdir(Path + "/Results/") if "Protons" in f]

    ElecMeV = 0
    ElecError = 0

    ProtMeV = 0
    ProtError = 0

    for File in ElecFiles:
        Res = readGrasCsv(Path + "/Results/" + File)
        ElecMeV += Res["Dose"]
        ElecError += Res["Error"]

    ElecMeV = ElecMeV / len(ElecFiles)
    ElecError = ElecError / len(ElecFiles)

    print("The total electron energy deposition per number of events is: " + str(ElecMeV) + " MeV")
    print("with " + str(100 * ElecError / ElecMeV) + " % error")

    for File in ProtFiles:
        Res = readGrasCsv(Path + "/Results/" + File)
        ProtMeV += Res["Dose"]
        ProtError += Res["Error"]
    
    ProtMeV = ProtMeV / len(ProtFiles)
    ProtError = ProtError / len(ProtFiles)

    print("The total proton energy deposition per number of events is: " + str(ProtMeV) + " MeV")
    print("with " + str(100 * ProtError / ProtMeV) + " % error")
