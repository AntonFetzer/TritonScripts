import numpy as np
import matplotlib.pyplot as plt
import csv


def readSpenvis_tri(file):

    PE = []
    PI = []
    PD = []
    EE = []
    EI = []
    ED = []
    ReadFlag = 0

    print("Reading in File: " + file)
    with open(file, 'r') as f:
        reader = csv.reader(f)
        for line in reader:
            if ReadFlag == 0 and "'Differential Flux'" in line:
                ReadFlag = 1
            elif ReadFlag == 1:
                if "'End of Block'" in line:
                    ReadFlag = 2
                else:
                    PE.append(float(line[0]))
                    PI.append(float(line[1]))
                    PD.append(float(line[2]))
            elif ReadFlag == 2 and "'Differential Flux'" in line:
                ReadFlag = 3
            elif ReadFlag == 3:
                if "'End of File'" in line:
                    ReadFlag = 4
                else:
                    EE.append(float(line[0]))
                    EI.append(float(line[1]))
                    ED.append(float(line[2]))

    ProtonData = np.asarray([PE, PI, PD])
    ElectronData = np.asarray([EE, EI, ED])

    return ProtonData, ElectronData



if __name__ == "__main__":

    file = "/l/triton_work/SuperGTO/spenvis_tri.txt"

    Protons, Electrons = readSpenvis_tri(file)

    print(np.shape(Protons))

    for l in range(np.shape(Protons)[1]):
        if Protons[1, l] == 0 or Protons[2, l] == 0:
            Protons = np.delete(Protons, l, 1)

    for l in range(np.shape(Electrons)[1]):
        if Electrons[1, l] == 0 or Electrons[2, l] == 0:
            Electrons = np.delete(Electrons, l, 1)


    plt.plot(Protons[0], Protons[1], label="Trapped Protons")
    #plt.plot(Protons[0], Protons[2], label="Proton Differential")
    plt.plot(Electrons[0], Electrons[1], label="Trapped Electrons")
    #plt.plot(Electrons[0], Electrons[2], label="Electron Differential")



    plt.yscale("log")
    plt.xscale("log")
    plt.grid(which="both")

    plt.title("AP-9/AE-9 and SAPPHIRE particle flux on Super-GTO")
    plt.xlabel("Kinetic energy [MeV]")
    plt.ylabel("Integral Flux [cm-2 s-1]")
    #plt.ylabel("Differential Flux [cm-2 s-1 MeV-1]")

    plt.show()
