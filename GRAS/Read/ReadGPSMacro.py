import numpy as np
import matplotlib.pyplot as plt
import csv


def readGPSMacro(file):

    Energy = []
    DiffFlux = []
    ReadFlag = 0

    print("Reading in File: " + file)
    with open(file, 'r') as f:
        reader = csv.reader(f)
        for line in reader:
            #print(line)
            if ReadFlag == 0:
                if '/gps/hist/type arb' in line:
                    ReadFlag = 1
                    #print(line)
            elif ReadFlag == 1:
                if '/gps/hist/inter Log' in line:
                    ReadFlag = 2
                    #print(line)
                else:
                    Text = line[0]
                    Text = Text.split()
                    #print(float(Text[1]), float(Text[2]))
                    Energy.append(float(Text[1]))
                    DiffFlux.append(float(Text[2]))

    return np.asarray([Energy, DiffFlux])



if __name__ == "__main__":

    file = "/home/anton/Desktop/triton_work/SuperGTO/SuperGTOElectrons.mac"
    Electrons = readGPSMacro(file)

    file = "/home/anton/Desktop/triton_work/SuperGTO/SuperGTOProtons.mac"
    Protons = readGPSMacro(file)

    plt.plot(Protons[0], Protons[1], label="Proton Differential")
    plt.plot(Electrons[0], Electrons[1], label="Electron Differential")

    plt.yscale("log")
    plt.xscale("log")
    plt.grid(which="both")
    plt.legend()
    plt.title("AP-9 and AE-9 Super-GTO spectra")
    plt.xlabel("Kinetic energy [MeV]")
    plt.ylabel("Differential Flux [cm-2 s-1 MeV-1]")

    #plt.show()
    plt.savefig("/home/anton/Desktop/triton_work/SuperGTO/SpectraForGeant4.eps", format='eps', bbox_inches="tight")
