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

    file = "/home/anton/Desktop/triton_work/Spectra/A9/AE9/AE9500keV.mac"
    Electrons = readGPSMacro(file)

    file = "/home/anton/Desktop/triton_work/Spectra/A9/AP9/AP910MeV.mac"
    Protons = readGPSMacro(file)

    plt.plot(Electrons[0], Electrons[1], label="AE-9 Electron Flux")
    plt.plot(Protons[0], Protons[1], label="AP-9 Proton Flux")

    plt.yscale("log")
    plt.xscale("log")
    plt.grid(which="both")
    plt.legend()
    plt.title("Differential AP-9 and AE-9 spectra on GTO")
    plt.xlabel("Kinetic energy [MeV]")
    plt.ylabel("Differential Flux [cm-2 s-1 MeV-1]")

    #plt.show()
    plt.savefig("/home/anton/Desktop/TritonPlots/Paper/SpectraBasic.pdf", format='pdf', bbox_inches="tight")
