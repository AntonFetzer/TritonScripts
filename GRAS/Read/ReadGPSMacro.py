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
            # print(line)
            if ReadFlag == 0:
                if '/gps/hist/type arb' in line:
                    ReadFlag = 1
                    # print(line)
            elif ReadFlag == 1:
                if "/gps/hist/inter " in line[0]:
                    ReadFlag = 2
                    # print(line)
                else:
                    Text = line[0]
                    Text = Text.split()
                    # print(float(Text[1]), float(Text[2]))
                    Energy.append(float(Text[1]))
                    DiffFlux.append(float(Text[2]))

    return np.asarray([Energy, DiffFlux])


if __name__ == "__main__":
    files = [
        "/l/triton_work/Spectra/Carrington/Electron/CarringtonElectronDiffPowTabelated.mac",
        "/l/triton_work/Spectra/A9-GTO/AE9/AE9Mission.mac",
        "/l/triton_work/Spectra/FS1/A9-FS1/AE9-FS1-mission.mac"
        ]
    
    labels = [
        "Carrington Electron Flux",
        "AE9-GTO",
        "AE9-LEO"
        ]
    
    for file in files:
        Data = readGPSMacro(file)
        plt.plot(Data[0], Data[1], label=labels[files.index(file)], marker=".")

    plt.yscale("log")
    plt.xscale("log")
    plt.grid(which="both")
    plt.legend()
    plt.title("Differential Electron Flux Comparison")
    plt.xlabel("Kinetic energy [MeV]")
    plt.ylabel("Differential Flux [cm-2 s-1 MeV-1]")
    plt.xlim(0.1, 20)
    plt.ylim(1e1, 1e13)

    #plt.show()
    plt.savefig("/u/02/fetzera1/unix/Desktop/TritonPlots/CarringtonPaper/ElectronFluxComparison.pdf", format='pdf', bbox_inches="tight")
