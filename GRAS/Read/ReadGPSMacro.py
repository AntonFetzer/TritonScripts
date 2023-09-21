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
    file1 = "/l/triton_work/Spectra/ISO-GTO/ISO-GTO-Fe-mission.mac"
    Data1 = readGPSMacro(file1)

    file2 = "/l/triton_work/Spectra/ISO-GTO/ISO-GTO-Fe.mac"
    Data2 = readGPSMacro(file2)

    plt.plot(Data1[0], Data1[1], label="ISO-GTO-Fe-mission")
    plt.plot(Data2[0], Data2[1], label="ISO-GTO-Fe", linestyle="--")

    plt.yscale("log")
    plt.xscale("log")
    plt.grid(which="both")
    plt.legend()
    plt.title("Differential particle flux")
    plt.xlabel("Kinetic energy [MeV]")
    plt.ylabel("Differential Flux [cm-2 s-1 MeV-1]")

    plt.show()
    # plt.savefig("/l/triton_work/Spectra/ISO-GTO/MacroComparison.pdf", format='pdf', bbox_inches="tight")
