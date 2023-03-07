import numpy as np
import matplotlib.pyplot as plt

def readSDQ2(fileName):
    print("Reading in", fileName)
    f = open(fileName, "r")

    SDData = []
    ReadFlag = 0

    for line in f:

        if ReadFlag == 0:
            if "Dose in Si" in line:
                ReadFlag = 1
        elif ReadFlag == 1:
            if "End of File" in line:
                ReadFlag = 2
            else:
                SDData.append(np.fromstring(line, dtype=np.float64, sep=','))

    # SDDataCollumns = ['Aluminium Thickness', 'Total Dose', 'Electrons', 'Bremsstrahlung', 'Protons']

    return np.asarray(SDData)


if __name__ == "__main__":
    Data = readSDQ2("/home/anton/Desktop/triton_work/Spectra/A9/Shieldose/spenvis_sqo.txt")

    print(np.shape(Data))

    print("Aluminium Thicknesses", Data[:, 0])
    print("Total Dose", Data[:, 1])
    print("Electrons", Data[:, 2])
    print("Bremsstrahlung", Data[:, 3])
    print("Protons", Data[:, 4])

    plt.plot(Data[:, 0], Data[:, 2] / 1000, '.', label="SHIELDOSE-2Q Electrons")
    plt.plot(Data[:, 0], Data[:, 3] / 1000, '.', label="SHIELDOSE-2Q Bremsstrahlung")
    plt.plot(Data[:, 0], Data[:, 4] / 1000, '.', label="SHIELDOSE-2Q Protons")
    plt.plot(Data[:, 0], Data[:, 1] / 1000, '.', label="SHIELDOSE-2Q Total Dose")

    plt.yscale("log")
    # plt.xscale("log")
    plt.grid(which='major')
    plt.title("Shielddose")
    plt.xlabel("Aluminium Shield Thickness [mm]")
    plt.ylabel("Total ionizing dose in silicon [krad]")
    plt.legend()

    plt.show()