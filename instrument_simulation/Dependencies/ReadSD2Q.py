import numpy as np


def readSDQ2(fileName):
    print(fileName)
    f = open(fileName, "r")

    ShieldoseTable = []
    ReadFlag = 0

    for line in f:
        if ReadFlag == 0:
            if "Dose in Si" in line:
                ReadFlag = 1
        elif ReadFlag == 1:
            if "End of File" in line:
                ReadFlag = 2
            else:
                ShieldoseTable.append(line)

    leng = len(ShieldoseTable)

    # SDDataCollumns = ['Aluminium Thickness', 'Total Dose', 'Electrons', 'Bremsstrahlung', 'Protons']
    SDData = np.zeros((len(ShieldoseTable), 5), dtype=float)

    for i in range(leng):
        SDData[i] = np.fromstring(ShieldoseTable[i], dtype=float, count=5, sep=',')

    return SDData


if __name__ == "__main__":
    print(readSDQ2("Shieldose/Slab/spenvis_sqo.txt"))
