import numpy as np
import matplotlib.pyplot as plt
import csv


def readGrasCsv(file):
    data = [[], []]

    ReadFlag = 0
    # print("Reading in File: " + file)
    with open(file, 'r') as f:
        reader = csv.reader(f)
        for line in reader:
            # print(line)
            if ReadFlag == 0:
                if "'TOTAL DOSE FOR INDIVIDUAL VOLUMES'" in line:
                    ReadFlag = 1
            elif ReadFlag == 1:
                if "'Non zero entries'" in line:
                    ReadFlag = 2
            elif ReadFlag == 2:
                if "'End of Block'" in line:
                    break
                else:
                    data[0].append(float(line[0]))
                    data[1].append(float(line[1]))
    return np.asarray(data)  # 2xNumTiles matrix
    # data[0]: MeV per particle for each volume
    # data[1]: Absolute Error on MeV per particle
    # data[:, x]: Tuple of MeVs with corresponding error


if __name__ == "__main__":
    File = "/home/anton/Desktop/triton_work/2LayerOpt/Al-Pb/Res/ElectronsA_966954_65703.csv"

    Results = readGrasCsv(File)

    print(np.shape(Results))
    print(Results)

    plt.plot(Results[0])
    plt.show()
