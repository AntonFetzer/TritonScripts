import numpy as np
import matplotlib.pyplot as plt
import csv


def readGrasCsv(file):
    data = [[], [], [], []]

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
                    data[2].append(float(line[2]))
                    data[3].append(float(line[3]))

    #Res = np.ndarray(np.shape(data), dtype=np.float64)
    Res = np.array(data, dtype=np.float64)

    #for c, column in enumerate(data):
    #    for i, item in enumerate(column):
    #       Res[c, i] = item

    return Res  # 2xNumTiles matrix
    # data[0]: MeV per particle for each volume
    # data[1]: Absolute Error on MeV per particle
    # data[:, x]: Tuple of MeVs with corresponding error


if __name__ == "__main__":
    File = "/home/anton/Desktop/triton_work/ShieldingCurves/100Al/Res2e5/Electrons_86777_662.csv"

    Results = readGrasCsv(File)

    print(np.shape(Results))
    print(Results)

    plt.plot(Results[3])
    #plt.yscale("log")
    plt.show()
