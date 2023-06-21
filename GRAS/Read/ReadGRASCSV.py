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

    Res = np.ndarray(np.shape(data), dtype=np.float64)
    #Res = np.array(data, dtype=np.float64)

    # Dose data is in rad/s --> multiply with number of seconds in a month to get to dose per months.
    # Dose is given per generated particle --> need to divide by the number of files
    # Dose is given in rad --> divide by 1000 to get to krad.
    ScaleFactor = 30 * 24 * 60 * 60 / 1000

    for i, tile in enumerate(data[0]):
        Res[0][i] = data[0][i] * ScaleFactor
        Res[1][i] = data[1][i] * ScaleFactor
        Res[2][i] = data[2][i]
        Res[3][i] = data[3][i]

    return Res  # 4xNumTiles matrix
    # data[0]: kRad per particle for each volume
    # data[1]: Absolute Error in kRad per particle
    # data[2]: Number of Entries
    # data[3]: Number of non zero entries


if __name__ == "__main__":
    File = "/home/anton/Desktop/triton_work/RadEx/RadEx0mm/Res/Electrons_981978_151919.csv"

    Results = readGrasCsv(File)

    #print(np.shape(Results))
    print(Results)

    plt.plot(100 * Results[1]/Results[0], '.')
    plt.title("Relative Error in %")
    #plt.yscale("log")
    plt.show()
