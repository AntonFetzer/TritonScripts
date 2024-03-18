import numpy as np
import matplotlib.pyplot as plt


def readGRASdose(file):

    keys = ['dose', 'error', 'entries', 'non-zeros']

    # Generating the dictionary from the list of keys
    TID = {key: [] for key in keys}

    ReadFlag = 0

    # Open the file and read line by line
    with open(file, 'r') as f:
        for line in f:
            # Process the file based on the current state of ReadFlag
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
                    values = [float(x) for x in line.split(',')]
                    for i, key in enumerate(keys):
                        TID[key].append(values[i])

    # Convert lists inside the dicts to numpy arrays
    for key in keys:
        TID[key] = np.array(TID[key])

    # Dose data is in rad/s --> multiply with number of seconds in a month to get to dose per months.
    # Dose is given per generated particle --> need to divide by the number of files
    # Dose is given in rad --> divide by 1000 to get to krad.
    ScaleFactor = 30 * 24 * 60 * 60 / 1000

    for key in ['dose', 'error']:
        TID[key] = TID[key] * ScaleFactor

    return TID  # 4xNumTiles matrix
    # data[0]: kRad per particle for each volume
    # data[1]: Absolute Error in kRad per particle
    # data[2]: Number of Entries
    # data[3]: Number of non zero entries


if __name__ == "__main__":
    File = "/l/triton_work/RadEx/RadEx0mm/Res/Electrons_981978_151919.csv"

    Results = readGRASdose(File)

    #print(np.shape(Results))
    print(Results)

    plt.plot(100 * Results[1]/Results[0], '.')
    plt.title("Relative Error in %")
    #plt.yscale("log")
    plt.show()
