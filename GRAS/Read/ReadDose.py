import numpy as np
import matplotlib.pyplot as plt


def readDose(file):

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
                    # Split the line into parts based on the comma
                    values = line.split(',')
                    TID['dose'].append(float(values[0]))
                    TID['error'].append(float(values[1]))
                    TID['entries'].append(int(float(values[2])))
                    TID['non-zeros'].append(int(float(values[3])))

    # Convert lists inside the dicts to numpy arrays
    for key in keys:
        TID[key] = np.array(TID[key])

    # The default dose unit in GRAS is rad.
    # Most spectra are flux spectra. The resulting dose is therefore given in rad/s.
    # --> divide by 1000 to get to krad.
    radToKrad = 1 / 1000

    for key in ['dose', 'error']:
        TID[key] = TID[key] * radToKrad

    return TID  # 4xNumTiles matrix
    # data[0]: kRad per particle for each volume
    # data[1]: Absolute Error in kRad per particle
    # data[2]: Number of Entries
    # data[3]: Number of non zero entries


if __name__ == "__main__":
    File = "/l/triton_work/ShieldingCurves/Carrington/CarringtonElectronDiffPowTabelated-10mm/Res/TID_757480_117187.csv"

    Results = readDose(File)

    print(np.shape(Results))
    #print(Results)

    # Plot the dose with error bars
    plt.figure(0)
    plt.errorbar(np.arange(len(Results['dose'])), Results['dose'], yerr=Results['error'], fmt=' ', capsize=5, elinewidth=1, capthick=1)
    plt.title('Dose per tile')
    plt.xlabel('Tile number')
    plt.ylabel('Dose [kRad]')
    plt.yscale('log')
    plt.minorticks_on()
    plt.grid(axis='x', which='both')
    plt.grid(axis='y', which='major')

    # Plot the relative error
    plt.figure(1)
    plt.plot(100 * Results['error'] / Results['dose'], '.')
    plt.title('Relative Error in %')
    plt.xlabel('Tile number')
    plt.ylabel('Relative Error [%]')
    plt.grid(which='both')

    # Plot the number of non-zero entries
    plt.figure(2)
    plt.plot(Results['non-zeros'], '.')
    # Add horizontal line at 1
    plt.axhline(y=1, color='r', linestyle='--')
    plt.title('Number of non-zero entries')
    plt.xlabel('Tile number')
    plt.ylabel('Number of non-zero entries')
    plt.yscale('log')
    plt.minorticks_on()
    plt.grid(axis='x', which='both')
    plt.grid(axis='y', which='major')
    plt.show()

