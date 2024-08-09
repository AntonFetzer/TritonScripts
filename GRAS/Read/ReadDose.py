import numpy as np
import matplotlib.pyplot as plt


def readDose(file):
    """ 
    Reads the TID values of all sensitive volumes of a GRAS CSV output file.

    Args:
        file (str): The path to the CSV file.
    
    Returns:
        dict: A dictionary containing the following keys:
            - 'dose': A numpy array containing the dose values in kRad
            - 'error': A numpy array containing the absolute dose error in kRad
            - 'entries': A numpy array containing the number of entries.
            - 'non-zeros': A numpy array containing the number of non-zero entries
    """

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
    # --> divide by 1000 to get to krad.
    radToKrad = 1 / 1000

    for key in ['dose', 'error']:
        TID[key] = TID[key] * radToKrad

    return TID


if __name__ == "__main__":
    File = "/l/triton_work/Shielding_Curves/Carrington/VAB-AE9-mission/Res/TID_891597_137936.csv"

    # Plot path is the parent folder of the file path
    PlotPath = File.split('/')
    PlotPath = '/'.join(PlotPath[:-2])

    Results = readDose(File)

    NumTiles = len(Results['dose'])
    # print("Number of tiles: ", NumTiles)

    # print("Results keys: ",Results.keys())
    # print("Results['dose'] type: ", type(Results['dose']))
    # print("Results['dose'] ndarray length: ", len(Results['dose']))

    # Plot the dose with error bars
    plt.figure(0)
    plt.errorbar(np.arange(NumTiles), Results['dose'], yerr=Results['error'], fmt=' ', capsize=5, elinewidth=1, capthick=1, label='Dose')
    # Add horizontal line at 1 kRad
    plt.axhline(y=1, color='r', linestyle='--', label='1 kRad')
    plt.title('Dose per tile')
    plt.xlabel('Tile number')
    plt.ylabel('Dose [kRad]')
    plt.yscale('log')
    plt.grid(which='both')
    plt.legend()

    plt.savefig(PlotPath + '/Dose.png', format='pdf', bbox_inches="tight")

    # Plot the relative error
    plt.figure(1)
    plt.plot(100 * Results['error'] / Results['dose'], '.', label='Relative Error')
    # Add horizontal line at 1%
    plt.axhline(y=1, color='r', linestyle='--', label='1% error')
    plt.title('Relative Error in %')
    plt.xlabel('Tile number')
    plt.ylabel('Relative Error [%]')
    plt.grid(which='both')
    plt.legend()

    plt.savefig(PlotPath + '/Error.png', format='pdf', bbox_inches="tight")

    # Plot the number of non-zero entries
    plt.figure(2)
    plt.plot(Results['non-zeros'], '.', label='Non-zero entries')
    # Add horizontal line at 1
    plt.axhline(y=1, color='r', linestyle='--', label='1 entry')
    plt.title('Number of non-zero entries')
    plt.xlabel('Tile number')
    plt.ylabel('Number of non-zero entries')
    plt.yscale('log')
    plt.grid(axis='x', which='both')
    plt.grid(axis='y', which='major')
    plt.legend()

    plt.savefig(PlotPath + '/NonZeros.png', format='pdf', bbox_inches="tight")
    
    # plt.show()

