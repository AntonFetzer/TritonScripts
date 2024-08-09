import os
import numpy as np
from GRAS.Read.ReadDose import readDose
import matplotlib.pyplot as plt
import sys


def totalDose(path):
    """
    Reads the TID values of all sensitive volumes of all GRAS CSV output files in a folder and calculates the total dose for each of the tiles.

    Args:
        path (str): The path to the folder containing the CSV files.
    
    Returns:
        dict: A dictionary containing the following keys:
            - 'dose': A numpy array containing the total dose values in kRad for each tile.
            - 'error': A numpy array containing the absolute total dose error in kRad for each tile.
            - 'entries': A numpy array containing the number of entries for each tile.
            - 'non-zeros': A numpy array containing the number of non-zero entries for each tile.    
    """

    print("\nReading in all csv files in folder:", path)

    # Get list of all csv files in Path
    Files = [f for f in os.listdir(path) if '.csv' in f]

    if not Files:
        # Raise an error if no files are found
        sys.exit("ERROR !!! No csv files found in folder: " + path)

    keys = ['dose', 'error', 'entries', 'non-zeros']
    TID = {key: [] for key in keys}
    # Ionising dose in kRad
    # Number of entries is the same for each tile, because each detected particle causes an entry at all tiles. Most of these entries are therefore empty.
    # Number of non-zero entries is different for each tile.

    # RawData is a list of dictionaries.
    # Each dictionary contains the dose, error, entries and non-zero entries for each file.
    RawData = [readDose(os.path.join(path, file)) for file in Files]

    # Initialize numpy arrays for dose, error, entries, and non-zeros using the first file to get the correct shape.
    for key in keys:
        TID[key] = np.zeros_like(RawData[0][key])

    for raw in RawData:
        #print(raw['dose'] * raw['entries'])
        TID['dose'] += raw['dose'] * raw['entries']  # Weighted average of all the dose results. Multiply dose by number of entries, then later divide by the total number of entries.
        TID['entries'] += raw['entries']  # Summed up number of entries from all files
        TID['non-zeros'] += raw['non-zeros']  # Summed up number of Non Zero Entries from all files

    # Calculate the weighted average of the dose.
    TID['dose'] = TID['dose'] / TID['entries']  # Divide by summed up number of entries to get the weighted average of the dose results

    # Error estimation based on the standard deviation of the population and distance from the mean.
    ErrTerm = np.zeros_like(TID['error'])
    for raw in RawData:
        ErrTerm += raw['entries'] ** 2 * (raw['error'] ** 2 + (raw['dose'] - TID['dose']) ** 2)

    TID['error'] = np.sqrt(ErrTerm) / TID['entries']

    # Entries is the same for all tiles, because each tile gets an entry if a particle is detected in any tile.
    TotalEntries = TID['entries'][0]
    # Check if the number of entries is the same for all tiles.
    if not np.all(TID['entries'] == TID['entries'][0]):
        sys.exit("ERROR !!! Number of entries is not the same for all tiles")

    NumTiles = len(TID['dose'])
    MinNZE = min(TID['non-zeros'])   # Minimum number of non-zero entries
    LowestTile = np.argmin(TID['non-zeros'])
    ParticlesPerTile = TID['entries'][0] / NumTiles

    print("Total number of particles:", f"{TotalEntries:.3g}")
    print("Average number of Particles per tile:", f"{ParticlesPerTile:.3g}")
    print("Minimum number of non Zero Entries:", f"{MinNZE:.3g}", "at Tile number", str(LowestTile))


    if MinNZE < 100:
        print("Warning !!! Tile", LowestTile, "has only", MinNZE, "Non-Zero entries !!!")

    RelativeError = TID['error'] / TID['dose']
    MaxRelativeError = max(RelativeError)
    MaxRelativeErrorTile = np.argmax(RelativeError)
    print("Maximum relative error =", f"{MaxRelativeError * 100:.2f} % at Tile number {MaxRelativeErrorTile}")

    PartNumRequired = (MaxRelativeError / 0.01) ** 2 * TotalEntries
    ParticleMultipleRequired = PartNumRequired / TotalEntries

    print(f"The number of particles required to achieve 1% error is: {PartNumRequired:.3g} or {ParticleMultipleRequired:.3g} times the number of particles in the run")

    if MaxRelativeError > 1:
        print(f"Warning !!! Tile {MaxRelativeErrorTile} has {MaxRelativeError * 100:.2f} % relative error!!!")

    # Print datatypes and shapes of the arrays
    #for key in keys:
    #    print(key, type(TID[key]), TID[key].shape)
        # print data type of the array entries
    #    print(key, TID[key].dtype)
        
    return TID


if __name__ == "__main__":
    Path = "/l/triton_work/Shielding_Curves/Carrington/GEO-AP9-mission/Res/"

    Results = totalDose(Path)

    NumTiles = len(Results['dose'])

    # Multiply dose and error with the number of seconds in a month
    # to get the dose in kRad per month
    # Results['dose'] *= 60 * 60 * 24 * 30
    # Results['error'] *= 60 * 60 * 24 * 30

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

    plt.savefig(Path + "..//Dose.pdf", format='pdf', bbox_inches="tight")

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

    plt.savefig(Path + "../Error.pdf", format='pdf', bbox_inches="tight")

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

    plt.savefig(Path + "../NonZeros.pdf", format='pdf', bbox_inches="tight")

    #plt.show()



