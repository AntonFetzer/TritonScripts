import numpy as np
import matplotlib.pyplot as plt
from GRAS.Dependencies.TotalDose import totalDose

def mergeTotalDose(List_of_Dose_Dicts: list) -> dict:
    """ 
    Merges the dose dicts in List_of_Dose_Dicts into one dose dict

    Args:
        List_of_Dose_Dicts (list): A list containing dictionaries of the format:
            - 'dose': A numpy array containing the total dose values in kRad for each tile.
            - 'error': A numpy array containing the absolute total dose error in kRad for each tile.
            - 'entries': A numpy array containing the number of entries for each tile.
            - 'non-zeros': A numpy array containing the number of non-zero entries for each tile.  
    
    Returns:
        dict: A dictionary containing the following keys:
            - 'dose': A numpy array containing the total dose values in kRad for each tile.
            - 'error': A numpy array containing the absolute total dose error in kRad for each tile.
            - 'entries': A numpy array containing the number of entries for each tile.
            - 'non-zeros': A numpy array containing the number of non-zero entries for each tile.

    Excepts:
        If the input is not a list of dicts, a TypeError is raised.
        If the number of tiles in the dose dictionaries is not the same, a ValueError is raised.

    """

    # Check if input is a list of dicts
    if not all(isinstance(d, dict) for d in List_of_Dose_Dicts):
        raise TypeError("Input must be a list of dictionaries")
    
    keys = ['dose', 'error', 'entries', 'non-zeros']

    # Check if all dicts in the list have the correct keys
    if not all(all(k in d for k in keys) for d in List_of_Dose_Dicts):
        raise ValueError("All dictionaries in the list must have the keys 'dose', 'error', 'entries', and 'non-zeros'")
    
    # Check if all numpy arrays that are the entries of the dicts have the same lenght
    if not all(all(len(d[k]) == len(List_of_Dose_Dicts[0][k]) for k in keys) for d in List_of_Dose_Dicts):
        raise ValueError("All dictionaries in the list must have the same number of tiles")
    
    # Generating the results dictionary from the list of keys
    Total = {key: None for key in keys}

    # Initialise dict for the TotalDose using the first entry of the list
    for key in keys:
        Total[key] = np.zeros_like(List_of_Dose_Dicts[0][key])

    # Calculate weighted average dose of all datasets using the number of entries
    for d in List_of_Dose_Dicts:
        Total['dose'] += d['dose']
        Total['error'] += np.square(d['error'])
        Total['entries'] += d['entries']
        Total['non-zeros'] += d['non-zeros']

    # Calculate the weighted average of the dose and error
    Total['dose'] = Total['dose']  # Divide by summed up number of entries to get the weighted average of the dose results
    Total['error'] = np.sqrt(Total['error'])

    return Total


if __name__ == "__main__":
    Elec_path = "/l/triton_work/RadEx/RadEx/LEO-Electrons"
    Prot_path = "/l/triton_work/RadEx/RadEx/LEO-Protons"

    Elec = totalDose(Elec_path)
    Prot = totalDose(Prot_path)

    # Multiply dose and error with the number of seconds in a month to get the dose in kRad per month
    # Elec['dose'] *= 60 * 60 * 24 * 30
    # Elec['error'] *= 60 * 60 * 24 * 30
    # Prot['dose'] *= 60 * 60 * 24 * 30
    # Prot['error'] *= 60 * 60 * 24 * 30

    Total = mergeTotalDose([Elec, Prot])

    NumTiles = len(Total['dose'])

    # Plot the dose with error bars
    plt.figure(0)
    plt.errorbar(np.arange(NumTiles), Elec['dose'], yerr=Elec['error'], fmt=' ', capsize=5, elinewidth=1, capthick=1, label='Electron')
    plt.errorbar(np.arange(NumTiles), Prot['dose'], yerr=Prot['error'], fmt=' ', capsize=5, elinewidth=1, capthick=1, label='Proton')
    plt.errorbar(np.arange(NumTiles), Total['dose'], yerr=Total['error'], fmt=' ', capsize=5, elinewidth=1, capthick=1, label='Total')
    # Add horizontal line at 1 kRad
    plt.axhline(y=1, color='r', linestyle='--', label='1 kRad')
    plt.title('Dose per tile')
    plt.xlabel('Tile number')
    plt.ylabel('Dose [kRad]')
    plt.yscale('log')
    plt.grid(which='both')
    plt.legend()

    plt.show()

    # Print Results as comma seperated table with collumns for total dose, total error, electron dose, electron error, proton dose and proton error
    print("Tile, Total Dose, Total Error, Electron Dose, Electron Error, Proton Dose, Proton Error")
    for i in range(NumTiles):
        print(f"{i+1}, {Total['dose'][i]}, {Total['error'][i]}, {Elec['dose'][i]}, {Elec['error'][i]}, {Prot['dose'][i]}, {Prot['error'][i]}")