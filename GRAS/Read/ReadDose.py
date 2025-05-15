import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def readDose(file):
    """ 
    Reads TID values from a GRAS CSV output file, supporting both single and multiple volume formats.

    Args:
        file (str): Path to the GRAS CSV file.
    
    Returns:
        dict: Dictionary with numpy arrays of dose info:
            - 'dose': Dose in kRad
            - 'error': Absolute error in kRad
            - 'entries': Number of entries
            - 'non-zeros': Number of non-zero entries
    """
    keys = ['dose', 'error', 'entries', 'non-zeros']
    TID = {key: [] for key in keys}
    
    ReadFlag = 0
    rad_to_krad = 1 / 1000  # Convert rad to kRad

    with open(file, 'r') as f:
        for line in f:
            if ReadFlag == 0:
                if "'TOTAL DOSE FOR INDIVIDUAL VOLUMES'" in line:
                    ReadFlag = 1  # Found the multi-volume block
            elif ReadFlag == 1:
                if "'Non zero entries'" in line:
                    ReadFlag = 2  # Prepare to read values
            elif ReadFlag == 2:
                if "'End of Block'" in line:
                    break  # Done reading relevant block
                else:
                    # Parse values
                    values = [value.strip() for value in line.split(',')]
                    TID['dose'].append(float(values[0]))
                    TID['error'].append(float(values[1]))
                    TID['entries'].append(int(float(values[2])))
                    TID['non-zeros'].append(int(float(values[3])))


    # Fallback to single-volume block if multi not found
    if ReadFlag == 0:
        with open(file, 'r') as f:
            for line in f:
                if "'Non zero entries'" in line:
                    ReadFlag = 1
                elif ReadFlag == 1 and "'End of Block'" not in line:
                    values = [v.strip() for v in line.split(',')]
                    if len(values) >= 4:
                        try:
                            TID['dose'] = [float(values[0])]
                            TID['error'] = [float(values[1])]
                            TID['entries'] = [int(float(values[2]))]
                            TID['non-zeros'] = [int(float(values[3]))]
                            break
                        except ValueError:
                            continue

    # Convert to numpy arrays
    for key in keys:
        TID[key] = np.array(TID[key])

    # Apply rad -> kRad conversion
    TID['dose'] *= rad_to_krad
    TID['error'] *= rad_to_krad

    return TID



def format_number(n):
    return f"{n:.6f}" if isinstance(n, float) or isinstance(n, np.floating) else str(n)


if __name__ == "__main__":
    path = "/l/triton_work/RadEx/1Tile/Res/"

    print("Reading in all files in folder:", path)

    # Get list of all CSV files in the path
    Files = [f for f in os.listdir(path) if f.endswith(".csv")]

    if not Files:
        sys.exit("ERROR !!! No files found")

    print("\n")
    # Header
    print(f"{'File':<25} {'|':<3} {'Dose (kRad/month)':<20} {'|':<3} {'Error (kRad/month)':<20} {'|':<3} {'Entries':<15} {'|':<3} {'Non Zero Entries':<25} {'|':<3} {'Rel. Error (%)':<20}")
    print("-" * 150)

    for File in Files:
        full_path = os.path.join(path, File)
        Results = readDose(full_path)

        # Multiply dose and error by number of seconds in a month
        # time_scale = 30 * 24 * 60 * 60  # seconds in a month
        # Results['dose'] *= time_scale
        # Results['error'] *= time_scale

        # Aggregate statistics
        dose_total = Results['dose'].sum()
        error_total = Results['error'].sum()
        entries_total = Results['entries'].sum()
        non_zeros_total = Results['non-zeros'].sum()

        # Avoid division by zero
        if dose_total != 0:
            relative_dose_error_percent = (error_total / dose_total) * 100
        else:
            relative_dose_error_percent = 0

        print(f"{File:<25} {'|':<3} {format_number(dose_total):<20} {'|':<3} {format_number(error_total):<20} {'|':<3} {entries_total:<15} {'|':<3} {non_zeros_total:<25} {'|':<3} {relative_dose_error_percent:<20.2f}")


    # Num of tiles
    NumTiles = len(Results['dose'])
    print(f"Number of tiles: {NumTiles}")

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

    # plt.savefig(PlotPath + '/Dose.pdf', format='pdf', bbox_inches="tight")

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

    # plt.savefig(PlotPath + '/Error.pdf', format='pdf', bbox_inches="tight")

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

    # plt.savefig(PlotPath + '/NonZeros.pdf', format='pdf', bbox_inches="tight")
    
    plt.show()

