import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

def readDoseSingle(file):
    data = pd.DataFrame(columns=['dose', 'error', 'entries', 'Non zero entries'])

    read_values = False

    with open(file, 'r') as f:
        for line in f:
            line = line.strip()

            if not read_values and "'Non zero entries'" in line:
                read_values = True
            elif read_values:
                if "'End of Block'" in line:
                    break

                values = line.split(',')
                data.at[0, 'dose'] = float(values[0])
                data.at[0, 'error'] = float(values[1])
                data.at[0, 'entries'] = float(values[2])
                data.at[0, 'Non zero entries'] = float(values[3])

    # Dose data is in rad/s --> multiply with number of seconds in a month to get to dose per months.
    # Dose is given per generated particle --> need to divide by the number of files
    # Dose is given in rad --> divide by 1000 to get to krad.
    # scale_factor = 30 * 24 * 60 * 60 / 1000
    scale_factor = 1 / 1000

    # Use in-place operations to update the data array
    data['dose'] *= scale_factor
    data['error'] *= scale_factor

    return data  # 4xNumTiles matrix
    # data[0]: kRad per particle for each volume
    # data[1]: Absolute Error in kRad per particle
    # data[2]: Number of Entries
    # data[3]: Number of non zero entries


if __name__ == "__main__":
    path = "/l/triton_work/RadEx/RadEx1mm/Res/"
    #particle = "Elec"
    particle = "Prot"

    print("Reading in all", particle, "files in folder:", path)

    # Get list of all csv files in Path
    Files = [f for f in os.listdir(path) if particle in f]

    if not Files:
        sys.exit("ERROR !!! No files found")

    RawData = []

    for File in Files:
        Results = readDoseSingle(path + File)

    print(Results)

    print("Total number of entries:", Results['entries'][0])
    print("Number of non Zero Entries:", Results['Non zero entries'][0])
    print("Dose Datatype:", type(Results['dose']))


    relative_dose_error_percent = (Results['error'][0] / Results['dose'][0]) * 100

    ax1 = Results[['dose']].plot(kind='bar', yerr=Results['error'], capsize=5, logy=True, legend=False)
    ax1.set_ylabel('Dose (krad)')
    ax1.set_title(f'Dose with Error (Relative Error: {relative_dose_error_percent:.2f}%)')
    ax1.set_xticklabels(['dose'], rotation=0)

    # Reshape the data for side-by-side bar plot
    reshaped_data = Results.melt(value_vars=['entries', 'Non zero entries'])
    ax2 = reshaped_data.plot(x='variable', y='value', kind='bar', logy=True, legend=False)
    ax2.set_ylabel('Count (log scale)')
    ax2.set_title('Entries and Non Zero Entries')
    ax2.set_xticklabels(reshaped_data['variable'], rotation=0)
    ax2.set_xlabel('')

    plt.show()
