import os
import numpy as np
from GRAS.Read.ReadGRASCSVTotal import readGrasCsvTotal
import matplotlib.pyplot as plt
import sys
import pandas as pd

def totalkRadGrasSingle(path):
    print("")
    print("Reading in all csv files in folder:", path)

    # Get list of all csv files in Path
    Files = os.listdir(path)

    if not Files:
        sys.exit("ERROR !!! No files found")

    # Initialize the combined DataFrame
    combined_data = pd.DataFrame(columns=['Dose', 'Error', 'Entries', 'Non zero entries'])

    for file in Files:
        data = readGrasCsvTotal(os.path.join(path, file))
        combined_data = pd.concat([combined_data, data], ignore_index=True)

    # Sum up all the entries in the 'Entries' column
    total_entries = combined_data['Entries'].sum()

    # Sum up all the Non zero entries in the 'Non zero entries' column
    total_Non_zero_entries = combined_data['Non zero entries'].sum()

    # Calculate the weighted average of the dose results
    combined_data['WeightedDose'] = combined_data['Dose'] * combined_data['Entries']
    Dose = combined_data['WeightedDose'].sum() / total_entries


    # The 'ErrTerm' is an intermediate value used to compute the overall error
    # It accounts for both the individual errors and the differences in the dose values
    # This is achieved by the following formula:
    # ErrTerm = Entries^2 * (Error^2 + (Dose - AvgDose)^2)

    # Step 1: Calculate the squared difference between each file's dose and the overall weighted average dose
    squared_dose_diff = (combined_data['Dose'] - Dose) ** 2

    # Step 2: Add squared error of each file's dose to the squared dose difference
    error_components = combined_data['Error'] ** 2 + squared_dose_diff

    # Step 3: Multiply the sum of squared error components by the square of the number of entries for each file
    combined_data['ErrTerm'] = combined_data['Entries'] ** 2 * error_components

    # To find the overall error, we need to sum up all the error terms
    # and divide by the total number of entries

    # Step 4: Sum up all the error terms in the 'ErrTerm' column
    total_error_term = combined_data['ErrTerm'].sum()

    # Step 5: Calculate the overall error by taking the square root of the total error term
    # and dividing by the total number of entries
    Error = np.sqrt(total_error_term) / total_entries

    print(f"Total number of entries: {total_entries:.2}")
    print(f"Number of non Zero Entries: {total_Non_zero_entries:.2}")

    if total_Non_zero_entries < 100:
        print("ERROR !!! Only", total_Non_zero_entries, "Non-Zero entries !!!")

    RelativeError = Error / Dose
    print("Relative error =", str(round(RelativeError * 100, 2)), "%")

    PartNumRequired = (RelativeError/0.01)**2 * total_entries

    print("The number of particles required to achieve 1% err is:", f"{PartNumRequired:.2}", "or", f"{PartNumRequired/total_entries:.2f}", "times the simulated number")

    if RelativeError > 1:
        print("ERROR !!! " + str(RelativeError) + " % relative error!!!")

    # Create a new DataFrame with the results
    Res = pd.DataFrame(data={'Dose': [Dose],
                             'Error': [Error],
                             'Entries': [total_entries],
                             'Non zero entries': [total_Non_zero_entries]})

    return Res



if __name__ == "__main__":
    Path = "/l/triton_work/1Tile/A9-GTO-FullSpectrum/Electrons-6mmAl/Res/"

    Results = totalkRadGrasSingle(Path)

    print(Results)
    '''
    relative_dose_error_percent = (Results['Error'][0] / Results['Dose'][0]) * 100

    ax1 = Results[['Dose']].plot(kind='bar', yerr=Results['Error'], capsize=5, legend=False)
    ax1.set_ylabel('Dose (krad)')
    ax1.set_title(f'Dose with Error (Relative Error: {relative_dose_error_percent:.2f}%)')
    ax1.set_xticklabels(['Dose'], rotation=0)
    ax1.set_ylim(bottom=0)  # Set the lower limit of the y-axis to 0

    # Reshape the data for side-by-side bar plot
    reshaped_data = Results.melt(value_vars=['Entries', 'Non zero entries'])
    ax2 = reshaped_data.plot(x='variable', y='value', kind='bar', legend=False)
    ax2.set_ylabel('Count (log scale)')
    ax2.set_title('Entries and Non Zero Entries')
    ax2.set_xticklabels(reshaped_data['variable'], rotation=0)
    ax2.set_xlabel('')
    ax2.set_ylim(bottom=0)  # Set the lower limit of the y-axis to 0

    plt.show()
    '''