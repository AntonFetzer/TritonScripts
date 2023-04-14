import pandas as pd
import os

file_name = "/home/anton/Desktop/triton_work/Permutations/5Layer/Analysis/5Layer-Raw.csv"

# Read the CSV file and create a dataframe
dataframe = pd.read_csv(file_name)

#dataframe['Electron Rank'] = 0
#dataframe['Proton Rank'] = 0
#dataframe['Total Rank'] = 0

# Calculate the rank based on "Total Dose [krad/Month]" and store it in the "Total Rank" column
dataframe['Electron Rank'] = dataframe['Electron Dose [krad/Month]'].rank(method='min').astype(int)
dataframe['Proton Rank'] = dataframe['Proton Dose [krad/Month]'].rank(method='min').astype(int)
dataframe['Total Rank'] = dataframe['Total Dose [krad/Month]'].rank(method='min').astype(int)

'''
# Add a new column "Total Rel Err [%]" containing the relative error calculation
dataframe['Total Rel Err [%]'] = (dataframe['Total Err [krad/Month]'] * 100 / dataframe['Total Dose [krad/Month]'])

# Add a new column with a boolean value depending on which error is larger
dataframe['Proton Err Larger'] = dataframe['Proton Err [krad/Month]'] > dataframe['Electron Err [krad/Month]']

print("Number of True entries in 'Proton Err Larger':", dataframe['Proton Err Larger'].sum())

# Sort the dataframe by 'Total Rank' in ascending order
dataframe.sort_values(by='Total Rank', inplace=True)

# Calculate the difference between consecutive rows for total dose and total error
dataframe['Total Dose Diff'] = dataframe['Total Dose [krad/Month]'].diff().shift(-1)

# Add a new column with a boolean value based on the comparison of the total error difference and total dose difference
dataframe['Total Err Diff Larger'] = dataframe['Total Err [krad/Month]'] > dataframe['Total Dose Diff']

# After all the columns have been filled, sort the DataFrame by the 'Combination #' column to restore the original order of rows
dataframe.sort_values(by='Combination #', inplace=True)
'''

# Replace "-Raw" with "-Processed" in the file name
new_file_name = file_name.replace("-Raw", "-Processed")

# Export the DataFrame to a CSV file
dataframe.to_csv(new_file_name, index=False)

print(dataframe)
