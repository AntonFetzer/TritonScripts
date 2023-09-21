import pandas as pd
import re

# Read the file into a list of lines
with open('/l/triton_work/LunarRadiaitonAnalysis/LunarGCR-NarrowHist/0mm/Res/Protons_3842_594.csv', 'r') as file:
    lines = file.readlines()

# Initialize the data list
data = []

# Iterate through the lines to find the data
for line in lines:
    # Check if the line contains numeric data
    if re.match(r'\s*[0-9]', line):
        # Split the line into values and convert them to float or int, if possible
        values = line.split(',')
        row = []
        for value in values:
            try:
                row.append(float(value))
            except ValueError:
                try:
                    row.append(int(value))
                except ValueError:
                    row.append(value.strip())
        data.append(row)

# Create a pandas DataFrame from the data list
columns = ['lower', 'upper', 'mean', 'value', 'error', 'entries']
df = pd.DataFrame(data, columns=columns)

print(df.head())
