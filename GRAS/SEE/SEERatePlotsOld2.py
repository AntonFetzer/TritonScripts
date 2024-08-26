import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd

plt.rc('axes', axisbelow=True)

# Load the CSV file into a dictionary with numerical data as 1D numpy arrays
file_path = '/l/triton_work/LET_Histograms/Carrington/SEERatesCorrectable.csv'

data = {}
with open(file_path, 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        for key, value in row.items():
            if key not in data:
                data[key] = []
            data[key].append(value)

# Convert relevant columns to numpy arrays of floats
for key in ['SEE_Rate', 'SEE_Error', 'EntriesContributingToSEE']:
    data[key] = np.array(data[key], dtype=float)

# Shielding collumn is integer
data['Shielding'] = np.array(data['Shielding'], dtype=int)

# Convert other non-numeric data to numpy arrays of strings
for key in ['Data', 'Crossection']:
    data[key] = np.array(data[key])

print(data)

DataNames = [
    "Carrington SEP EVT",
    "Carrington SEP +2 Sigma",
    "Carrington SEP -2 Sigma",
    # "Carrington Electron",
    # "ISS AP9",
    # "ISS Solar Proton",
    # "VAB AP9",
    # "VAB AE9",
    # "VAB Solar Proton",
    # "GEO AP9",
    # "GEO Solar Proton",
    # "GEO Solar Proton 5min Peak Flux",
]

Colours = ['C1', 'C0', 'C2', 'C8', 'C3', 'C9', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']

# Extract the crossection name
CrossectionName = data['Crossection'][0]

# Check that the crossection name is the same for all entries
if not all(data['Crossection'] == CrossectionName):
    raise ValueError("Crossection name is not the same for all entries")

# Plotting the data
plt.figure(0)

# Convert integer shielding values to string labels for categorical plotting
data['Shielding_str'] = np.array([f'{val}mm' for val in data['Shielding']])

for i, DataName in enumerate(DataNames):
    mask = np.array(data['Data']) == DataName
    plt.errorbar(data['Shielding_str'][mask], data['SEE_Rate'][mask], yerr=data['SEE_Error'][mask], label=DataName, color=Colours[i])

plt.axhline(1/8e+3, linestyle='--', label="1 Bitflip per kB per Second", color='black')
plt.axhline(1/8e+6, linestyle=':', label="1 Bitflip per MByte per second", color='black')
plt.axhline(1/8e+9, linestyle='-.', label="1 Bitflip per GByte per second", color='black')
plt.axhline(1/8e+12, linestyle='--', label="1 Bitflip per TByte per second", color='black')

# Set the x-axis to integer values representing the shielding thickness
plt.xticks([0, 1, 2, 4, 8, 16], ['0mm', '1mm', '2mm', '4mm', '8mm', '16mm'])


# plt.ylim(1e-15, 1e-3)
plt.yscale("log")
plt.grid()
plt.title(CrossectionName)

plt.xlabel("Aluminium Shielding Thickness")
plt.ylabel("Single Bit Upset Rate [s-1 bit-1]")
plt.legend(loc='lower left')

plt.savefig("/l/triton_work/LET_Histograms/Carrington/" + CrossectionName + " Rates.pdf", format='pdf', bbox_inches="tight")
# plt.show()