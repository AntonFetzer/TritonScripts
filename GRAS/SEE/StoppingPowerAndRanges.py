import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

File = '/home/anton/Desktop/triton_work/Chess1 GNSS SEE Analysis/Stopping Powers and Ranges/StoppingPower-Ranges.csv'

# Column names
columns = ['Proton Energy', 'Proton LET', 'Proton Range', 'Electron Energy', 'Electron LET', 'Electron Range']

# Read the CSV file into a pandas DataFrame
data = pd.read_csv(File, skiprows=5, names=columns)

C = 23.30  # to convert from MeV/cm to Mev cm2 g-1 in silicon with silicion density being 2.33 g/cm3
# Additional data
G4Electron_Energy = np.array([0.001, 0.01, 0.1, 1, 2.5, 5, 10, 20, 40, 80, 160, 320, 640, 1000, 10000])
G4Electron_LET = np.array([2000, 401.81, 77.47, 35.79, 37.293, 40.854, 46.85, 57.539, 78.11, 119.26, 202.57, 370.71, 708.86, 1090.7, 10874])

G4Electron_LET = G4Electron_LET / C

G4Proton_Energy = np.array([0.01, 0.1, 1, 2, 5, 10, 25, 50, 100, 150, 200, 500, 1000, 10000, 100000, 1e6, 1e7])
G4Proton_LET = np.array([7786.6, 11782, 4088, 2602.7, 1359.4, 805.05, 392.61, 229.38, 136, 102.14, 84.543, 52.008, 42.068, 42.592, 51.965, 64.612, 137.44])

G4Proton_LET = G4Proton_LET / C

fig, ax1 = plt.subplots()

# Plot the LET data
ax1.loglog(data['Electron Energy'], data['Electron LET'], label='ESTAR Electron LET')
ax1.loglog(data['Proton Energy'], data['Proton LET'], label='PSTAR Proton LET')
ax1.loglog(G4Electron_Energy, G4Electron_LET, label='G4Electron LET', linestyle='dotted')
ax1.loglog(G4Proton_Energy, G4Proton_LET, label='G4Proton LET', linestyle='dotted')

# Formatting the first y-axis
ax1.set_xlabel('Energy (MeV)')
ax1.set_ylabel('LET (MeV cm2/g)')
ax1.grid(True)

# Create a second y-axis
ax2 = ax1.twinx()

# Plot the Range data on the second y-axis
ax2.loglog(data['Electron Energy'], data['Electron Range'], label='ESTAR Electron Range', linestyle='dashed')
ax2.loglog(data['Proton Energy'], data['Proton Range'], label='PSTAR Proton Range', linestyle='dashed')

# Formatting the second y-axis
ax2.set_ylabel('Range (g/cm2)')

# Adding the legend
fig.legend(loc="upper left", bbox_to_anchor=[0, 1], bbox_transform=ax1.transAxes)

# Save the plot
plt.savefig("/home/anton/Desktop/triton_work/Chess1 GNSS SEE Analysis/Stopping Powers and Ranges/Stopping Power and Ranges in Silicon.pdf", format='pdf', bbox_inches="tight")
