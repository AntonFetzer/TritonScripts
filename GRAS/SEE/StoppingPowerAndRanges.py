import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

File = '/l/triton_work/LET/Stopping Powers and Ranges/StoppingPower-Ranges.csv'

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

# Convert LET from MeV cm2/g to MeV cm2/mg
data['Electron LET'] = data['Electron LET'] / 1000
data['Proton LET'] = data['Proton LET'] / 1000
G4Electron_LET = G4Electron_LET / 1000
G4Proton_LET = G4Proton_LET / 1000

# Create a new figure with 2 subplots
fig, (ax1, ax2) = plt.subplots(2)

xlims = [1e0, 1e3]

# Plot the LET data with a linear y-axis
ax1.plot(data['Electron Energy'], data['Electron LET'], label='ESTAR Electron LET')
ax1.plot(data['Proton Energy'], data['Proton LET'], label='PSTAR Proton LET')
ax1.plot(G4Electron_Energy, G4Electron_LET, label='G4Electron LET', linestyle='dotted')
ax1.plot(G4Proton_Energy, G4Proton_LET, label='G4Proton LET', linestyle='dotted')

# Formatting the first subplot
ax1.set_xlabel('Energy (MeV)')
ax1.set_ylabel('LET (MeV cm2/mg)')
ax1.grid(True)
ax1.legend()
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(xlims)


# Convert Range from g/cm2 to mm of Aluminium
data['Electron Range'] = data['Electron Range'] / 0.27
data['Proton Range'] = data['Proton Range'] / 0.27

# Plot the Range data on the second subplot
ax2.plot(data['Electron Energy'], data['Electron Range'], label='ESTAR Electron Range', linestyle='dashed')
ax2.plot(data['Proton Energy'], data['Proton Range'], label='PSTAR Proton Range', linestyle='dashed')

# Formatting the second subplot
ax2.set_xlabel('Energy (MeV)')
ax2.set_ylabel('Range (mm of Aluminium)')  # Update the y-label
ax2.grid(True)
ax2.legend()
ax2.set_xscale('log')
#ax2.set_yscale('log')
ax2.set_xlim(xlims)
ax2.set_ylim(0, 6)

# Adjust the space between the subplots
plt.tight_layout()

# Save the plot
plt.savefig("/l/triton_work/LET/Stopping Powers and Ranges/Stopping Power and Ranges in Silicon.pdf", format='pdf', bbox_inches="tight")

# Print maximum LET and corresponding energy of PSTAR data
max_LET = data['Proton LET'].max()
max_LET_energy = data['Proton Energy'][data['Proton LET'].idxmax()]
RangeAtMaxLET = data['Proton Range'][data['Proton LET'].idxmax()]
print('Maximum LET of Proton data:', max_LET, 'MeV cm2/mg')
print('Energy of maximum LET of Proton data:', max_LET_energy, 'MeV')
print('Range at maximum LET of Proton data:', RangeAtMaxLET, 'mm of Aluminium')