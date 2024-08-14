import matplotlib.pyplot as plt
import numpy as np
from GRAS.Read.ReadLETHistos import readLETHistos
from GRAS.Dependencies.TotalLETHistos import totalLETHistos

File = '/l/triton_work/LET_Histograms/Stopping Powers and Ranges/StoppingPower-Ranges.csv'

# Column names
columns = ['Proton Energy', 'Proton LET', 'Proton Range', 'Electron Energy', 'Electron LET', 'Electron Range']

# Initialize the dictionary to store the data
data = {col: [] for col in columns}

# Read the CSV file manually and populate the dictionary
with open(File, 'r') as f:
    for _ in range(5):  # Skip the first 5 rows
        next(f)
    
    for line in f:
        values = line.strip().split(',')
        for i, value in enumerate(values):
            if value.strip():  # Check if the value is not an empty string
                data[columns[i]].append(float(value))
            else:
                data[columns[i]].append(np.nan)  # Append NaN for missing values

# Convert lists to numpy arrays
for key in data:
    data[key] = np.array(data[key])

# Convert LET from MeV cm2/g to MeV cm2/mg
data['Electron LET'] = data['Electron LET'] / 1000
data['Proton LET'] = data['Proton LET'] / 1000


########## Read in GRAS data ##########
ProtonPath = '/l/triton_work/LET_Histograms/Mono/Protons/'
ElectronPath = '/l/triton_work/LET_Histograms/Mono/Electrons/'

Proton_LET = {}
Proton_LET['Energies'] = np.zeros(11, dtype=int)
Proton_LET['Mean'] = np.zeros(11)

Electron_LET = {}
Electron_LET['Energies'] = np.zeros(11, dtype=int)
Electron_LET['Mean'] = np.zeros(11)

for i in range(11):
    Proton_LET['Energies'][i] = 2**i
    LET, _ = totalLETHistos(ProtonPath + f"{Proton_LET['Energies'][i]}MeV/Res/")
    # Calculate average mean LET
    AverageMeanLET = np.sum(LET['mean'] * LET['entries']) / np.sum(LET['entries'])
    #print(f"Average Mean LET for {Proton_LET['Energies'][i]} MeV Protons: {AverageMeanLET:.3f} MeV cm2/mg")
    Proton_LET['Mean'][i] = AverageMeanLET

for i in range(11):
    Electron_LET['Energies'][i] = 2**i
    LET, _ = totalLETHistos(ElectronPath + f"{Electron_LET['Energies'][i]}MeV/Res/")
    # Calculate average mean LET
    AverageMeanLET = np.sum(LET['mean'] * LET['entries']) / np.sum(LET['entries'])
    #print(f"Average Mean LET for {Electron_LET['Energies'][i]} MeV Electrons: {AverageMeanLET:.3f} MeV cm2/mg")
    Electron_LET['Mean'][i] = AverageMeanLET


 
# Create a new figure with 2 subplots
fig, (ax1, ax2) = plt.subplots(2)

xlims = [1e1, 1e4]

# Plot the LET data
ax1.plot(data['Electron Energy'], data['Electron LET'], label='ESTAR Electron LET')
ax1.plot(data['Proton Energy'], data['Proton LET'], label='PSTAR Proton LET')

# Add GRAS data to the plot
ax1.plot(Proton_LET['Energies'], Proton_LET['Mean'], '+', label='GRAS Proton LET')
ax1.plot(Electron_LET['Energies'], Electron_LET['Mean'], '+', label='GRAS Electron LET')

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
# ax2.set_yscale('log')
ax2.set_xlim(xlims)
ax2.set_ylim(0, 10)

# Adjust the space between the subplots
plt.tight_layout()

# Save the plot
plt.savefig("/l/triton_work/LET_Histograms/Stopping Powers and Ranges/Stopping Power and Ranges in Silicon.pdf", format='pdf', bbox_inches="tight")
# plt.show()

# Print maximum LET and corresponding energy of PSTAR data
max_LET_index = np.nanargmax(data['Proton LET'])  # Use np.nanargmax to ignore NaNs
max_LET = data['Proton LET'][max_LET_index]
max_LET_energy = data['Proton Energy'][max_LET_index]
RangeAtMaxLET = data['Proton Range'][max_LET_index]
print('Maximum LET of Proton data:', max_LET, 'MeV cm2/mg')
print('Energy of maximum LET of Proton data:', max_LET_energy, 'MeV')
print('Range at maximum LET of Proton data:', RangeAtMaxLET*1000, 'µm of Aluminium')
