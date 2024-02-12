import os
import numpy as np
from GRAS.Dependencies.TotalKRadGrasSingle import totalkRadGrasSingle
import matplotlib.pyplot as plt
import pandas as pd
from natsort import natsorted
from GRAS.Dependencies.TotalKRadGras import totalkRadGras

Path = "/l/triton_work/RadEx/"

Folders = [f for f in os.listdir(Path) if f.endswith('mm')]
Folders = natsorted(Folders)

print(Folders)

# Initialize the combined DataFrames
combined_data_prot = pd.DataFrame(columns=['Dose', 'Error', 'Entries', 'Non zero entries'])
combined_data_elec = pd.DataFrame(columns=['Dose', 'Error', 'Entries', 'Non zero entries'])

for folder in Folders:
    temp_prot = totalkRadGrasSingle(os.path.join(Path, folder, "Res"), "Prot")
    temp_prot.index = [folder[5:]]  # Set the index of the temp DataFrame to the folder name without the "RadEx" prefix
    combined_data_prot = pd.concat([combined_data_prot, temp_prot])

    temp_elec = totalkRadGrasSingle(os.path.join(Path, folder, "Res"), "Elec")
    temp_elec.index = [folder[5:]]  # Set the index of the temp DataFrame to the folder name without the "RadEx" prefix
    combined_data_elec = pd.concat([combined_data_elec, temp_elec])

#print("Proton Data:")
#print(combined_data_prot)
#print("\nElectron Data:")
#print(combined_data_elec)

# Initialize the combined TID DataFrame
combined_data_tid = pd.DataFrame(columns=['Dose', 'Error'])

# Add the 'Dose' values of the proton and electron DataFrames
combined_data_tid['Dose'] = combined_data_prot['Dose'] + combined_data_elec['Dose']

# Use the error propagation rule for the 'Error' column
combined_data_tid['Error'] = np.sqrt(combined_data_prot['Error'] ** 2 + combined_data_elec['Error'] ** 2)

print("\nTID Data:")
print(combined_data_tid)

ThickList = np.array([int(folder[5:].split("mm")[0]) for folder in Folders])
print(ThickList)

# Import Shiedling Curve Data

Path2 = "/l/triton_work/ShieldingCurves"

Electrons = totalkRadGras(Path2 + "/AE9-500keV-10mm/Res/", "Elec")
Protons = totalkRadGras(Path2 + "/AP9-10MeV-10mm/Res/", "Prot")

x = np.linspace(0, 10, num=101, dtype=float)

# Dose with Error plot
fig, ax = plt.subplots(figsize=(6, 5))  # Adjust figure size

####### Plot 10kRad line #########
CriticalDose = [1] * 101
ax.plot(x, CriticalDose, color='k', linewidth=2, label='1 krad per month')

# ----------------------------------------- Plot electrons500kev ---------------------------------------------------------
ax.errorbar(x, Electrons[0], Electrons[1], fmt='C0 ', capsize=4, label="Electron Dose Planar Shielding", elinewidth=1, capthick=1)

# ----------------------------------------- Plot protons10MeV -------------------------------------------------------
ax.errorbar(x, Protons[0], Protons[1], fmt='C1 ', capsize=4, label="Proton Dose Planar Shielding", elinewidth=1, capthick=1)

Total = Electrons
Total[0] = Electrons[0] + Protons[0]
Total[1] = np.sqrt(Electrons[1] * Electrons[1] + Protons[1] * Protons[1])

# ----------------------------------------- Plot Total -------------------------------------------------------
ax.errorbar(x, Total[0], Total[1], fmt='C2 ', capsize=4, label="Total Ionizing Dose Planar Shielding", elinewidth=1, capthick=1)

ax.errorbar(ThickList, combined_data_elec['Dose'], yerr=combined_data_elec['Error'], capsize=10, capthick=2, elinewidth=2, fmt='sC0', markersize=5, label='RadEx Electrons')
ax.errorbar(ThickList, combined_data_prot['Dose'], yerr=combined_data_prot['Error'], capsize=10, capthick=2, elinewidth=2, fmt='oC1', markersize=5, label='RadEx Protons')
ax.errorbar(ThickList, combined_data_tid['Dose'], yerr=combined_data_tid['Error'], capsize=10, capthick=2, elinewidth=2, fmt='DC2', markersize=5, label='RadEx Total Ionising Dose')

ax.set_ylabel('Ionising Dose [krad per month]')
ax.set_yscale('log')
ax.set_title('Ionising Dose behind shielding cutouts in RadEx')
ax.set_xticks(ThickList)
ax.set_xticklabels(combined_data_prot.index, rotation=0)
ax.set_xlabel('Aluminium shielding thickness')
ax.legend(loc='upper right')
ax.set_xlim(-0.2, 8)
ax.set_ylim(1e-1, 1e3)

ax.grid(True, which='both', linestyle='--', alpha=0.5)  # Add gridlines

plt.show()
#plt.savefig(os.path.join(Path, "Plots", "RadEx-Total.pdf"), format='pdf', bbox_inches="tight")





'''
# Entries and Non Zero Entries plot
ax2 = combined_data[['Entries', 'Non zero entries']].plot(kind='bar', legend=True)
ax2.set_ylabel('Count')
ax2.set_title('Entries and Non Zero Entries')
ax2.set_xticklabels(combined_data.index, rotation=45)
ax2.set_xlabel('')
ax2.set_ylim(bottom=0)  # Set the lower limit of the y-axis to 0
'''

#plt.show()


