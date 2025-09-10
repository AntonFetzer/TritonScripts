import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rc('axes', axisbelow=True)

df = pd.read_csv('/l/triton_work/LET_Histograms/Carrington/SEERatesCorrectable.csv')

DataNames = {
    "Carrington SEP EVT": "Carrington SEP EVT",
    "Carrington SEP +2 Sigma": "Carrington SEP +2 Sigma",
    "Carrington SEP -2 Sigma": "Carrington SEP -2 Sigma",
    # "Carrington Electron", # Electons dont seem to cause SEEs
    "ISS AP9": "ISS LEO Trapped Protons",
    # "ISS Solar Proton",  # Much lower than any other data
    "VAB AP9": "Van-Allen-Belt Probes Trapped Protons",
    # "VAB AE9", # Electons dont seem to cause SEEs
    "VAB Solar Proton": "VAB Solar Proton",
    # "GEO AP9", # On GEO trapped protons are not a problem
    "GEO Solar Proton": "GEO Solar Proton",
    "GEO Solar Proton 5min Peak Flux": "GEO Solar Proton 5min Peak Flux"
}

Colours = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11']

# print(df.to_string())

# print(df.keys())

# Extract the crossection name
CrossectionName = df['Crossection'][0]

# Check that the crossection name is the same for all entries
if not all(df['Crossection'] == CrossectionName):
    raise ValueError("Crossection name is not the same for all entries")

# Plot +2 Sigma and -2 Sigma as fill between
plt.fill_between(df.Shielding[df.Data == 'Carrington SEP EVT'], df.SEE_Rate[df.Data == 'Carrington SEP EVT'], 
                 df.SEE_Rate[df.Data == 'Carrington SEP +2 Sigma'], color='C1', alpha=0.5, label="Carrington SEP +2 Sigma")
plt.fill_between(df.Shielding[df.Data == 'Carrington SEP EVT'], df.SEE_Rate[df.Data == 'Carrington SEP EVT'], 
                 df.SEE_Rate[df.Data == 'Carrington SEP -2 Sigma'], color='C2', alpha=0.5, label="Carrington SEP -2 Sigma")

# Create a new list without +-2 Sigma
DataNames = {
    "Carrington SEP EVT": "Carrington SEP EVT",
    # "Carrington Electron", # Electons dont seem to cause SEEs
    "ISS AP9": "ISS LEO Trapped Protons",
    # "ISS Solar Proton",  # Much lower than any other data
    "VAB AP9": "Van-Allen-Belt Probes Trapped Protons",
    # "VAB AE9", # Electons dont seem to cause SEEs
    "VAB Solar Proton": "VAB Solar Proton",
    # "GEO AP9", # On GEO trapped protons are not a problem
    "GEO Solar Proton": "GEO Solar Proton",
    "GEO Solar Proton 5min Peak Flux": "GEO Solar Proton 5min Peak Flux"
}

# Plot the error bars for the remaining data
for i, DataName in enumerate(DataNames):
    # Select the data with the correct name
    Cd = df[df['Data'] == DataName]
    plt.errorbar(Cd.Shielding, Cd.SEE_Rate, yerr=Cd.SEE_Error, capsize=5, elinewidth=1, capthick=2, label=DataNames[DataName], color=Colours[i])

# plt.axhline(1/8e+3, linestyle='--', label="1 Bitflip per kB per Second", color='black')
# plt.axhline(1/8e+6, linestyle=':', label="1 Bitflip per MByte per second", color='black')
plt.axhline(1/8e+9, linestyle='-.', label="1 Bitflip per GByte per second", color='black')
plt.axhline(1/8e+12, linestyle='--', label="1 Bitflip per TByte per second", color='black')

#plt.ylim(1e-15, 1e-3)
plt.yscale("log")
plt.grid()
plt.title(CrossectionName)

plt.xlabel("Aluminium Shielding Thickness")
plt.ylabel("Single Bit Upset Rate [s-1 bit-1]")

# Manually add legend entries in the desired order
handles, labels = plt.gca().get_legend_handles_labels()

# Create a new order for the legend
order = [
    labels.index("1 Bitflip per GByte per second"),
    labels.index("1 Bitflip per TByte per second"),
    labels.index("Carrington SEP +2 Sigma"), 
    labels.index("Carrington SEP EVT"), 
    labels.index("Carrington SEP -2 Sigma"),
    labels.index("ISS LEO Trapped Protons"),
    labels.index("Van-Allen-Belt Probes Trapped Protons"),
    labels.index("VAB Solar Proton"),
    labels.index("GEO Solar Proton"),
    labels.index("GEO Solar Proton 5min Peak Flux")
    ]

plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='lower left')


plt.savefig("/l/triton_work/LET_Histograms/Carrington/" + CrossectionName + " Rates.pdf", format='pdf', bbox_inches="tight")
# plt.show()