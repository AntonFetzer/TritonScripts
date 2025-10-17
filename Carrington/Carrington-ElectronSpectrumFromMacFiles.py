import numpy as np
import matplotlib.pyplot as plt

from Read.ReadGPSMacro import readGPSMacro

Datasets = ['EVT', 'Plus', 'Minus', 'LEO', 'MEO', 'VAP', 'GEO']

# Main dict
S = {}

# initialise dataset ditcts
for d in Datasets:
    S[d] = {}

# Define colors
S["EVT"]["color"] = 'blue'   # Blue
S["Plus"]["color"] = 'C1'    # Green
S["Minus"]["color"] = 'C2'   # Orange
S["LEO"]["color"] = 'C8'     # yellow
S["MEO"]["color"] = 'C9'     # turquoise
S["VAP"]["color"] = 'C3'     # Red
S["GEO"]["color"] = 'C7'     # grey

# File paths
S["EVT"]["file"] = "/l/triton_work/Spectra/Carrington/Electron/CarringtonElectronDiffPowTabelated.mac"
S["LEO"]["file"] = "/l/triton_work/Spectra/Carrington/LEO/LEO-electron.mac"
S["GEO"]["file"] = "/l/triton_work/Spectra/Carrington/GEO/GEO-electron.mac"
S["VAP"]["file"] = "/l/triton_work/Spectra/Carrington/VAP/VAP-electron.mac"
S["MEO"]["file"] = "/l/triton_work/Spectra/Carrington/MEO/MEO-electron.mac"

# Read in the data
for d in Datasets:
    if "file" in S[d]:
        Data = readGPSMacro(S[d]["file"])
        S[d]["Energy"] = Data['Energy']
        S[d]["Flux"] = Data['Flux']

# Plotting
plt.figure(1)

for d in Datasets:
    if "Energy" in S[d]:
        plt.plot(S[d]["Energy"], S[d]["Flux"], '.-', label=d, linewidth=2.5, color=S[d]["color"])

plt.yscale("log")
plt.xscale("log")
plt.grid(which="both")
plt.legend()
plt.title("Carrington Electron Spectra")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Integral Flux [cm$^{-2}$ s$^{-1}$ sr$^{-1}$]")

plt.show()