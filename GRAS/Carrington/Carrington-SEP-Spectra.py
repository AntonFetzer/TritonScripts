import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from GRAS.Read.ReadSpenvis_tri import readSpenvis_tri
from GRAS.Read.ReadGPSMacro import readGPSMacro
from GRAS.Read.ReadSpenvis_sef import readSpenvis_sef

## Raw Data from Adnanes E-mail on 13/02/2023 15:53
Energies = np.array([10, 30, 60, 100, 200])      # MeV

ExpectedInt = np.array([4.6e4, 8.5e3, 2.3e3, 1.1e3, 2.9e2])     # cm-2 s-1 str-1
Minus = np.array([2.2e4, 5.6e3, 1.4e3, 4.3e2, 1.1e2])       # cm-2 s-1 str-1
Plus = np.array([1.3e5, 1.6e4, 5.4e3, 4.4e3, 1.4e3])        # cm-2 s-1 str-1

# Convert from str-1 to omnidirectional
ExpectedInt *= 4 * np.pi
Minus *= 4 * np.pi
Plus *= 4 * np.pi

Colours = ['C1', 'C0', 'C2', 'C7']

plt.figure(1)

plt.fill_between(Energies, Plus, ExpectedInt, color='C1', alpha=0.5)
plt.fill_between(Energies, ExpectedInt, Minus, color='C2', alpha=0.5)

plt.plot(Energies, Plus, label="Carrington SEP EVT +2 Sigma", color='C1')
plt.plot(Energies, ExpectedInt, label="Carrington SEP EVT Expected", color='C0', linewidth=4)
plt.plot(Energies, Minus, label="Carrington SEP EVT -2 Sigma", color='C2')









## Numbers from paper "The magnitude and effects of extreme solar particle events" DOI: 10.1051/swsc/2014017 Figue 9
# SPE2003 = np.array([1.5678e+10, 3.4589e+9, 7.0664e+8, 1.5325e+8, 1.9409e+7])    # cm-2 Event Fluence

## Numbers from Townsend Carrington paper 
# TownsendCarrington = np.array([1e11, 2e10, 4e9, 9e8, 6e7])


# plt.plot(Energies, TownsendCarrington, ':', label="Townsend Carrington Estimate", color='C9', linewidth=4)
# plt.plot(Energies, SPE2003, '--', label="SPE 2003", color='C8', linewidth=4)





### A9 Spectra ###

## Read in ISS A9 spectrum
ISS_file = "/l/triton_work/Spectra/ISS/spenvis_tri.txt"

Protons, Electrons = readSpenvis_tri(ISS_file)

plt.plot(Protons['Energy'], Protons['Integral'], '.-', label="AE9 LEO Proton Flux")

## Read in Geostationary A9 spectrum
GEO_file = "/l/triton_work/Spectra/GEO/spenvis_tri.txt"

Protons, Electrons = readSpenvis_tri(GEO_file)

plt.plot(Protons['Energy'], Protons['Integral'], '.-.', label="AE9 GEO Proton Flux")

## Read in Van-Allen Belt Probes A9 spectrum
VAB_file = "/l/triton_work/Spectra/Van-Allen-Belt-Probes/spenvis_tri.txt"

Protons, Electrons = readSpenvis_tri(VAB_file)

plt.plot(Protons['Energy'], Protons['Integral'], '.:', label="AE9 Van-Allen-Belt Proton Flux")


## SAPPHIRE Solar Proton Spectra

ISS



## Read in GPS Macro
# ExpectedInt = "/l/triton_work/Spectra/Carrington/SEP-Final/Carrington-SEP-Expected-Int.mac"
# ExpectedIntWith0 = "/l/triton_work/Spectra/Carrington/SEP-Final/Carrington-SEP-Expected-Int-With0.mac"
# Minus2SigmaIntWith0 = "/l/triton_work/Spectra/Carrington/SEP-Final/Carrington-SEP-Minus2Sigma-Int-With0.mac"
# Plus2SigmaIntWith0 = "/l/triton_work/Spectra/Carrington/SEP-Final/Carrington-SEP-Plus2Sigma-Int-With0.mac"

# ExpectedInt = readGPSMacro(ExpectedInt)
# ExpectedIntWith0 = readGPSMacro(ExpectedIntWith0)
# Minus2SigmaIntWith0 = readGPSMacro(Minus2SigmaIntWith0)
# Plus2SigmaIntWith0 = readGPSMacro(Plus2SigmaIntWith0)

# plt.plot(ExpectedInt['Energy'], ExpectedInt['Flux'], label="Carrington SEP EVT Expected", color='r', linewidth=4)
# plt.plot(ExpectedIntWith0['Energy'], ExpectedIntWith0['Flux'], label="Carrington SEP EVT Expected with 0", color='r', linewidth=4, linestyle='--')
# plt.plot(Minus2SigmaIntWith0['Energy'], Minus2SigmaIntWith0['Flux'], label="Carrington SEP EVT -2 Sigma with 0", color='r', linestyle='--')
# plt.plot(Plus2SigmaIntWith0['Energy'], Plus2SigmaIntWith0['Flux'], label="Carrington SEP EVT +2 Sigma with 0", color='r', linestyle='--')


# plt.xlim(9, 250)
# plt.ylim(1e2, 2e7)
plt.yscale("log")
plt.xscale("log")
plt.title("Carrington Solar Energetic Proton Comparison")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Integral Flux [cm-2 s-1]")
plt.legend()
plt.grid(which='both')

# plt.savefig("/l/TritonPlots/Carrington/SpectrumComparison.pdf", format='pdf', bbox_inches="tight")
plt.show()