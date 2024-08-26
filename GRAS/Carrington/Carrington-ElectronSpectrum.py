import numpy as np
import matplotlib.pyplot as plt
# import sympy as sp
from GRAS.Read.ReadSpenvis_tri import readSpenvis_tri

Expected = 'blue' # Blue
PlusColor = 'C1'        # Green
MinusColor = 'C2'       # Orange
ISSColor = 'C8'         # yellow
GEOColor = 'C7'         # grey
VABColor = 'C3'         # Red

## Carrington Integral Electron Spectrum Parameters from EVT Analysis of Adnane Osmane
f0 = 10 ** 10   # cm-2 s-1 sr-1 
E0 = 0.13       # MeV
a = 3.7

# Define symbolic variables and function
# Energy = sp.Symbol('Energy')
# f = 4 * sp.pi * f0 * (Energy / E0) ** (-a)
# print('Funciton f:', f)

# Calculate the derivative of the function with respect to Energy
# fdiff = sp.diff(f, Energy)
# print('Function fdiff:', fdiff)

# Convert the symbolic functions to numerical functions
# f = sp.lambdify(Energy, f, "numpy")
# fdiff = sp.lambdify(Energy, fdiff, "numpy")

def f(Energy):
    return 4 * np.pi * f0 * (Energy / E0) ** (-a)

Energies = np.geomspace(0.13, 10, num=2)

# print("Differential")
# for E in Energies:
#     print(f"{E:.4}", f"{-fdiff(E):.4}")
    #print("/gps/hist/point", f"{E:.4}", f"{-fdiff(E):.4}")

print("Integral")
for E in Energies:
    print(f"{E:.4}", f"{f(E):.4}")


## Plotting
plt.figure(1)

plt.plot(Energies, f(Energies), '.-', label="Carrington Peak Electron Flux", linewidth=2.5, color=Expected)


## Read in ISS A9 spectrum
ISS_file = "/l/triton_work/Spectra/ISS/spenvis_tri.txt"

Protons, Electrons = readSpenvis_tri(ISS_file)

plt.plot(Electrons['Energy'], Electrons['Integral'], '.--', label="AE9 LEO Electron Flux", color=ISSColor)

## Read in Geostationary A9 spectrum
GEO_file = "/l/triton_work/Spectra/GEO/spenvis_tri.txt"

Protons, Electrons = readSpenvis_tri(GEO_file)

plt.plot(Electrons['Energy'], Electrons['Integral'], '.-.', label="AE9 GEO Electron Flux", color=GEOColor)

## Read in Van-Allen Belt Probes A9 spectrum
VAB_file = "/l/triton_work/Spectra/VAB/spenvis_tri.txt"

Protons, Electrons = readSpenvis_tri(VAB_file)

plt.plot(Electrons['Energy'], Electrons['Integral'], '.:', label="AE9 Van-Allen-Belt Electron Flux", color=VABColor)

""" 
## Electron fluxes from paper:
# "Proton, helium, and electron spectra during the large solar particle events of October–November 2003"
# doi:10.1029/2005JA011038
# Event 10/28/03 Table 7
# EPAM (0.04 to 0.32 MeV)
EPAMNorm = 6.75e8   # cm-2 sr-1 MeV-1
EPAMSlope = -1.90
# PET (1.6 – 8 MeV)
PETNorm = 1.46e8    # cm-2 sr-1 MeV-1
PETSlope = -4.27

# Convert to from sr-1 to omnidirectional
EPAMNorm *= 4 * np.pi  
PETNorm *= 4 * np.pi

# Convert from Fluence per day to Flux per second
EPAMNorm /= 24 * 3600
PETNorm /= 24 * 3600

# Energies
EPAMEnergy = np.geomspace(0.04, 0.32, num=10)
PETEnergy = np.geomspace(1.6, 8, num=10)

IntegralEPAMFlux = -EPAMNorm / (EPAMSlope + 1) * (EPAMEnergy ** (EPAMSlope + 1))
IntegralPETFlux = -PETNorm / (PETSlope + 1) * (PETEnergy ** (PETSlope + 1))

plt.plot(EPAMEnergy, IntegralEPAMFlux, label="10/28/03 EPAM Electron Flux")
plt.plot(PETEnergy, IntegralPETFlux, label="10/28/03 PET Electron Flux")
 """


## Plot formatting
plt.legend()
plt.yscale("log")
plt.xscale("log")
plt.title("Integral Electron Flux Comparison")
plt.xlabel("Electron Kinetic Energy [MeV]")
plt.ylabel("Integral Flux [cm-2 s-1]")
plt.xlim(0.1, 10)
plt.ylim(1e1, 1e12)
plt.minorticks_on()
plt.grid(axis='x', which='both')
plt.grid(axis='y', which='both')
# New code to create logarithmically spaced Y-axis ticks
y_ticks = np.logspace(1, 12, num=12)  # Generates 12 points between 10^1 and 10^12
plt.yticks(y_ticks)

x_ticks = [0.1, 0.2, 0.5, 1, 2, 5, 10]
plt.xticks(x_ticks, labels=[str(tick) for tick in x_ticks])

plt.savefig("/l/triton_work/Spectra/Carrington/Electron/ElectronSpectrumComparison.pdf", format='pdf', bbox_inches="tight")
#plt.show()
