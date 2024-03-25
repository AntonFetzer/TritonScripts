import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from GRAS.Read.ReadSpenvis_tri import readSpenvis_tri

## Carrington Integral Electron Spectrum Parameters from EVT Analysis of Adnane Osmane
f0 = 10 ** 10   # cm-2 s-1 sr-1 
E0 = 0.13       # MeV
a = 3.7

# Define symbolic variables and function
Energy = sp.Symbol('Energy')
f = 4 * sp.pi * f0 * (Energy / E0) ** (-a)
print('Funciton f:', f)

# Calculate the derivative of the function with respect to Energy
# fdiff = sp.diff(f, Energy)
# print('Function fdiff:', fdiff)

# Convert the symbolic functions to numerical functions
f = sp.lambdify(Energy, f, "numpy")
# fdiff = sp.lambdify(Energy, fdiff, "numpy")

Evals = np.geomspace(0.13, 100, num=10)

# print("Differential")
# for E in Evals:
#     print(f"{E:.4}", f"{-fdiff(E):.4}")
    #print("/gps/hist/point", f"{E:.4}", f"{-fdiff(E):.4}")

# print("Integral")
# for E in Evals:
#     print(f"{E:.4}", f"{f(E):.4}")


## Plotting
plt.figure(1)

plt.plot(Evals, f(Evals), label="Carrington Peak Electron Flux", linewidth=2.5)
# plt.plot(Evals, -fdiff(Evals), label="Carrington Differential Electron Flux")

## Read in ISS A9 spectrum
ISS_file = "/l/triton_work/Spectra/ISS/spenvis_tri.txt"

Protons, Electrons = readSpenvis_tri(ISS_file)

#plt.plot(Protons[0], Protons[1], label="ISS AP9 Integral Protons")
plt.plot(Electrons[0], Electrons[1], '.-', label="AE9 LEO Electron Flux")

## Read in Geostationary A9 spectrum
GEO_file = "/l/triton_work/Spectra/GEO/spenvis_tri.txt"

Protons, Electrons = readSpenvis_tri(GEO_file)

#plt.plot(Protons[0], Protons[1], label="GEO AP9 Integral Protons")
plt.plot(Electrons[0], Electrons[1], '.-.', label="AE9 GEO Electron Flux")

## Read in Van-Allen Belt Probes A9 spectrum
VAB_file = "/l/triton_work/Spectra/Van-Allen-Belt-Probes/spenvis_tri.txt"

Protons, Electrons = readSpenvis_tri(VAB_file)

#plt.plot(Protons[0], Protons[1], label="VAB AP9 Integral Protons")
plt.plot(Electrons[0], Electrons[1], '.:', label="AE9 Van-Allen-Belt Electron Flux")

## Plot formatting
plt.legend()
plt.yscale("log")
plt.xscale("log")
plt.title("Integral Electron Flux Comparison")
plt.xlabel("Electron Kinetic Energy [MeV]")
plt.ylabel("Integral Flux [cm-2 s-1]")
plt.xlim(0.1, 20)
plt.ylim(1e1, 1e12)
plt.minorticks_on()
plt.grid(axis='x', which='both')
plt.grid(axis='y', which='both')
# New code to create logarithmically spaced Y-axis ticks
y_ticks = np.logspace(1, 12, num=12)  # Generates 12 points between 10^1 and 10^12
plt.yticks(y_ticks)

plt.savefig("/l/triton_work/Spectra/Carrington/Electron/SpectrumComparison.pdf", format='pdf', bbox_inches="tight")
#plt.show()
