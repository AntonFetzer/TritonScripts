import numpy as np
import matplotlib.pyplot as plt
# import sympy as sp
from Read.ReadSpenvis_tri import readSpenvis_tri
from Read.ReadGPSMacro import readGPSMacro

Expected = 'blue' # Blue
PlusColor = 'C1'  # Green
MinusColor = 'C2' # Orange
LEOColor = 'C8'   # yellow
GEOColor = 'C7'   # grey
VAPColor = 'C3'   # Red

'''
## Carrington Integral Electron Spectrum Parameters from EVT Analysis of Adnane Osmane
f0 = 10 ** 10   # cm-2 s-1 sr-1 
E0 = 0.13       # MeV
a = 3.7

Define symbolic variables and function
Energy = sp.Symbol('Energy')
f = 4 * sp.pi * f0 * (Energy / E0) ** (-a)
print('Funciton f:', f)

Calculate the derivative of the function with respect to Energy
fdiff = sp.diff(f, Energy)
print('Function fdiff:', fdiff)

Convert the symbolic functions to numerical functions
f = sp.lambdify(Energy, f, "numpy")
fdiff = sp.lambdify(Energy, fdiff, "numpy")

def f(Energy):
    return 4 * np.pi * f0 * (Energy / E0) ** (-a)

Energies = np.geomspace(0.13, 10, num=2)

print("Differential")
for E in Energies:
    print(f"{E:.4}", f"{-fdiff(E):.4}")
    print("/gps/hist/point", f"{E:.4}", f"{-fdiff(E):.4}")

print("Integral")
for E in Energies:
    print(f"{E:.4}", f"{f(E):.4}")
'''

## Plotting
plt.figure(1)

## Read in Carrington EVT spectrum
EVT_file = "/l/triton_work/Spectra/Carrington/Electron/CarringtonElectronINTEGRALPowTabelated.mac"
EVT_Data = readGPSMacro(EVT_file)
plt.plot(EVT_Data["Energy"], EVT_Data["Flux"], '-', label="Carrington peak Electron Flux", linewidth=2.5, color=Expected)

## FLUMIC spectrum
# FLUMIC Energy in MeV
FLUMIC_Energy = [2.0000E-01, 5.0000E-01, 8.0000E-01, 1.1000E+00, 1.4000E+00, 1.7000E+00, 2.0000E+00, 2.3000E+00, 2.6000E+00, 2.9000E+00, 3.2000E+00, 3.5000E+00, 3.8000E+00, 4.1000E+00, 4.4000E+00, 4.7000E+00, 5.0000E+00, 5.3000E+00, 5.6000E+00, 5.9000E+00]
# Flumic Flux in Flux (m-2 s-1 sr-1)
FLUMIC_Average_Flux = [1.0016E+10, 4.3943E+09, 1.9280E+09, 8.4589E+08, 3.7113E+08, 1.6283E+08, 7.1441E+07, 3.1344E+07, 1.3752E+07, 6.0337E+06, 2.6473E+06, 1.1615E+06, 5.0959E+05, 2.2358E+05, 9.8094E+04, 4.3038E+04, 1.8883E+04, 8.2847E+03, 3.6349E+03, 1.5948E+03]
# Flumic Max Flux in Flux (m-2 s-1 sr-1)
FLUMIC_Max_Flux = [14921000000, 6546400000, 2872200000, 1260200000, 552890000, 242580000, 106430000, 46695000, 20487000, 8988700, 3943800, 1730300, 759160, 333080, 146140, 64116, 28131, 12342, 5415.1, 2375.8]

# Convert from m-2 to cm-2 and from sr-1 to omnidirectional
FLUMIC_Average_Flux = np.array(FLUMIC_Average_Flux) * 1e-4 * 4 * np.pi
FLUMIC_Max_Flux = np.array(FLUMIC_Max_Flux) * 1e-4 * 4 * np.pi

plt.plot(FLUMIC_Energy, FLUMIC_Max_Flux, 'x-', label="FLUMIC GEO max Electron Flux", color='magenta', linewidth=1)

## Read in LEO A9 spectrum
LEO_file = "/l/triton_work/Spectra/Carrington/LEO/spenvis_tri.txt"
Protons, Electrons = readSpenvis_tri(LEO_file)
plt.plot(Electrons['Energy'], Electrons['Integral'], '+--', label="AE9 LEO mean Electron Flux", color=LEOColor)

## Read in Geostationary A9 spectrum
GEO_file = "/l/triton_work/Spectra/Carrington/GEO/spenvis_tri.txt"
Protons, Electrons = readSpenvis_tri(GEO_file)
plt.plot(Electrons['Energy'], Electrons['Integral'], '1-.', label="AE9 GEO mean Electron Flux", color=GEOColor)

## Read in Van-Allen Belt Probes A9 spectrum
VAP_file = "/l/triton_work/Spectra/Carrington/VAP/spenvis_tri.txt"
Protons, Electrons = readSpenvis_tri(VAP_file)
plt.plot(Electrons['Energy'], Electrons['Integral'], '2:', label="AE9 VAP mean Electron Flux", color=VAPColor)





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
plt.title("Electron Spectra Comparison")
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

plt.savefig("/l/triton_work/Spectra/Carrington/Electron/ElectronSpectrumComparison_NEW.pdf", format='pdf', bbox_inches="tight")
#plt.show()
