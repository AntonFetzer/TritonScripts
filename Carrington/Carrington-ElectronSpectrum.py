import sys, os
# Ensure the package root (the "Python" folder) is on sys.path so "from Read.* import ..." works
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from Read.ReadSpenvis_tri import readSpenvis_tri
from Read.ReadGPSMacro import readGPSMacro

Path = "/scratch/work/fetzera1/"

Expected = 'blue' # Blue
# PlusColor = 'C1'  # Green
# MinusColor = 'C2' # Orange
LEOColor = 'C8'   # yellow
GEOColor = 'C7'   # grey
# VAPColor = 'C3'   # Red
MEOColor = 'C3'   # Red
FlumicColor = 'magenta' # Magenta
GOESColor = 'cyan' # Cyan

'''  Deprecated 
## Carrington Integral Electron Spectrum Parameters from EVT Analysis of Adnane Osmane
f0 = 10 ** 10   # cm-2 s-1 sr-1 
E0 = 0.13       # MeV
a = 3.7

# Define symbolic variables and function
Energy = sp.Symbol('Energy')
f = 4 * sp.pi * f0 * (Energy / E0) ** (-a)
print('Funciton f:', f)

# Calculate the derivative of the function with respect to Energy
fdiff = sp.diff(f, Energy)
print('Function fdiff:', fdiff)

# Convert the symbolic functions to numerical functions
f = sp.lambdify(Energy, f, "numpy")
fdiff = sp.lambdify(Energy, fdiff, "numpy")

def f(Energy):
    return 4 * np.pi * f0 * (Energy / E0) ** (-a)
'''



## Carrington differential electron return levels (1-in-150 year), Mikko Borg, 28 October 2026
E1, E2 = 0.13, 1.2           # MeV
j1, j2 = 1.07e8, 2.11e3      # keV-1 cm-2 s-1 sr-1
j1, j2 = j1 * 1e3, j2 * 1e3  # MeV-1 cm-2 s-1 sr-1
Emax = 1000                  # MeV, upper end of the production macro
g = np.log(j1 / j2) / np.log(E2 / E1)   # differential slope, 4.874551

print("Differential slope g:", g)

def f(Energy):
    # omnidirectional integral flux above Energy: 4 pi * int_Energy^Emax j dE, in cm-2 s-1
    return 4 * np.pi * j1 * E1 / (g - 1) * ((Energy / E1) ** (1 - g) - (Emax / E1) ** (1 - g))

Energies = [0.13, 1.2, 100.0]  # MeV
# Diff = []
Int = []

# print("Differential")
# for E in Energies:
#     print(f"{E:.4}", f"{-fdiff(E):.4}")
#     print("/gps/hist/point", f"{E:.4}", f"{-fdiff(E):.4}")

print("Integral")
for E in Energies:
    print(f"{E:.4}", f"{f(E):.4}")
    Int.append(f(E))





## Plotting
plt.figure(1, figsize=(6, 6))

# Plot the data points from the function
plt.plot(Energies, Int, 'o-', label="GEO EVT 1-in-150 year peak 6-min", linewidth=2.5, color=Expected)


## Read in Carrington EVT spectrum  -- NO! The new spectrum is in differential units. Cannot be plotted in the integral plot. The function above is the integral spectrum.
# EVT_file = Path + "/Spectra/Carrington/Electron-EVT-2026/Carrington-Electron-Expected-EVT2026-Arb-130keVto1000MeV.mac"
# EVT_Data = readGPSMacro(EVT_file)
# plt.plot(EVT_Data["Energy"], EVT_Data["Flux"], '-', linewidth=2.5, color=Expected)

# Add data point from Meredith 2015 1-in-100 year GOES-W and GOES-E electron fluxes
# GOES-W 1-in-100 year 2 MeV electron flux: 7.68e5 cm−2 s−1 sr−1
# GOES-E 1-in-100 year 2 MeV electron flux: 3.25e5 cm−2 s−1 sr−1
GOES_W_100 = 7.68e5 * 4 * np.pi  # Convert from sr-1 to omnidirectional
GOES_E_100 = 3.25e5 * 4 * np.pi  # Convert from sr-1 to omnidirectional
# plt.plot(2, GOES_W_100, 'x', label="GEO GOES-W 1-in-100 year 1-day avg", color=GOESColor, markersize=8)
# plt.plot(2, GOES_E_100, '+', label="GEO GOES-E 1-in-100 year 1-day avg", color=GOESColor, markersize=8)

# GOES-W 1-in-150 year 2 MeV electron flux: 9.86e5 cm−2 s−1 sr−1
# GOES-E 1-in-150 year 2 MeV electron flux: 4.35e5 cm−2 s−1 sr−1
GOES_W_150 = 9.86e5 * 4 * np.pi  # Convert from sr-1 to omnidirectional
GOES_E_150 = 4.35e5 * 4 * np.pi  # Convert from sr-1 to omnidirectional
plt.plot(2, GOES_W_150, 'x', label="GEO GOES-W 1-in-150 year 1-day avg", color=GOESColor, markersize=8)
plt.plot(2, GOES_E_150, '+', label="GEO GOES-E 1-in-150 year 1-day avg", color=GOESColor, markersize=8)


## FLUMIC spectrum
# FLUMIC Energy in MeV
FLUMIC_Energy = [2.0000E-01, 5.0000E-01, 8.0000E-01, 1.1000E+00, 1.4000E+00, 1.7000E+00, 2.0000E+00, 2.3000E+00, 2.6000E+00, 2.9000E+00, 3.2000E+00, 3.5000E+00, 3.8000E+00, 4.1000E+00, 4.4000E+00, 4.7000E+00, 5.0000E+00, 5.3000E+00, 5.6000E+00, 5.9000E+00]
# Flumic Flux in Flux (m-2 s-1 sr-1)
FLUMIC_Average_Flux = [1.0016E+10, 4.3943E+09, 1.9280E+09, 8.4589E+08, 3.7113E+08, 1.6283E+08, 7.1441E+07, 3.1344E+07, 1.3752E+07, 6.0337E+06, 2.6473E+06, 1.1615E+06, 5.0959E+05, 2.2358E+05, 9.8094E+04, 4.3038E+04, 1.8883E+04, 8.2847E+03, 3.6349E+03, 1.5948E+03]
# Flumic Max Flux in Flux (m-2 s-1 sr-1)
FLUMIC_Max_Flux = [14921000000, 6546400000, 2872200000, 1260200000, 552890000, 242580000, 106430000, 46695000, 20487000, 8988700, 3943800, 1730300, 759160, 333080, 146140, 64116, 28131, 12342, 5415.1, 2375.8]
# Estimated FLUMIC 3 worst case at L = 6.6 (fsc 0.802, foy 0.219), from the FLUMIC 3 equations; illustration only, not a SPENVIS run
FLUMIC_Worst_Flux = [2.8025E+10, 1.5405E+10, 8.4681E+09, 4.6548E+09, 2.5587E+09, 1.4065E+09, 7.7314E+08, 4.2499E+08, 2.3361E+08, 1.2841E+08, 7.0588E+07, 3.8801E+07, 2.1329E+07, 1.1724E+07, 6.4447E+06, 3.5426E+06, 1.9473E+06, 1.0704E+06, 5.8840E+05, 3.2344E+05]

# Convert from m-2 to cm-2 and from sr-1 to omnidirectional
FLUMIC_Average_Flux = np.array(FLUMIC_Average_Flux) * 1e-4 * 4 * np.pi
FLUMIC_Max_Flux = np.array(FLUMIC_Max_Flux) * 1e-4 * 4 * np.pi
FLUMIC_Worst_Flux = np.array(FLUMIC_Worst_Flux) * 1e-4 * 4 * np.pi

# plt.plot(FLUMIC_Energy, FLUMIC_Average_Flux, 'o-', label="GEO FLUMIC 1-day avg", color=FlumicColor, linewidth=1)
plt.plot(FLUMIC_Energy, FLUMIC_Max_Flux, 'x-', label="GEO FLUMIC worst-case 1-day avg", color=FlumicColor, linewidth=1)
# plt.plot(FLUMIC_Energy, FLUMIC_Worst_Flux, 's-', label="GEO FLUMIC worst-case 1-day avg", color=FlumicColor, linewidth=1)






## Read in LEO A9 spectrum
LEO_file = Path + "/Spectra/Carrington/LEO/spenvis_tri.txt"
Protons, Electrons = readSpenvis_tri(LEO_file)
plt.plot(Electrons['Energy'], Electrons['Integral'], '+--', label="LEO AE9 11-year mean", color=LEOColor)

## Read in MEO A9 spectrum
MEO_file = Path + "/Spectra/Carrington/MEO/spenvis_tri.txt"
Protons, Electrons = readSpenvis_tri(MEO_file)
plt.plot(Electrons['Energy'], Electrons['Integral'], '2:', label="MEO AE9 11-year mean", color=MEOColor)

# ## Read in Van-Allen Belt Probes A9 spectrum
# VAP_file = Path + "/Spectra/Carrington/VAP/spenvis_tri.txt"
# Protons, Electrons = readSpenvis_tri(VAP_file)
# plt.plot(Electrons['Energy'], Electrons['Integral'], '2:', label="VAP AE9 11-year mean", color=VAPColor)

## Read in Geostationary A9 spectrum
GEO_file = Path + "/Spectra/Carrington/GEO/spenvis_tri.txt"
Protons, Electrons = readSpenvis_tri(GEO_file)
plt.plot(Electrons['Energy'], Electrons['Integral'], '1-.', label="GEO AE9 11-year mean", color=GEOColor)









""" 
## electron fluxes from paper:
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

plt.plot(EPAMEnergy, IntegralEPAMFlux, label="10/28/03 EPAM electron flux")
plt.plot(PETEnergy, IntegralPETFlux, label="10/28/03 PET electron flux")
 """




## Plot formatting
plt.legend()
plt.yscale("log")
plt.xscale("log")
#plt.title("Electron Flux Spectra")
plt.xlabel("Electron kinetic energy [MeV]")
plt.ylabel("Integral electron flux [cm-2 s-1]")
plt.xlim(0.1, 10)
plt.ylim(1e1, 1e13)
plt.minorticks_on()
plt.grid(axis='x', which='both')
plt.grid(axis='y', which='both')
# New code to create logarithmically spaced Y-axis ticks
y_ticks = np.logspace(1, 12, num=12)  # Generates 12 points between 10^1 and 10^12
plt.yticks(y_ticks)

x_ticks = [0.1, 0.2, 0.5, 1, 2, 5, 10]
plt.xticks(x_ticks, labels=[str(tick) for tick in x_ticks])

plt.savefig(Path + "/Spectra/Carrington/Electron-EVT-2026/Electron-EVT-2026.pdf", format='pdf', bbox_inches="tight")
#plt.show()
