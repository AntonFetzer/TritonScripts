import numpy as np
import matplotlib.pyplot as plt
#import sympy as sp
from GRAS.Read.ReadSpenvis_tri import readSpenvis_tri
from GRAS.Read.ReadGPSMacro import readGPSMacro
from GRAS.Read.ReadSpenvis_sef import readSpenvis_sef
from GRAS.Read.ReadSpenvis_gcf import readSpenvis_gcf
from GRAS.Read.ReadSpenvis_sefflare import readSpenvis_sefflare

ExpectedIntColor = 'blue' # Blue
PlusColor = 'C1'        # Green
MinusColor = 'C2'       # Orange
ISSColor = 'C8'         # yellow
GEOColor = 'C7'         # grey
VABColor = 'C3'         # Red

## Raw Data from Adnanes E-mail on 13/02/2023 15:53
Energies = np.array([10, 30, 60, 100, 200])      # MeV

ExpectedInt = np.array([4.6e4, 8.5e3, 2.3e3, 1.1e3, 2.9e2])     # cm-2 s-1 str-1
Minus = np.array([2.2e4, 5.6e3, 1.4e3, 4.3e2, 1.1e2])       # cm-2 s-1 str-1
Plus = np.array([1.3e5, 1.6e4, 5.4e3, 4.4e3, 1.4e3])        # cm-2 s-1 str-1

ExtrapolatedEnergies = np.array([2, 10])      # MeV
Slope = (np.log(ExpectedInt[1]) - np.log(ExpectedInt[0])) / (np.log(Energies[1]) - np.log(Energies[0]))
Intercept = np.log(ExpectedInt[0]) - Slope * np.log(Energies[0])
ExtrapolatedInt = np.exp(Slope * np.log(ExtrapolatedEnergies) + Intercept)
print(ExtrapolatedInt)

# Convert from str-1 to omnidirectional
ExpectedInt *= 4 * np.pi
Minus *= 4 * np.pi
Plus *= 4 * np.pi



plt.figure(1, figsize=(8, 8))

plt.fill_between(Energies, Plus, ExpectedInt, color=PlusColor, alpha=0.5, label="Carrington SEP EVT +2 Sigma")
#plt.plot(Energies, Plus, '.-', label="Carrington SEP EVT +2 Sigma", color=PlusColor)
plt.plot(Energies, ExpectedInt, '-', color=ExpectedIntColor, linewidth=2, markersize=10)
# Plot the expected integral flux as error bars with the plus and minus 2 sigma.
plt.errorbar(Energies, ExpectedInt, yerr=[ExpectedInt - Minus, Plus - ExpectedInt], fmt=' ', label="Carrington SEP EVT Expected", color=ExpectedIntColor, linewidth=2, markersize=10, capsize=5, capthick=2, zorder=2)

#plt.plot(Energies, Minus, '.-', label="Carrington SEP EVT -2 Sigma", color=MinusColor)
plt.fill_between(Energies, ExpectedInt, Minus, color=MinusColor, alpha=0.5, label="Carrington SEP EVT -2 Sigma")

# Plot the extrapolated values
plt.plot(ExtrapolatedEnergies, ExtrapolatedInt, '-o', label="Extrapolated Values", color='C4', markersize=10)

""" 
### 2003 SPE Data ###

# Numbers from Mewaldt et. al.
# "Proton, helium, and electron spectra during the large solar particle events of October–November 2003
# doi:10.1029/2005JA011038
SPE2003_Mewaldt_Fluence = np.array([11424063484.6895, 2705029335.12195, 667289524.189153, 184869107.24167, 10435956.1620884])    # cm-2 Event Fluence 
# Convert to flux estimate by dividing by the duration of the event, which is roughly 1 day
SPE2003_Mewaldt_Flux = SPE2003_Mewaldt_Fluence / (24 * 3600)    # cm-2 s-1
# Plot the flux
plt.plot(Energies, SPE2003_Mewaldt_Flux, '--', label="SPE 2003 Mewaldt et. al.", color='C8', linewidth=4)


## Numbers from Jiggens et. al. 
# "The magnitude and effects of extreme solar particle events"
# DOI: 10.1051/swsc/2014017 Figue 9
SPE2003_Jiggens_Fluence = np.array([1.5678e+10, 3.4589e+9, 7.0664e+8, 1.5325e+8, 1.9409e+7])    # cm-2 Event Fluence
# Convert to flux estimate by dividing by the duration of the event, which is roughly 1 day
SPE2003_Jiggens_Flux = SPE2003_Jiggens_Fluence / (24 * 3600)    # cm-2 s-1
# Plot the flux
plt.plot(Energies, SPE2003_Jiggens_Flux, '.-', label="SPE 2003 Jiggens et. al.", color='C3', linewidth=2)


## Numbers from Townsend et. al. 
#"Carrington Flare of 1859 as a Prototypical Worst-Case Solar Energetic Particle Event"
# DOI: 10.1109/TNS.2003.821602
Carrington_Townsend_Mar91 = np.array([1.0504e+11, 1.8805e+10, 4.0681e+9, 8.2153e+8, 6.6666e+7])    # cm-2 Event Fluence
# Convert to flux estimate by dividing by the duration of the event, which is roughly 1 day
Carrington_Townsend_Flux = Carrington_Townsend_Mar91 / (24 * 3600)    # cm-2 s-1
# Plot the flux
plt.plot(Energies, Carrington_Townsend_Flux, ':', label="Carrington Estimate Townsend et. al.", color='C9', linewidth=4)

 """
### A9 Spectra ###

## Read in ISS A9 spectrum
ISS_AE9 = "/l/triton_work/Spectra/ISS/spenvis_tri.txt"
ISS_AE9_Protons, ISS_AE9_Electrons = readSpenvis_tri(ISS_AE9)
plt.plot(ISS_AE9_Protons['Energy'], ISS_AE9_Protons['Integral'], '.-', label="AP9 LEO Trapped Proton Flux", color=ISSColor)

## Read in Geostationary A9 spectrum
# GEO_AE9 = "/l/triton_work/Spectra/GEO/spenvis_tri.txt"
# GEO_AE9_Protons, GEO_AE9_Electrons = readSpenvis_tri(GEO_AE9)
# plt.plot(GEO_AE9_Protons['Energy'], GEO_AE9_Protons['Integral'], '.-', label="AP9 GEO Trapped Proton Flux", color=GEOColor)
# Flux at low energy is high, but then drops to zero at 8 MeV, therefore not plotted

## Read in Van-Allen Belt Probes A9 spectrum
VAB_AE9 = "/l/triton_work/Spectra/VAB/spenvis_tri.txt"
VAB_AE9_Protons, VAB_AE9_Electrons = readSpenvis_tri(VAB_AE9)
plt.plot(VAB_AE9_Protons['Energy'], VAB_AE9_Protons['Integral'], '.-', label="AP9 Van-Allen-Belt Trapped Proton Flux", color=VABColor)


### SAPPHIRE Solar Proton Spectra ###

## Read in ISS SEP spectrum
# ISS_SEP = "/l/triton_work/Spectra/ISS/spenvis_sef.txt"
# ISS_SEP_Data = readSpenvis_sef(ISS_SEP)
# plt.plot(ISS_SEP_Data[0, :, 0], ISS_SEP_Data[0, :, 1], '+--', label="SAPPHIRE LEO Solar Proton Flux", color=ISSColor)
# Earth mangetosphere shields LEO from SEP, therefore flux is below 1e-3 cm-2 s-1, therefore not plotted

# ## Read in Geostationary SEP spectrum
GEO_SEP = "/l/triton_work/Spectra/GEO/spenvis_sef.txt"
GEO_SEP_Data = readSpenvis_sef(GEO_SEP)
plt.plot(GEO_SEP_Data[0, :, 0], GEO_SEP_Data[0, :, 1], '+--', label="SAPPHIRE GEO Solar Proton Flux", color=GEOColor)

# ## Read in Van-Allen Belt Probes SEP spectrum
# VAB_SEP = "/l/triton_work/Spectra/Van-Allen-Probes/spenvis_sef.txt"
# VAB_SEP_Data = readSpenvis_sef(VAB_SEP)
# plt.plot(VAB_SEP_Data[0, :, 0], VAB_SEP_Data[0, :, 1], '+--', label="SAPPHIRE Van-Allen-Belt Solar Proton Flux", color=VABColor)
# Very similar to the GEO SEP, therefore not plotted


### ISO Cosmic Spectra ###

## Read in ISS Cosmic spectrum
# ISS_Cosmic = "/l/triton_work/Spectra/ISS/spenvis_gcf.txt"
# ISS_Cosmic_Data = readSpenvis_gcf(ISS_Cosmic)
# plt.plot(ISS_Cosmic_Data[0, :, 0], ISS_Cosmic_Data[0, :, 1], '*:', label="ISO ISS Cosmic Proton Flux", color=ISSColor)
# Cosmics are below 1 cm-2 s-1, therefore not plotted

# ## Read in Geostationary Cosmic spectrum
GEO_Cosmic = "/l/triton_work/Spectra/GEO/spenvis_gcf.txt"
GEO_Cosmic_Data = readSpenvis_gcf(GEO_Cosmic)
plt.plot(GEO_Cosmic_Data[0, :, 0], GEO_Cosmic_Data[0, :, 1], '*:', label="ISO GEO Cosmic Proton Flux", color=GEOColor)

# ## Read in Van-Allen Belt Probes Cosmic spectrum
# VAB_Cosmic = "/l/triton_work/Spectra/Van-Allen-Probes/spenvis_gcf.txt"
# VAB_Cosmic_Data = readSpenvis_gcf(VAB_Cosmic)
# plt.plot(VAB_Cosmic_Data[0, :, 0], VAB_Cosmic_Data[0, :, 1], '*:', label="ISO Van-Allen-Belt Cosmic Proton Flux", color=VABColor)
# Very similar to the GEO Cosmic, therefore not plotted


### CREME96 Peak 5min Solar Proton Flux ###
# ## Read in GEO CREME96 spectrum
GEO_Flare = "/l/triton_work/Spectra/GEO/spenvis_sefflare.txt"
GEO_Flare_Data = readSpenvis_sefflare(GEO_Flare)
plt.plot(GEO_Flare_Data['Energy'], GEO_Flare_Data['IFlux'], 'o', label="CREME96 GEO Peak 5min Flux", color=GEOColor, linewidth=5, zorder=3)

for energy, flux in zip(GEO_Flare_Data['Energy'], GEO_Flare_Data['IFlux']):
    print(f"{energy} {flux}")


### Carrington Integral GPS Macros ###
# This was just to verify that the GPS macros agree with the raw data from Adnanes E-mail
# They agree, therefore do not need to be plotted again

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


plt.xlim(1, 250)
plt.ylim(1, 1e7)
plt.yscale("log")
plt.xscale("log")
plt.title("Carrington Solar Energetic Proton Comparison")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Integral Flux [cm-2 s-1]")


# Adjust legend to ensure proper order
handles, labels = plt.gca().get_legend_handles_labels()
order = [0, 7, 1, 2, 3, 4, 5, 6]  # Adjust this list to reorder as needed
plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='lower right')

#plt.legend(loc='lower right')



plt.grid(which='both')

plt.savefig("/l/triton_work/Spectra/Carrington/SEP-Final/ProtonSpectrumComparison.pdf", format='pdf', bbox_inches="tight")
# plt.show()