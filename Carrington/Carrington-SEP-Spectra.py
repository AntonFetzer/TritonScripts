import numpy as np
import matplotlib.pyplot as plt
#import sympy as sp
from Read.ReadSpenvis_tri import readSpenvis_tri
from Read.ReadGPSMacro import readGPSMacro
from Read.ReadSpenvis_sef import readSpenvis_sef_protons
from Read.ReadSpenvis_gcf import readSpenvis_gcf
from Read.ReadSpenvis_sefflare import readSpenvis_sefflare

Expected = 'blue' # Blue
PlusColor = 'C1'  # Green
MinusColor = 'C2' # Orange
LEOColor = 'C8'   # yellow
GEOColor = 'C7'   # grey
VAPColor = 'C3'   # Red
CREME96Color = 'C6' # Pink
SAP_PeakColor = 'C9' # Turquoise

## Raw Data from Adnanes E-mail on 13/02/2023 15:53
# Energies = np.array([10, 30, 60, 100, 200])      # MeV

# ExpectedInt = np.array([4.6e4, 8.5e3, 2.3e3, 1.1e3, 2.9e2])     # cm-2 s-1 str-1
# Minus = np.array([2.2e4, 5.6e3, 1.4e3, 4.3e2, 1.1e2])       # cm-2 s-1 str-1
# Plus = np.array([1.3e5, 1.6e4, 5.4e3, 4.4e3, 1.4e3])        # cm-2 s-1 str-1

# ExtrapolatedEnergies = np.array([2, 10])      # MeV
# Slope = (np.log(ExpectedInt[1]) - np.log(ExpectedInt[0])) / (np.log(Energies[1]) - np.log(Energies[0]))
# Intercept = np.log(ExpectedInt[0]) - Slope * np.log(Energies[0])
# ExtrapolatedInt = np.exp(Slope * np.log(ExtrapolatedEnergies) + Intercept)
# print(ExtrapolatedInt)

# # Convert from str-1 to omnidirectional
# ExpectedInt *= 4 * np.pi
# Minus *= 4 * np.pi
# Plus *= 4 * np.pi

### Carrington Integral GPS Macros ###
# This was just to verify that the GPS macros agree with the raw data from Adnanes E-mail
# They agree, therefore do not need to be plotted again

## Read in GPS Macro
ExpectedInt = "/l/triton_work/Spectra/Carrington/SEP-Final/Carrington-SEP-Expected-Int-With0.mac"
ExpectedIntWith0 = "/l/triton_work/Spectra/Carrington/SEP-Final/Carrington-SEP-Expected-Int-With0.mac"
Minus2SigmaIntWith0 = "/l/triton_work/Spectra/Carrington/SEP-Final/Carrington-SEP-Minus2Sigma-Int-With0.mac"
Plus2SigmaIntWith0 = "/l/triton_work/Spectra/Carrington/SEP-Final/Carrington-SEP-Plus2Sigma-Int-With0.mac"

ExpectedInt = readGPSMacro(ExpectedInt)
ExpectedIntWith0 = readGPSMacro(ExpectedIntWith0)
Minus2SigmaIntWith0 = readGPSMacro(Minus2SigmaIntWith0)
Plus2SigmaIntWith0 = readGPSMacro(Plus2SigmaIntWith0)

Energies = ExpectedInt['Energy']
ExpectedInt = ExpectedInt['Flux']
Minus = Minus2SigmaIntWith0['Flux']
Plus = Plus2SigmaIntWith0['Flux']


plt.figure(1, figsize=(8, 8))

plt.fill_between(Energies, Plus, ExpectedInt, color=PlusColor, alpha=0.5, label="Carrington SEP EVT +2 Sigma Peak Flux")
#plt.plot(Energies, Plus, '.-', label="Carrington SEP EVT +2 Sigma", color=PlusColor)
#plt.plot(Energies, ExpectedInt, '-', color=Expected, linewidth=2, markersize=10)
# Plot the expected integral flux as error bars with the plus and minus 2 sigma.
plt.errorbar(Energies, ExpectedInt, yerr=[ExpectedInt - Minus, Plus - ExpectedInt], fmt=' ', label="Carrington SEP EVT Expected Peak Flux", color=Expected, linewidth=2, markersize=10, capsize=5, capthick=2, zorder=2)

#plt.plot(Energies, Minus, '.-', label="Carrington SEP EVT -2 Sigma", color=MinusColor)
plt.fill_between(Energies, ExpectedInt, Minus, color=MinusColor, alpha=0.5, label="Carrington SEP EVT -2 Sigma Peak Flux")

# Plot the extrapolated values
#plt.plot(ExtrapolatedEnergies, ExtrapolatedInt, '-o', label="Extrapolated Values", color='C4', markersize=10)

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

## Read in LEO A9 spectrum
LEO_AE9 = "/l/triton_work/Spectra/Carrington/LEO/spenvis_tri.txt"
LEO_AE9_Protons, LEO_AE9_Electrons = readSpenvis_tri(LEO_AE9)
plt.plot(LEO_AE9_Protons['Energy'], LEO_AE9_Protons['Integral'], '.-', label="AP9 mean LEO Trapped Proton Flux", color=LEOColor)

## Read in GEO A9 spectrum
# GEO_AE9 = "/l/triton_work/Spectra/Carrington/GEO/spenvis_tri.txt"
# GEO_AE9_Protons, GEO_AE9_Electrons = readSpenvis_tri(GEO_AE9)
# plt.plot(GEO_AE9_Protons['Energy'], GEO_AE9_Protons['Integral'], '.-', label="AP9 mean GEO Trapped Proton Flux", color=GEOColor)
# Flux at low energy is high, but then drops to zero at 8 MeV, therefore not plotted

## Read in VAP A9 spectrum
VAP_AE9 = "/l/triton_work/Spectra/Carrington/VAP/spenvis_tri.txt"
VAP_AE9_Protons, VAP_AE9_Electrons = readSpenvis_tri(VAP_AE9)
plt.plot(VAP_AE9_Protons['Energy'], VAP_AE9_Protons['Integral'], '.-', label="AP9 mean VAP Trapped Proton Flux", color=VAPColor)


### SAPPHIRE Solar Proton Spectra ###

## Read in LEO SEP spectrum
# LEO_SEP = "/l/triton_work/Spectra/Carrington/LEO/spenvis_sef.txt"
# LEO_SEP_Data = readSpenvis_sef_protons(LEO_SEP)
# LEO_SEP_Fluence = LEO_SEP_Data['IFluence']
# Convert from Fluence to Flux. Mission duration is 4015 days
# LEO_SEP_Flux = LEO_SEP_Fluence / (4015 * 24 * 3600)

# plt.plot(LEO_SEP_Data['Energy'], LEO_SEP_Flux, '+--', label="SAPPHIRE mean LEO Solar Proton Flux", color=LEOColor)
# Earth mangetosphere shields LEO from SEP, therefore flux is below 1e-3 cm-2 s-1, therefore not plotted

# ## Read in Geostationary SEP spectrum
GEO_SEP = "/l/triton_work/Spectra/Carrington/GEO/spenvis_sef.txt"
GEO_SEP_Data = readSpenvis_sef_protons(GEO_SEP)
GEO_SEP_Fluence = GEO_SEP_Data['IFluence']
# Convert from Fluence to Flux. Mission duration is 4015 days
GEO_SEP_Flux = GEO_SEP_Fluence / (4015 * 24 * 3600)

plt.plot(GEO_SEP_Data['Energy'], GEO_SEP_Flux, '+--', label="SAPPHIRE mean GEO Solar Proton Flux", color=GEOColor)

# ## Read in VAP SEP spectrum
# VAP_SEP = "/l/triton_work/Spectra/Carrington/VAP/spenvis_sef.txt"
# VAP_SEP_Data = readSpenvis_sef_protons(VAP_SEP)
# VAP_SEP_Fluence = VAP_SEP_Data['IFluence']
# Convert from Fluence to Flux. Mission duration is 4015 days
# VAP_SEP_Flux = VAP_SEP_Fluence / (4015 * 24 * 3600)
# plt.plot(VAP_SEP_Data['Energy'], VAP_SEP_Flux, '+--', label="SAPPHIRE mean VAP Solar Proton Flux", color=VAPColor)
# Very similar to the GEO SEP, therefore not plotted


### CREME96 Peak 5min Solar Proton Flux ###
# Read in GEO CREME96 spectrum
GEO_Flare = "/l/triton_work/Spectra/Carrington/GEO-Extreme/CREME96/spenvis_sefflare.txt"
GEO_Flare_Data = readSpenvis_sefflare(GEO_Flare)
plt.plot(GEO_Flare_Data['Energy'], GEO_Flare_Data['IFlux'], 'o', label="CREME96 Peak 5 min Flux", color=CREME96Color, linewidth=5, zorder=3)

## SAPPHIRE 1 in n year Solar Proton Event Peak Fluxes
# Read in GEO 1 in 100 year SEP Peak Flux
GEO_SEP_100yr = "/l/triton_work/Spectra/Carrington/GEO-Extreme/SAPPHIRE100yearPeakFlux/spenvis_sefflare.txt"
GEO_SEP_100yr_Data = readSpenvis_sefflare(GEO_SEP_100yr)
plt.plot(GEO_SEP_100yr_Data['Energy'], GEO_SEP_100yr_Data['IFlux'], 'P', label="SAPPHIRE 1-in-100 year Proton Peak Flux", color=SAP_PeakColor, zorder=4)

# Read in GEO 1 in 300 year SEP Peak Flux
GEO_SEP_300yr = "/l/triton_work/Spectra/Carrington/GEO-Extreme/SAPPHIRE300yearPeakFlux/spenvis_sefflare.txt"
GEO_SEP_300yr_Data = readSpenvis_sefflare(GEO_SEP_300yr)
plt.plot(GEO_SEP_300yr_Data['Energy'], GEO_SEP_300yr_Data['IFlux'], 'X', label="SAPPHIRE 1-in-300 year Proton Peak Flux", color=SAP_PeakColor, zorder=4)


### ISO Cosmic Spectra ###

## Read in LEO Cosmic spectrum
# LEO_Cosmic = "/l/triton_work/Spectra/Carrington/LEO/spenvis_gcf.txt"
# LEO_Cosmic_Data = readSpenvis_gcf(LEO_Cosmic)
# plt.plot(LEO_Cosmic_Data[0, :, 0], LEO_Cosmic_Data[0, :, 1], '*:', label="ISO mean LEO Cosmic Proton Flux", color=LEOColor)
# Cosmics are below 1 cm-2 s-1, therefore not plotted

# ## Read in Geostationary Cosmic spectrum
GEO_Cosmic = "/l/triton_work/Spectra/Carrington/GEO/spenvis_gcf.txt"
GEO_Cosmic_Data = readSpenvis_gcf(GEO_Cosmic)
plt.plot(GEO_Cosmic_Data[0, :, 0], GEO_Cosmic_Data[0, :, 1], '*:', label="ISO mean GEO Cosmic Proton Flux", color=GEOColor)

# ## Read in Van-Allen Belt Probes Cosmic spectrum
# VAP_Cosmic = "/l/triton_work/Spectra/Carrington/VAP/spenvis_gcf.txt"
# VAP_Cosmic_Data = readSpenvis_gcf(VAP_Cosmic)
# plt.plot(VAP_Cosmic_Data[0, :, 0], VAP_Cosmic_Data[0, :, 1], '*:', label="ISO mean VAP Cosmic Proton Flux", color=VAPColor)
# Very similar to the GEO Cosmic, therefore not plotted

plt.xlim(8, 250)
plt.ylim(1, 2e6)
plt.yscale("log")
plt.xscale("log")
plt.title("Carrington Solar Proton Flux Comparison")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Integral Flux [cm-2 s-1]")
plt.grid(which='both')

# Adjust legend to ensure proper order
handles, labels = plt.gca().get_legend_handles_labels()
order = [0, 9, 1, 5, 6, 7, 2, 3, 4, 8]  # Adjust this list to reorder as needed
plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='lower right')

plt.savefig("/l/triton_work/Spectra/Carrington/SEP-Final/ProtonSpectrumComparison.pdf", format='pdf', bbox_inches="tight")
# plt.show()