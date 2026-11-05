import numpy as np
import matplotlib.pyplot as plt
#import sympy as sp

import sys, os
# Ensure the package root (the "Python" folder) is on sys.path so "from Read.* import ..." works
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from Read.ReadSpenvis_tri import readSpenvis_tri
from Read.ReadGPSMacro import readGPSMacro
from Read.ReadSpenvis_sef import readSpenvis_sef_protons
from Read.ReadSpenvis_gcf import readSpenvis_gcf
from Read.ReadSpenvis_sefflare import readSpenvis_sefflare
from Dependencies.CarringtonColors import Expected, PlusColor, MinusColor, LEOColor, GEOColor, MEOColor, CREME96Color, SAP_PeakColor

Path = "/scratch/work/fetzera1/Spectra/Carrington"

## Raw Data from Adnanes E-mail on 13/02/2023 15:53
Energies = np.array([10, 30, 60, 100, 200])      # MeV

ExpectedInt = np.array([4.6e4, 8.5e3, 2.3e3, 1.1e3, 2.9e2])     # cm-2 s-1 str-1
Minus = np.array([2.2e4, 5.6e3, 1.4e3, 4.3e2, 1.1e2])       # cm-2 s-1 str-1
Plus = np.array([1.3e5, 1.6e4, 5.4e3, 4.4e3, 1.4e3])        # cm-2 s-1 str-1

# Convert from str-1 to omnidirectional
ExpectedInt *= 4 * np.pi
Minus *= 4 * np.pi
Plus *= 4 * np.pi



plt.figure(1, figsize=(8, 7))

plt.plot(Energies, ExpectedInt, '-', color=Expected, linewidth=2)
# Plot the expected integral flux as error bars with the plus and minus 2 sigma.
plt.errorbar(Energies, ExpectedInt, yerr=[ExpectedInt - Minus, Plus - ExpectedInt], fmt=' ', label="GEO EVT 1-in-150 year peak SEP", color=Expected, linewidth=2, markersize=10, capsize=5, capthick=2, zorder=2)

#plt.plot(Energies, Minus, '.-', label="Carrington SEP EVT -2 Sigma", color=MinusColor)
plt.fill_between(Energies, ExpectedInt, Minus, color=MinusColor, alpha=0.5, label="GEO EVT 1-in-150 year peak SEP -2\u03C3")
#plt.plot(Energies, Plus, '.-', label="Carrington SEP EVT +2 Sigma", color=PlusColor)
plt.fill_between(Energies, Plus, ExpectedInt, color=PlusColor, alpha=0.5, label="GEO EVT 1-in-150 year peak SEP +2\u03C3")



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
LEO_AE9 = Path + "/LEO/spenvis_tri.txt"
LEO_AE9_Protons, LEO_AE9_Electrons = readSpenvis_tri(LEO_AE9)
plt.plot(LEO_AE9_Protons['Energy'], LEO_AE9_Protons['Integral'], '.-', 
         label="LEO AP9 11-year mean trapped proton", color=LEOColor)

## Read in MEO A9 spectrum
# MEO_AE9 = Path + "/MEO/spenvis_tri.txt"
# MEO_AE9_Protons, MEO_AE9_Electrons = readSpenvis_tri(MEO_AE9)
# plt.plot(MEO_AE9_Protons['Energy'], MEO_AE9_Protons['Integral'], '.-', 
#          label="MEO AP9 11-year mean", color=MEOColor)

## Read in GEO A9 spectrum
# GEO_AE9 = Path + "/GEO/spenvis_tri.txt"
# GEO_AE9_Protons, GEO_AE9_Electrons = readSpenvis_tri(GEO_AE9)
# plt.plot(GEO_AE9_Protons['Energy'], GEO_AE9_Protons['Integral'], '.-', 
#          label="GEO AP9 11-year mean", color=GEOColor)
# Flux at low energy is high, but then drops to zero at 8 MeV, therefore not plotted

## Read in VAP A9 spectrum
# VAP_AE9 = Path + "/VAP/spenvis_tri.txt"
# VAP_AE9_Protons, VAP_AE9_Electrons = readSpenvis_tri(VAP_AE9)
# plt.plot(VAP_AE9_Protons['Energy'], VAP_AE9_Protons['Integral'], '.-', 
# label="VAP AP9 11-year mean", color=VAPColor)




### SAPPHIRE Solar Proton Spectra ###

## Read in LEO SEP spectrum
# LEO_SEP = Path + "/LEO/spenvis_sef.txt"
# LEO_SEP_Data = readSpenvis_sef_protons(LEO_SEP)
# plt.plot(LEO_SEP_Data['Energy'], LEO_SEP_Data['IFlux'], '+--',
#           label="LEO SAPPHIRE 11-year mean SEP", color=LEOColor)
# Earth mangetosphere shields LEO from SEP, therefore flux is below 1e-3 cm-2 s-1, therefore not plotted

# ## Read in MEO SEP spectrum
MEO_SEP = Path + "/MEO/spenvis_sef.txt"
MEO_SEP_Data = readSpenvis_sef_protons(MEO_SEP)
plt.plot(MEO_SEP_Data['Energy'], MEO_SEP_Data['IFlux'], '+--',
          label="MEO SAPPHIRE 11-year mean SEP", color=MEOColor)

# ## Read in Geostationary SEP spectrum
GEO_SEP = Path + "/GEO/spenvis_sef.txt"
GEO_SEP_Data = readSpenvis_sef_protons(GEO_SEP)
plt.plot(GEO_SEP_Data['Energy'], GEO_SEP_Data['IFlux'], '+--',
          label="GEO SAPPHIRE 11-year mean SEP", color=GEOColor)

# ## Read in VAP SEP spectrum
# VAP_SEP = Path + "/VAP/spenvis_sef.txt"
# VAP_SEP_Data = readSpenvis_sef_protons(VAP_SEP)
# plt.plot(VAP_SEP_Data['Energy'], VAP_SEP_Data['IFlux'], '+--',
#  label="VAP SAPPHIRE 11-year mean SEP", color=VAPColor)
# Very similar to the GEO SEP, therefore not plotted




### CREME96 Peak 5min Solar Proton Flux ###
# Read in GEO CREME96 spectrum
GEO_Flare = Path + "/GEO-Extreme/CREME96/spenvis_sefflare.txt"
GEO_Flare_Data = readSpenvis_sefflare(GEO_Flare)
plt.plot(GEO_Flare_Data['Energy'], GEO_Flare_Data['IFlux'], 'o', label="GEO CREME96 peak 5 min SEP", color=CREME96Color, linewidth=5, zorder=3)

## SAPPHIRE 1 in n year Solar Proton Event Peak Fluxes
# Read in GEO 1 in 100 year SEP Peak Flux
GEO_SEP_100yr = Path + "/GEO-Extreme/SAPPHIRE100yearPeakFlux/spenvis_sefflare.txt"
GEO_SEP_100yr_Data = readSpenvis_sefflare(GEO_SEP_100yr)
plt.plot(GEO_SEP_100yr_Data['Energy'], GEO_SEP_100yr_Data['IFlux'], 'P', label="GEO SAPPHIRE 1-in-100 year peak SEP", color=SAP_PeakColor, zorder=4)

# Read in GEO 1 in 300 year SEP Peak Flux
GEO_SEP_300yr = Path + "/GEO-Extreme/SAPPHIRE300yearPeakFlux/spenvis_sefflare.txt"
GEO_SEP_300yr_Data = readSpenvis_sefflare(GEO_SEP_300yr)
plt.plot(GEO_SEP_300yr_Data['Energy'], GEO_SEP_300yr_Data['IFlux'], 'X', label="GEO SAPPHIRE 1-in-300 year peak SEP", color=SAP_PeakColor, zorder=4)


### ISO Cosmic Spectra ###

## Read in LEO Cosmic spectrum
# LEO_Cosmic = Path + "/LEO/spenvis_gcf.txt"
# LEO_Cosmic_Data = readSpenvis_gcf(LEO_Cosmic)
# plt.plot(LEO_Cosmic_Data[0, :, 0], LEO_Cosmic_Data[0, :, 1], '*:', label="LEO ISO 11-year mean cosmic proton", color=LEOColor)
# Cosmics are below 1 cm-2 s-1, therefore not plotted

# ## Read in MEO Cosmic spectrum
# MEO_Cosmic = Path + "/MEO/spenvis_gcf.txt"
# MEO_Cosmic_Data = readSpenvis_gcf(MEO_Cosmic)
# plt.plot(MEO_Cosmic_Data[0, :, 0], MEO_Cosmic_Data[0, :, 1], '*:', label="MEO ISO 11-year mean cosmic proton", color=MEOColor)

# ## Read in Geostationary Cosmic spectrum
GEO_Cosmic = Path + "/GEO/spenvis_gcf.txt"
GEO_Cosmic_Data = readSpenvis_gcf(GEO_Cosmic)
plt.plot(GEO_Cosmic_Data[0, :, 0], GEO_Cosmic_Data[0, :, 1], '*:', label="GEO ISO 11-year mean cosmic proton", color=GEOColor)

# ## Read in Van-Allen Belt Probes Cosmic spectrum
# VAP_Cosmic = Path + "/VAP/spenvis_gcf.txt"
# VAP_Cosmic_Data = readSpenvis_gcf(VAP_Cosmic)
# plt.plot(VAP_Cosmic_Data[0, :, 0], VAP_Cosmic_Data[0, :, 1], '*:', label="VAP ISO 11-year mean cosmic proton", color=VAPColor)
# Very similar to the GEO Cosmic, therefore not plotted



plt.xlim(8, 250)
plt.ylim(1, 2e6)
plt.yscale("log")
plt.xscale("log")
plt.title("Proton Fluxes")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Integral Flux [cm-2 s-1]")
plt.grid(which='both')

# # Adjust legend to ensure proper order
handles, labels = plt.gca().get_legend_handles_labels()
order = [0, 9, 1, 5, 6, 7, 4, 3, 8, 2]  # Adjust this list to reorder as needed
plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='lower right')
# plt.legend()


plt.savefig(Path + "/SEP-EVT/ProtonSpectra.pdf", format='pdf', bbox_inches="tight")
# plt.show()