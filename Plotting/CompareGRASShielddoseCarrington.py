import numpy as np
import matplotlib.pyplot as plt
from Read.ReadSD2Q import readSDQ2
from Dependencies.TotalDose import totalDose
from Dependencies.MergeTotalDose import mergeTotalDose


plt.figure(0, [4, 6])
# Aluminium Shielding Thickness [mm]
# 101 Tiles from 0 to 10 mm
x = np.linspace(0, 10, num=101, dtype=float, endpoint=True)
# print(x)
####### Plot 1kRad line #########
CriticalDose = [1 for i in x]
plt.plot(x, CriticalDose, color='k', linewidth=2, label='1 krad')
CriticalDose = [10 for i in x]
plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='10 krad')
# CriticalDose = [100 for i in x]
# plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='100 krad')


##################################### Carrington Spectra ########################################

Path = "/l/triton_work/Shielding_Curves/Carrington/"
# The carrington fluxes are per second, but we are interested in the event fluence which corresponds to 1.23 d
PerStoEvent = 1.23 * 24 * 60 * 60

#### Carrington SEP fluxes ######
Ex = totalDose(Path + "Carrington-SEP-Expected-Int/Res")
Minus = totalDose(Path + "Carrington-SEP-Minus2Sigma-Int/Res")
Plus = totalDose(Path + "Carrington-SEP-Plus2Sigma-Int/Res")

# Plot Carrington SEP fluxes
plt.errorbar(x, Ex['dose'] * PerStoEvent, Ex['error'] * PerStoEvent, fmt='', label='Carrington SEP Event Dose', 
             color=Expected, linestyle='', capsize=1.5, capthick=2)
plt.fill_between(x, Plus['dose'] * PerStoEvent, Ex['dose'] * PerStoEvent, color='C1', alpha=0.5, 
                 label='Carrington SEP + 2\u03C3')
plt.fill_between(x, Ex['dose'] * PerStoEvent, Minus['dose'] * PerStoEvent, color='C2', alpha=0.5, 
                 label='Carrington SEP − 2\u03C3')

# Alternative SEP spectra
# Carrington_SEP_Expected_Diff = totalDose(Path + "Carrington-SEP-Expected-Diff/Res")
# plt.errorbar(x, Carrington_SEP_Expected_Diff['dose'] * PerStoEvent, Carrington_SEP_Expected_Diff['error'] * PerStoEvent,
#               markersize=2, capsize=2, label='Carrington SEP Expected Diff')

# Carrington_SEP_Expected_Int_With0 = totalDose(Path + "Carrington-SEP-Expected-Int-With0/Res")
# plt.errorbar(x, Carrington_SEP_Expected_Int_With0['dose'] * PerStoEvent, 
#              Carrington_SEP_Expected_Int_With0['error'] * PerStoEvent, markersize=2, capsize=2, 
#              label='Carrington SEP Expected Int With 0')

# MinusWith0 = totalDose(Path + "Carrington-SEP-Minus2Sigma-Int-With0/Res")
# plt.errorbar(x, MinusWith0['dose'] * PerStoEvent, MinusWith0['error'] * PerStoEvent, fmt='', markersize=2, capsize=2, 
#              label='Carrington SEP − 2\u03C3 With 0', color='orange')

# PlusWith0 = totalDose(Path + "Carrington-SEP-Plus2Sigma-Int-With0/Res")
# plt.errorbar(x, PlusWith0['dose'] * PerStoEvent, PlusWith0['error'] * PerStoEvent, fmt='', markersize=2, capsize=2, 
#              label='Carrington SEP + 2\u03C3 With 0', color='green')



#### Carrington Electron fluxes ######
# CarringtonElectronINTEGRALPowTabelated = totalDose(Path + "CarringtonElectronINTEGRALPowTabelated/Res")
# Alternative Electron spectra
# CarringtonElectronDiffPowTabelated = totalDose(Path + "CarringtonElectronDiffPowTabelated/Res")
# CarringtonElectronDiffPow = totalDose(Path + "CarringtonElectronDiffPow/Res")

# Plot Carrington Electron fluxes
# plt.errorbar(x, CarringtonElectronINTEGRALPowTabelated['dose'] * PerStoEvent, 
            #  yerr=CarringtonElectronINTEGRALPowTabelated['error'] * PerStoEvent, label='Carrington Electron\nEvent Dose', 
            #  color=Expected, linestyle='', capsize=1.5, capthick=2)



################################# Other extreme Spectra ########################################

## CREME96 Peak Solar Proton Flux ##
# solar_proton_5minPeakFlux = totalDose(Path + "Solar-Proton-5minPeakFlux/Res")
# solar_proton_5minPeakFlux_New = totalDose(Path + "Solar-Proton-5minPeakFlux_New/Res")

# Plot Solar Proton 5min Peak Flux
# The flux is in 1/s cm2, therefore we need to multiply by to get the event dose in krad
# plt.errorbar(x, solar_proton_5minPeakFlux['dose'] * 60 * 60 * 24 * 1.23, 
#              yerr=solar_proton_5minPeakFlux['error'] * 60 * 60 * 24 * 1.23, capsize=3, linestyle='', color='C4', 
#              label="1 month peak flux on GEO")
# plt.errorbar(x, solar_proton_5minPeakFlux_New['dose'] * 60 * 60 * 24 * 1.23, 
#              yerr=solar_proton_5minPeakFlux_New['error'] * 60 * 60 * 24 * 1.23, capsize=3, linestyle='', color='C9', 
#              label="1 month peak flux on GEO")


################################## Nominal Mission Spectra ########################################

#### LEO ######
LEO_electron = totalDose(Path + "LEO-electron/Res")
LEO_trapped_proton = totalDose(Path + "LEO-trapped-proton/Res")
LEO_solar_proton = totalDose(Path + "LEO-solar-proton/Res") # Insigngificant on LEO
LEO_cosmic_proton = totalDose(Path + "LEO-cosmic-proton/Res") # Negiligible
LEO_cosmic_iron = totalDose(Path + "LEO-cosmic-iron/Res") # Negiligible

# plt.errorbar(x, LEO_electron['dose'], yerr=LEO_electron['error'], label="LEO Electrons", 
#              linestyle='', capsize=1, color='C0')
# plt.errorbar(x, LEO_trapped_proton['dose'], yerr=LEO_trapped_proton['error'], label="LEO Trapped Protons", 
#              linestyle='', capsize=1, color='C1')
# plt.errorbar(x, LEO_solar_proton['dose'], yerr=LEO_solar_proton['error'], label="LEO Solar Protons", 
#              linestyle='', capsize=1, color='C2')
# plt.errorbar(x, LEO_cosmic_proton['dose'], yerr=LEO_cosmic_proton['error'], label="LEO Cosmic Protons", 
#              linestyle='', capsize=1, color='C3')
# plt.errorbar(x, LEO_cosmic_iron['dose'], yerr=LEO_cosmic_iron['error'], label="LEO Cosmic Iron", 
#              linestyle='', capsize=1, color='C4')

LEO_Total = mergeTotalDose([LEO_electron, LEO_trapped_proton, LEO_solar_proton, LEO_cosmic_proton, LEO_cosmic_iron])
plt.errorbar(x, LEO_Total['dose'] * 100/11, yerr=LEO_Total['error'] * 100/11, color=LEOColor, 
             label="LEO 100 Year Dose", linestyle='', capsize=1.5, capthick=2)

# LEO_SD = readSDQ2("/l/triton_work/Spectra/Carrington/LEO/spenvis_sqo.txt")
# plt.plot(LEO_SD['Thickness'], ( LEO_SD['Electrons'] + LEO_SD['Bremsstrahlung'] ), 
#          label="SHIELDOSE-2Q LEO Electrons + Bremsstrahlung", color='C0')
# plt.plot(LEO_SD['Thickness'], LEO_SD['Trapped Protons'], label="SHIELDOSE-2Q LEO Trapped Protons", color='C1')
# plt.plot(LEO_SD['Thickness'], LEO_SD['Total'] * 100/11, label="SHIELDOSE-2Q LEO Total Dose", color=LEOColor)


#### MEO ######
# MEO_electron = totalDose(Path + "MEO-electron/Res")
# MEO_trapped_proton = totalDose(Path + "MEO-trapped-proton/Res")  # Effectively no trapped protons on MEO
# MEO_solar_proton = totalDose(Path + "MEO-solar-proton/Res")  
# MEO_cosmic_proton = totalDose(Path + "MEO-cosmic-proton/Res")
# MEO_cosmic_iron = totalDose(Path + "MEO-cosmic-iron/Res")

# plt.errorbar(x, MEO_electron['dose'], yerr=MEO_electron['error'], label="MEO Electrons", 
#              linestyle='', capsize=1, color='C0')
# plt.errorbar(x, MEO_trapped_proton['dose'], yerr=MEO_trapped_proton['error'], label="MEO Trapped Protons", 
#              linestyle='', capsize=1, color='C1')
# plt.errorbar(x, MEO_solar_proton['dose'], yerr=MEO_solar_proton['error'], label="MEO Solar Protons", 
#              linestyle='', capsize=1, color='C2')
# plt.errorbar(x, MEO_cosmic_proton['dose'], yerr=MEO_cosmic_proton['error'], label="MEO Cosmic Protons", 
#              linestyle='', capsize=1, color='C3')
# plt.errorbar(x, MEO_cosmic_iron['dose'], yerr=MEO_cosmic_iron['error'], label="MEO Cosmic Iron", 
#              linestyle='', capsize=1, color='C4')

# MEO_Total = mergeTotalDose([MEO_electron, MEO_trapped_proton, MEO_solar_proton, MEO_cosmic_proton, MEO_cosmic_iron])
# plt.errorbar(x, MEO_Total['dose'], yerr=MEO_Total['error'], capsize=1, linestyle='', color=MEOColor, 
#              label="MEO Total", capthick=2)

# MEO_SD = readSDQ2("/l/triton_work/Spectra/Carrington/MEO/spenvis_sqo.txt")
# plt.plot(MEO_SD['Thickness'], ( MEO_SD['Electrons'] + MEO_SD['Bremsstrahlung'] ), 
#          label="SHIELDOSE-2Q MEO Electrons + Bremsstrahlung", color='C0')
# plt.plot(MEO_SD['Thickness'], MEO_SD['Trapped Protons'], label="SHIELDOSE-2Q MEO Trapped Protons", color='C1')
# plt.plot(MEO_SD['Thickness'], MEO_SD['Total'], label="SHIELDOSE-2Q MEO Total Dose", color=MEOColor)


#### VAP ######
# VAP_electron = totalDose(Path + "VAP-electron/Res")
# VAP_trapped_proton = totalDose(Path + "VAP-trapped-proton/Res")
# VAP_trapped_proton_integral = totalDose(Path + "VAP-trapped-proton-integral/Res")
# VAP_solar_proton = totalDose(Path + "VAP-solar-proton/Res")
# VAP_cosmic_proton = totalDose(Path + "VAP-cosmic-proton/Res")
# VAP_cosmic_rron = totalDose(Path + "VAP-cosmic-iron/Res")

# plt.errorbar(x, VAP_electron['dose'], yerr=VAP_electron['error'], label="VAP Electrons", 
#              linestyle='', capsize=1, color='C0')
# plt.errorbar(x, VAP_trapped_proton['dose'], yerr=VAP_trapped_proton['error'], label="VAP Trapped Protons", 
            #  linestyle='', capsize=1, color='C1')
# plt.errorbar(x, VAP_trapped_proton_integral['dose'], yerr=VAP_trapped_proton_integral['error'], label="VAP Trapped Protons Integral", 
            #  linestyle='--', capsize=1, color='black')
# plt.errorbar(x, VAP_solar_proton['dose'], yerr=VAP_solar_proton['error'], label="VAP Solar Protons", 
#              linestyle='', capsize=1, color='C2')
# plt.errorbar(x, VAP_cosmic_proton['dose'], yerr=VAP_cosmic_proton['error'], label="VAP Cosmic Protons", 
#              linestyle='', capsize=1, color='C3')
# plt.errorbar(x, VAP_cosmic_rron['dose'], yerr=VAP_cosmic_rron['error'], label="VAP Cosmic Iron", 
#              linestyle='', capsize=1, color='C4')

# VAP_Total = mergeTotalDose([VAP_electron, VAP_trapped_proton, VAP_solar_proton, VAP_cosmic_proton, VAP_cosmic_rron])
# plt.errorbar(x, VAP_Total['dose'] /11, yerr=VAP_Total['error'] /11, label="VAP 1 Year Dose\n(GRAS)", color=VAPColor, linestyle='', capsize=1.5, capthick=2)

# VAP_SD = readSDQ2("/l/triton_work/Spectra/Carrington/VAP/spenvis_sqo.txt")
# plt.plot(VAP_SD['Thickness'], ( VAP_SD['Electrons'] + VAP_SD['Bremsstrahlung'] ), label="SHIELDOSE-2Q VAP Electrons + Bremsstrahlung", color='C0')
# plt.plot(VAP_SD['Thickness'], VAP_SD['Trapped Protons'], label="SHIELDOSE-2Q VAP Trapped Protons", color='C1')
# plt.plot(VAP_SD['Thickness'], VAP_SD['Total'] /11, '.', label="VAP 1 Year Dose\n(SHIELDOSE-2Q)", color=VAPColor)


# #### GEO ######
GEO_electron = totalDose(Path + "GEO-electron/Res") # Similar but slightly less than VAP-electron, therefore not interesting for shielding
GEO_trapped_proton = totalDose(Path + "GEO-trapped-proton/Res") # The trapped proton spectrum on GEO ends at 8MeV, therefore not intetrsting for shielding
GEO_solar_proton = totalDose(Path + "GEO-solar-proton/Res")
GEO_cosmic_proton = totalDose(Path + "GEO-cosmic-proton/Res") # Cosmics have super low fluxes --> negligible TID
GEO_cosmic_iron = totalDose(Path + "GEO-cosmic-iron/Res")     # Cosmics have super low fluxes --> negligible TID

# plt.errorbar(x, GEO_electron['dose'], yerr=GEO_electron['error'], label="GEO Electrons", linestyle='', capsize=1, color='C0')
# # plt.errorbar(x, GEO_trapped_proton['dose'], yerr=GEO_trapped_proton['error'], label="GEO Trapped Protons", linestyle='', capsize=1, color='C1')
# plt.errorbar(x, GEO_solar_proton['dose'], yerr=GEO_solar_proton['error'], label="GEO Solar Protons", linestyle='', capsize=1, color='C2')
# plt.errorbar(x, GEO_cosmic_proton['dose'], yerr=GEO_cosmic_proton['error'], label="GEO Cosmic Protons", linestyle='', capsize=1, color='C3')
# plt.errorbar(x, GEO_cosmic_iron['dose'], yerr=GEO_cosmic_iron['error'], label="GEO Cosmic Iron", linestyle='', capsize=1, color='C4')

GEO_Total = mergeTotalDose([GEO_electron, GEO_trapped_proton, GEO_solar_proton, GEO_cosmic_proton, GEO_cosmic_iron])
plt.errorbar(x, GEO_Total['dose'] * 10/11, yerr=GEO_Total['error'] * 10/11, color=GEOColor, label="GEO 10 year Dose", linestyle='', capsize=1.5, capthick=2)

# GEO_SD = readSDQ2("/l/triton_work/Spectra/Carrington/GEO/spenvis_sqo.txt")
# plt.plot(GEO_SD['Thickness'], ( GEO_SD['Electrons'] + GEO_SD['Bremsstrahlung'] ), label="SHIELDOSE-2Q GEO Electrons + Bremsstrahlung", color='C0')
# plt.plot(GEO_SD['Thickness'], GEO_SD['Trapped Protons'], label="SHIELDOSE-2Q GEO Trapped Protons", color='C1')
# plt.plot(GEO_SD['Thickness'], GEO_SD['Total'] * 10/11, label="SHIELDOSE-2Q GEO Total Dose", color=GEOColor)


######################## Plot formatting ############################

plt.title("Total Ionising Dose Behind Shielding")
plt.xlabel("Aluminium Shielding Thickness [mm]")
plt.ylabel("Total Ionising Dose [krad]")
plt.yscale("log")
# plt.xscale("log")
# Increase the number of ticks on the y-axis
# plt.yticks(np.logspace(-4, 4, num=9))

plt.xlim(0.5, 10)
plt.ylim(1e-1, 1e+5)
plt.grid(which='both')


################### Re-ordering the legend labels ############################
# Store the handles and labels when plotting
handles, labels = plt.gca().get_legend_handles_labels()

# Reorder the handles and labels to your desired order
# Example: Suppose you want to reorder them based on some criteria or the original plotting order
order = [0, 1, 6, 5, 2, 4, 3]  # PROTON Plot specify the desired order of indices 
# order = [0, 1, 3, 4, 6, 5, 2]  # ELECTRON specify the desired order of indices

# Apply the reordering
handles = [handles[i] for i in order]
labels = [labels[i] for i in order]

# Then pass them to plt.legend()
plt.legend(handles, labels)

# External legend
# plt.legend()


plt.savefig(Path + "CompareGRASShieldose.pdf", format='pdf', bbox_inches="tight")
# plt.show()
