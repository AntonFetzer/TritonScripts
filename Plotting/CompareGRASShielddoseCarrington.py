import numpy as np
import matplotlib.pyplot as plt
from GRAS.Read.ReadSD2Q import readSDQ2
from GRAS.Dependencies.TotalDose import totalDose
from GRAS.Dependencies.MergeTotalDose import mergeTotalDose


##################################### Reading in GRAS data ########################################

Path = "/l/triton_work/Shielding_Curves/Carrington/"

#### Carrington SEP fluxes ######
Ex = totalDose(Path + "Carrington-SEP-Expected-Int-With0/Res")
Minus = totalDose(Path + "Carrington-SEP-Minus2Sigma-Int-With0/Res")
Plus = totalDose(Path + "Carrington-SEP-Plus2Sigma-Int-With0/Res")
# Alternative SEP spectra
# Carrington_SEP_Expected_Diff = totalDose(Path + "Carrington-SEP-Expected-Diff/Res")
# Carrington_SEP_Expected_Int = totalDose(Path + "Carrington-SEP-Expected-Int/Res")

#### Carrington Electron fluxes ######
CarringtonElectronINTEGRALPowTabelated = totalDose(Path + "CarringtonElectronINTEGRALPowTabelated/Res")
# Alternative Electron spectra
# CarringtonElectronDiffPowTabelated = totalDose(Path + "CarringtonElectronDiffPowTabelated/Res")
# CarringtonElectronDiffPow = totalDose(Path + "CarringtonElectronDiffPow/Res")

#### ISS ######
ISS_AE9_mission = totalDose(Path + "ISS-AE9-mission/Res")
ISS_AP9_mission = totalDose(Path + "ISS-AP9-mission/Res")

ISS_SolarProton_mission = totalDose(Path + "ISS-SolarProton-mission/Res")

ISS_A9_Total = mergeTotalDose([ISS_AE9_mission, ISS_AP9_mission, ISS_SolarProton_mission])

#### VAB ######
VAB_AE9_mission = totalDose(Path + "VAB-AE9-mission/Res")
VAB_AP9_mission = totalDose(Path + "VAB-AP9-mission/Res")

VAB_SolarProton_mission = totalDose(Path + "VAB-SolarProton-mission/Res")

VAB_A9_Total = mergeTotalDose([VAB_AE9_mission, VAB_AP9_mission, VAB_SolarProton_mission])


# #### GEO ######
GEO_AE9_mission = totalDose(Path + "GEO-AE9-mission/Res") # Similar but slightly less than VAB-AE9-mission, therefore not interesting for shielding
GEO_AP9_mission = totalDose(Path + "GEO-AP9-mission/Res") # The trapped proton spectrum on GEO ends at 8MeV, therefore not intetrsting for shielding

GEO_SolarProton_mission = totalDose(Path + "GEO-SolarProton-mission/Res")
GEO_SolarProton_5minPeakFlux = totalDose(Path + "GEO-SolarProton-5minPeakFlux/Res") # Similar Flux as the carrington protons, but onlf for 5 minutes, this makes the TID negligible

# GEO_CosmicProton_mission = totalDose(Path + "GEO-CosmicProton-mission/Res") # Cosmics have super low fluxes --> negligible TID
# GEO_CosmicIron_mission = totalDose(Path + "GEO-CosmicIron-mission/Res")     # Cosmics have super low fluxes --> negligible TID

GEO_A9_Total = mergeTotalDose([GEO_AE9_mission, GEO_AP9_mission, GEO_SolarProton_mission])




############################# Plotting GRAS data ########################################
plt.figure(0, [4, 6])
# Aluminium Shielding Thickness [mm]
# 101 Tiles from 0 to 10 mm
x = np.linspace(0, 10, num=101, dtype=float, endpoint=True)
# print(x)

# Plot ISS total
plt.errorbar(x, ISS_A9_Total['dose'] * 12 * 100, yerr=ISS_A9_Total['error'] * 12 * 100, capsize=3, linestyle='', color='C8', label="100 year dose on ISS LEO")

# Plot GEO total
plt.errorbar(x, GEO_A9_Total['dose'] * 120, yerr=GEO_A9_Total['error'] * 120, capsize=3, linestyle='', color='C7', label="10 year dose on GEO")

# Plot VAB total
plt.errorbar(x, VAB_A9_Total['dose'], yerr=VAB_A9_Total['error'], capsize=3, linestyle='', color='C3', label="1 month dose\nVan-Allen-Belt Probes")

# Plot Solar Proton 5min Peak Flux
# The flux is in 1/s cm2, therefore we need to multiply by to get the event dose in krad
plt.errorbar(x, GEO_SolarProton_5minPeakFlux['dose'] * 60 * 60 * 24 * 1.23, yerr=GEO_SolarProton_5minPeakFlux['error'] * 60 * 60 * 24 * 1.23, capsize=3, linestyle='', color='C4', label="1 month peak flux on GEO")

####### Plot 1kRad line #########
CriticalDose = [1 for i in x]
plt.plot(x, CriticalDose, color='k', linewidth=2, label='1 krad')
CriticalDose = [10 for i in x]
plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='10 krad')
# CriticalDose = [100 for i in x]
# plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='100 krad')


##### SHIELDOSE #####

# ------------------------------- Import and Plot SHIELDOSE Data -------------------------------------------------------
ConvertKrad = 1 / 1000  #  1/1000 to convert to krad
# mmtogcm2 = 0.27  # 0.27 to convert mm to g/cm2 of aliminium
# SDData = readSDQ2("/l/triton_work/Spectra/Van-Allen-Probes/Shieldose/spenvis_sqo.txt")
# SDDataCollumns = ['Aluminium Thickness', 'Total Dose', 'Electrons', 'Bremsstrahlung', 'Protons']
#plt.plot(SDData[:, 0], ( SDData[:, 2]+SDData[:, 3] ) * ConvertKrad, 'C0--', label="SHIELDOSE-2Q trapped Electrons")
#plt.plot(SDData[:, 0] * mmtogcm2, SDData[:, 3] * ConvertKrad, 'C1-..', label="SHIELDOSE-2Q Bremsstrahlung")
#plt.plot(SDData[:, 0], SDData[:, 4] * ConvertKrad, 'C1-', label="SHIELDOSE-2Q trapped Protons")
# plt.plot(SDData[:, 0], SDData[:, 1] * ConvertKrad * 12, '.', label="1 year dose\nVan-Allen-Belt Probes\n(SHIELDOSE-2Q)", color='C3')



# Plot Carrington SEP fluxes
# The carrington fluxes are per second, but we are interested in the event fluence which corresponds to 1.23 d
PerStoEvent = 1.23 * 24 * 60 * 60
plt.errorbar(x, Ex['dose'] * PerStoEvent, Ex['error'] * PerStoEvent, fmt='', markersize=2, capsize=2, label='Carrington SEP Event Dose', color='blue', linestyle='')
plt.fill_between(x, Plus['dose'] * PerStoEvent, Ex['dose'] * PerStoEvent, color='C1', alpha=0.5, label='Carrington SEP + 2\u03C3')
plt.fill_between(x, Ex['dose'] * PerStoEvent, Minus['dose'] * PerStoEvent, color='C2', alpha=0.5, label='Carrington SEP − 2\u03C3')
# Plot Carrington Electron fluxes
plt.errorbar(x, CarringtonElectronINTEGRALPowTabelated['dose'] * PerStoEvent, yerr=CarringtonElectronINTEGRALPowTabelated['error'] * PerStoEvent, fmt='', markersize=3, capsize=3, label='Carrington Electron Event Dose', color='blue', linestyle='')



plt.title("Ionising Dose Behind Shielding")
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
#order = [0, 1, 4, 5, 2, 6, 7, 3]  # specify the desired order of indices
order = [0, 1, 2, 7, 3, 4, 5, 6, 8]  # specify the desired order of indices

# Apply the reordering
handles = [handles[i] for i in order]
labels = [labels[i] for i in order]

# Then pass them to plt.legend()
plt.legend(handles, labels)
# plt.legend()


plt.savefig(Path + "CompareGRASShieldose.pdf", format='pdf', bbox_inches="tight")
# plt.show()
