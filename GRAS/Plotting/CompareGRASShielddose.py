import numpy as np
import matplotlib.pyplot as plt
from GRAS.Read.ReadSD2Q import readSDQ2
from GRAS.Dependencies.TotalKRadGras import totalkRadGras

Path = "/l/triton_work/ShieldingCurves"

Electrons = totalkRadGras(Path + "/MultilayerPaper/AE9-GTO/Res/", "")
Protons = totalkRadGras(Path + "/MultilayerPaper/AP9-GTO/Res/", "")

x = np.linspace(0, 2.5, num=101, dtype=float)

fig1 = plt.figure(1, [4.5, 5.5])
# ------------------------------- Import and Plot SHIELDOSE Data -------------------------------------------------------
ConvertKrad = 1 / 1000  #  1/1000 to convert to krad
mmtogcm2 = 0.27  # 0.27 to convert mm to g/cm2 of aliminium

SDData = readSDQ2("/l/triton_work/Spectra/A9-GTO/Shieldose/spenvis_sqo.txt")
# SDDataCollumns = ['Aluminium Thickness', 'Total Dose', 'Electrons', 'Bremsstrahlung', 'Protons']
plt.plot(SDData[:, 0] * mmtogcm2, ( SDData[:, 2]+SDData[:, 3] ) * ConvertKrad, 'C0--', label="SHIELDOSE-2Q trapped Electrons")
#plt.plot(SDData[:, 0] * mmtogcm2, SDData[:, 3] * ConvertKrad, 'C1-..', label="SHIELDOSE-2Q Bremsstrahlung")
plt.plot(SDData[:, 0] * mmtogcm2, SDData[:, 4] * ConvertKrad, 'C1-', label="SHIELDOSE-2Q trapped Protons")
plt.plot(SDData[:, 0] * mmtogcm2, SDData[:, 1] * ConvertKrad, 'C2:', label="SHIELDOSE-2Q Total trapped particles")

####### Plot 10kRad line #########
CriticalDose = [1] * 101
plt.plot(x, CriticalDose, color='k', linewidth=2, label='1 krad per month')

# ----------------------------------------- Plot electrons500kev ---------------------------------------------------------
plt.errorbar(x, Electrons[0], Electrons[1], fmt='C0 ', capsize=4, label="Geant4 trapped Electrons", elinewidth=1, capthick=1)

# ----------------------------------------- Plot protons10MeV -------------------------------------------------------
plt.errorbar(x, Protons[0], Protons[1], fmt='C1 ', capsize=4, label="Geant4 trapped Protons", elinewidth=1, capthick=1)

Total = Electrons
Total[0] = Electrons[0] + Protons[0]
Total[1] = np.sqrt(Electrons[1] * Electrons[1] + Protons[1] * Protons[1])

# ----------------------------------------- Plot Total -------------------------------------------------------
plt.errorbar(x, Total[0], Total[1], fmt='C2 ', capsize=4, label="Geant4 Total trapped particles", elinewidth=1, capthick=1)



## Solar Energetic Particles

SEP_Paths = [
    "/l/triton_work/ShieldingCurves/MultilayerPaper/SolarGTO-Protons/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/SolarGTO-He/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/SolarGTO-C/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/SolarGTO-O/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/SolarGTO-Ne/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/SolarGTO-Mg/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/SolarGTO-Si/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/SolarGTO-Fe/Res/"
]

SEP_Data = []
for i, path in enumerate(SEP_Paths):
    SEP_Data.append(totalkRadGras(path, ""))

# Calculate total dose and error for SEP particles
Total_SEP_Dose = sum([data[0] for data in SEP_Data])
Total_SEP_Error = np.sqrt(sum([data[1]**2 for data in SEP_Data]))

# ----------------------------------------- Plot Solar Particles -------------------------------------------------------
plt.errorbar(x, Total_SEP_Dose, Total_SEP_Error, fmt='C8. ', capsize=4, label="Geant4 Solar particles", elinewidth=1, capthick=1)


## Cosmic Particles

Cosmic_Paths = [
    "/l/triton_work/ShieldingCurves/MultilayerPaper/ISO-GTO-Protons-mission/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/ISO-GTO-He-mission/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/ISO-GTO-C-mission/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/ISO-GTO-N-mission/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/ISO-GTO-O-mission/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/ISO-GTO-Mg-mission/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/ISO-GTO-Si-mission/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/ISO-GTO-Fe-mission/Res/"
]

Cosmic_Data = []
for i, path in enumerate(Cosmic_Paths):
    Cosmic_Data.append(totalkRadGras(path, ""))

seconds_in_a_month = 60 * 60 * 24 * 30.44  # number of seconds in a month

# Calculate total dose and error for Cosmic particles
Total_Cosmic_Dose = sum([data[0] for data in Cosmic_Data])/ seconds_in_a_month
Total_Cosmic_Error = np.sqrt(sum([data[1]**2 for data in Cosmic_Data])) / seconds_in_a_month
# Create a new dataset for Cosmic with the same structure
Total_Cosmic_Data = [Total_Cosmic_Dose, Total_Cosmic_Error]

# ----------------------------------------- Plot Cosmic Particles -------------------------------------------------------
plt.errorbar(x, Total_Cosmic_Dose, Total_Cosmic_Error, fmt='C7* ', capsize=4, label="Geant4 Cosmic particles", elinewidth=1, capthick=1)



plt.xlim(0.25, 2.5)
plt.ylim(1e-4, 2e2)
plt.yscale("log")

plt.grid(which='both', linestyle='-', linewidth=0.5)
plt.title("Comparison between Geant4 and SHIELDOSE-2Q", x=0.425)
plt.xlabel("Aluminium Shielding Depth [g/cm2]")
plt.ylabel("Ionizing dose in silicon [krad per month]")
plt.legend(loc='lower left', framealpha=0.75)
#plt.savefig("/l/TritonPlots/Paper/ShielddoseComparison.pdf", format='pdf', bbox_inches="tight")
'''


fig2 = plt.figure(2)
plt.title(Folder)
plt.plot(x, Electrons[3], 'C0.', label="GRAS Electrons > " + str(ElecMeV) + " MeV")
plt.plot(x, Protons[3], 'C1.', label="GRAS Protons > " + str(ProtMeV) + " MeV")
plt.legend()
plt.xlabel("Aluminium Shield Thickness [g/cm2]")
plt.ylabel("Number of Non Zero Entries")
plt.grid(which='both')
plt.yscale("log")
#plt.savefig(Path + Folder + "/NonZeroEntries.eps", format='eps', bbox_inches="tight")

fig3 = plt.figure(3)
plt.title(Folder)
plt.plot(x, 100 * Electrons[1] / Electrons[0], 'C0.', label="GRAS Electrons > " + str(ElecMeV) + " MeV")
plt.plot(x, 100 * Protons[1] / Protons[0], 'C1.', label="GRAS Protons > " + str(ProtMeV) + " MeV")
plt.legend()
plt.xlabel("Aluminium Shield Thickness [g/cm2]")
plt.ylabel("Relative Error in %")
plt.grid(which='both')
#plt.yscale("log")
#plt.savefig(Path + Folder + "/RelativeError.eps", format='eps', bbox_inches="tight")
'''
plt.show()
