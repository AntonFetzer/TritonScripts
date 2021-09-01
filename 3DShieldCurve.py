import numpy as np
import matplotlib.pyplot as plt
from ReadSD2Q import readSDQ2
import os
import csv

Path = "/home/anton/Desktop/triton_work/3D/MultiChipTest/AluVault1cm2Hole/csv/"

# Get list of all root files in that folder
CSVFiles = [f for f in os.listdir(Path) if f.endswith('.csv')]
print(CSVFiles)

ThickList = [1, 2, 4, 8, 16]

CSVFilesContent = []  # 2D List to store the whole CSV file

i = 0
for file in CSVFiles:
    # print(Path + file)
    CSVFilesContent.append([])
    with open(Path + file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            CSVFilesContent[i].append(row)  # Dump the whole CSV in the 2D list
    i += 1

# print(CSVFilesContent[0][0][0][:])  #

ElectronFull = np.zeros((4, len(ThickList), 4))  # [Sivol, Thickness, Variable]
SolarProton = np.zeros((4, len(ThickList), 4))
TrapProton = np.zeros((4, len(ThickList), 4))
CosmicProton = np.zeros((4, len(ThickList), 4))

for j in range(4):
    for i in range(len(ThickList)):
        SolarProton[j][i] = CSVFilesContent[0][2 + j + i * 7][1:5]
        CosmicProton[j][i] = CSVFilesContent[1][2 + j + i * 7][1:5]
        ElectronFull[j][i] = CSVFilesContent[2][2 + j + i * 7][1:5]
        TrapProton[j][i] = CSVFilesContent[3][2 + j + i * 7][1:5]

print(ElectronFull[0][:, 0])

Total = ElectronFull + SolarProton + TrapProton + CosmicProton

# ------------------------------- Import and Plot SHIELDOSE Data -------------------------------------------------------
SDData = readSDQ2("spenvis_sqo.txt")
# SDDataCollumns = ['Aluminium Thickness', 'Total Dose', 'Electrons', 'Bremsstrahlung', 'Protons']
plt.plot(SDData[:, 0], SDData[:, 4] / 1000, label="SHIELDOSE-2Q trapped Protons")
plt.plot(SDData[:, 0], (SDData[:, 2] + SDData[:, 3]) / 1000, label="SHIELDOSE-2Q trapped Electrons")
plt.plot(SDData[:, 0], SDData[:, 1] / 1000, label="SHIELDOSE-2Q Total trapped Dose")

# ----------------------------------------- Plot Total Dose ------------------------------------------------------------
plt.errorbar(ThickList, Total[0][:, 0], Total[0][:, 1], fmt='+', capsize=5, markersize=10, label="Geant-4 Total Dose")

# ----------------------------------------- Plot ElectronsFull -------------------------------------------------------------
plt.errorbar(ThickList, ElectronFull[0][:, 0], ElectronFull[0][:, 1], fmt='+', capsize=5, markersize=10, label="Geant-4 Trapped Electrons full spectrum")

# ----------------------------------------- Plot Trapped Protons -------------------------------------------------------
plt.errorbar(ThickList, TrapProton[0][:, 0], TrapProton[0][:, 1], fmt='+', capsize=5, markersize=10, label="Geant-4 Trapped Protons")

# ----------------------------------------- Plot Solar Protons ---------------------------------------------------------
plt.errorbar(ThickList, SolarProton[0][:, 0], SolarProton[0][:, 1], fmt='+', capsize=5, markersize=10, label="Geant-4 Solar Protons")

# ----------------------------------------- Plot Cosmic Protons --------------------------------------------------------
plt.errorbar(ThickList, CosmicProton[0][:, 0], CosmicProton[0][:, 1], fmt='+', capsize=5, markersize=10, label="Geant-4 Cosmic Protons")

####### Plot 10kRad line #########
CriticalDose = [10] * SDData.shape[0]
plt.plot(SDData[:, 0], CriticalDose, color='k', linewidth=2, label='Critical Dose of 10 krad')

plt.xlim([0, 17])
plt.ylim([1e-4, 1e7])
plt.yscale("log")
# plt.xscale("log")
plt.grid(which='major')
plt.title("Trapped particle spectra shielded by aluminium vault")
plt.xlabel("Aluminium Vault Thickness [mm]")
plt.ylabel("Dose in silicon detector [krad]")
plt.legend()
plt.show()
