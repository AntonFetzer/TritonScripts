import numpy as np
import matplotlib.pyplot as plt
from Dependencies.ReadSD2Q import readSDQ2
from natsort import natsorted
import os
import csv

Path = "/home/anton/Desktop/triton_work/3D/AluVault/csv/"

# Get list of all csv files in that folder
CSVFiles = [f for f in os.listdir(Path) if f.endswith('.csv')]
print(CSVFiles)
CSVFiles = natsorted(CSVFiles)
print(CSVFiles)

ThickList = [0, 1, 2, 4, 8, 16]

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

Data = np.zeros((5, 5, len(ThickList), 4))  # [Particle, Sivol, Thickness, Variable]

for particle in range(4):
    for sivol in range(4):
        for thick in range(len(ThickList)):
            Data[particle][sivol][thick] = CSVFilesContent[particle][2 + sivol + thick * 7][1:5]

EdepList = np.zeros(4)
EsecList = np.zeros(4)

# Sivol #4 average of the Sivols 0 to 3
for particle in range(4):
    for thick in range(6):
        for Sivol in range(4):
            EdepList[Sivol] = Data[particle][Sivol][thick][0]
            EsecList[Sivol] = Data[particle][Sivol][thick][2]
        Data[particle][4][thick][0] = np.mean(EdepList)
        Data[particle][4][thick][1] = np.max(EdepList)-np.min(EdepList)
        Data[particle][4][thick][2] = np.mean(EsecList)
        Data[particle][4][thick][3] = np.max(EsecList)-np.min(EsecList)

EdepErrList = np.zeros(4)

# Particle #4 is the total dose
for thick in range(6):
    for particle in range(4):
        EdepList[particle] = Data[particle][4][thick][0]
        EdepErrList[particle] = Data[particle][4][thick][1]
    Data[4][4][thick][0] = np.sum(EdepList)
    Data[4][4][thick][1] = np.sqrt(np.sum(np.square(EdepErrList)))

print(Data[0][0][:, 0], Data[0][0][:, 1])

# ------------------------------- Import and Plot SHIELDOSE Data -------------------------------------------------------
SDData = readSDQ2("../Dependencies/spenvis_sqo.txt")
# SDDataCollumns = ['Aluminium Thickness', 'Total Dose', 'Electrons', 'Bremsstrahlung', 'Protons']
plt.plot(SDData[:, 0], (SDData[:, 2] + SDData[:, 3]) / 1000, '-.', label="SHIELDOSE-2Q trapped Electrons")
plt.plot(SDData[:, 0], SDData[:, 4] / 1000, '--', label="SHIELDOSE-2Q trapped Protons")
plt.plot(SDData[:, 0], SDData[:, 1] / 1000, '-', label="SHIELDOSE-2Q Total trapped Dose")

# ----------------------------------------- Plot Total Dose ------------------------------------------------------------
plt.errorbar(ThickList, Data[4][4][:, 0], Data[4][4][:, 1], fmt='ks', capsize=5, markersize=5, label="Geant-4 Total Dose")

# ----------------------------------------- Plot ElectronsFull ---------------------------------------------------------
plt.errorbar(ThickList, Data[0][4][:, 0], Data[0][4][:, 1], fmt='bD', capsize=5, markersize=4, label="Geant-4 Trapped Electrons")

# ----------------------------------------- Plot Trapped Protons -------------------------------------------------------
plt.errorbar(ThickList, Data[1][4][:, 0], Data[1][4][:, 1], fmt='rp', capsize=5, markersize=5, label="Geant-4 Trapped Protons")

# ----------------------------------------- Plot Solar Protons ---------------------------------------------------------
plt.errorbar(ThickList, Data[2][4][:, 0], Data[2][4][:, 1], fmt='yo', capsize=5, markersize=5, label="Geant-4 Solar Protons")

# ----------------------------------------- Plot Cosmic Protons --------------------------------------------------------
plt.errorbar(ThickList, Data[3][4][:, 0], Data[3][4][:, 1], fmt='mX', capsize=5, markersize=5, label="Geant-4 Cosmic Protons")


####### Plot 10kRad line #########
CriticalDose = [10] * SDData.shape[0]
plt.plot(SDData[:, 0], CriticalDose, color='k', linewidth=2, label='Critical Dose of 10 krad')

plt.xlim([-0.5, 17])
plt.ylim([1e-4, 1e7])
plt.yscale("log")
# plt.xscale("log")
plt.grid(which='major')
plt.title("GTO particle spectra shielded by 1U aluminium vault")
plt.xlabel("Aluminium Vault Thickness [mm]")
plt.ylabel("Dose in silicon detector [krad]")
plt.legend()
plt.show()
#plt.savefig(Path + "3DShieldCurve.eps", format='eps', bbox_inches="tight")
