import numpy as np
import matplotlib.pyplot as plt
from Dependencies.ReadSD2Q import readSDQ2
import os
import csv

Path = "/l/triton_work/3D/AluVaultGap/csv/"

# Get list of all csv files in that folder
CSVFiles = [f for f in os.listdir(Path) if f.endswith('.csv')]
print(CSVFiles)

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

Data = np.zeros((3, 4, 5, 4))  # [Particle, Sivol, Config, Variable]
# 0   Proton      0   small   Edep
# 1   Electrons   1   back    EdepStd
# 2   Sum         2   large   Esec
# 3               3   right   EsecStd
# 4               4    NO
for particle in range(2):
    for sivol in range(4):
        for config in range(5):
            Data[particle][sivol][config] = CSVFilesContent[particle][2 + sivol + config * 7][1:5]

EdepList = np.zeros(4)
EsecList = np.zeros(4)
EdepErrList = np.zeros(4)

# Particle #2 is the total dose
for sivol in range(4):
    for config in range(5):
        for particle in range(2):
            EdepList[particle] = Data[particle][sivol][config][0]
            EdepErrList[particle] = Data[particle][sivol][config][1]
        Data[2][sivol][config][0] = np.sum(EdepList)
        Data[2][sivol][config][1] = np.sqrt(np.sum(np.square(EdepErrList)))

print(Data[2, :, 4, 0], Data[2, :, 4, 1])

plt.hlines(10, -1, 5, colors='k', linewidth=3, label='Critical Dose of 10 krad', zorder=2)

for i in range(4):
    plt.vlines(0.5+i, ymin=0, ymax=150, colors='gray', zorder=1)

plt.grid(axis='y', color='gray')

####### Plot 10kRad line #########
CriticalDose = [10] * 5
plt.plot(["No Gap", "Small Gap \n in the front", "Large Gap \n in the front", "Side Gap", "Back Gap"], CriticalDose, 'k', zorder=2)

x = [-0.3, -0.1, +0.1, +0.3]
plt.bar(x, Data[1, :, 4, 0], yerr=Data[1, :, 4, 1], color='blue', width=0.18, align='center', ecolor='black', capsize=6, label="Trapped Electrons", zorder=3)
plt.bar(x, Data[0, :, 4, 0], yerr=Data[0, :, 4, 1], color='red', width=0.18, align='center', ecolor='black', capsize=6, label="Trapped Protons", zorder=3)
# plt.errorbar(x, Data[0, :, 4, 0], Data[0, :, 4, 1], fmt='ro', capsize=5, markersize=5, label="Geant-4 Trapped Protons")
# plt.errorbar(x, Data[1, :, 4, 0], Data[1, :, 4, 1], fmt='bD', capsize=5, markersize=5, label="Geant-4 Trapped Electrons")
# plt.errorbar(x, Data[2, :, 4, 0], Data[2, :, 4, 1], fmt='ks', capsize=5, markersize=5, label="Geant-4 Total Dose")

i = 1
for config in [0, 2, 3, 1]:
    x = [i - 0.3, i - 0.1, i + 0.1, i + 0.3]
    plt.bar(x, Data[1, :, config, 0], yerr=Data[1, :, config, 1], color='blue', width=0.18, align='center', ecolor='black', capsize=6, zorder=3)
    plt.bar(x, Data[0, :, config, 0], yerr=Data[0, :, config, 1], color='red', width=0.18, align='center', ecolor='black', capsize=6, zorder=3)
    i += 1
# plt.errorbar(x, Data[0, :, config, 0], Data[0, :, config, 1], fmt='ro', capsize=5, markersize=5)
# plt.errorbar(x, Data[1, :, config, 0], Data[1, :, config, 1], fmt='bD', capsize=5, markersize=5)
# plt.errorbar(x, Data[2, :, config, 0], Data[2, :, config, 1], fmt='ks', capsize=5, markersize=5)





plt.xlim([-0.5, 4.5])
#plt.ylim([2, 300])
#plt.yscale("log")
# plt.xscale("log")

plt.title("Four detectors shielded by 6mm aluminium vault with \n gaps against GTO trapped electron and proton spectra")
plt.xlabel("Vault Configuration")
plt.ylabel("Dose in silicon detector [krad]")
plt.legend()
#plt.show()
plt.savefig(Path + "VaultGaps.eps", format='eps', bbox_inches="tight")
