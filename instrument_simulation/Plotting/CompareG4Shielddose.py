import numpy as np
import matplotlib.pyplot as plt
from instrument_simulation.Dependencies.ReadSD2Q import readSDQ2
from instrument_simulation.Dependencies.MeVtokRad_2D import MeVtokRad_2D
from natsort import natsorted
import os
import csv

Path = "/home/anton/Desktop/TritonOffline/MulasTest2/"

#
Npart = 2e9
#                          [electrons500kev, electronsfull, protons10mev, protonsfull]
NORM_FACTOR_SPECTRUM_List = [5.886798E+14, 6.159454E+15, 3.381390E+11, 8.003046E+14]

# Get list of all csv files in that folder
CSVFiles = [f for f in os.listdir(Path) if f.endswith('.txt')]
# print(CSVFiles)
CSVFiles = natsorted(CSVFiles)
# print(CSVFiles)

ThickList = [0, 1, 2, 4, 8, 16]

CSVFilesContent = []  # 2D List to store the whole CSV file

i = 0
for file in CSVFiles:
    # print(Path + file)
    with open(Path + file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            CSVFilesContent.append(row)  # Dump the whole CSV in the 2D list
    CSVFilesContent[i][0] = file.split(".")[0]
    i += 1

#print([x[0] for x in CSVFilesContent])

#print(CSVFilesContent)


Data = np.zeros((len(ThickList), 4, 4))  # [Thickness, Spectrum ,Variable]
# Thickness = [0, 1, 2, 4, 8, 16]
# Spectrum = [electrons500kev, electronsfull, protons10mev, protonsfull]
# Variable = [Edep in MeV, Error in MeV, Edep in krad, Error in krad]

for thick in range(len(ThickList)):
    for spec in range(4):
        for var in range(2):
            MeV = float(CSVFilesContent[thick*4 + spec][var+1])
            Data[thick][spec][var] = MeV

            Data[thick][spec][var+2] = MeVtokRad_2D(MeV, NORM_FACTOR_SPECTRUM_List[spec], Npart)

print(Data)

# ------------------------------- Import and Plot SHIELDOSE Data -------------------------------------------------------
SDData = readSDQ2("../Dependencies/spenvis_sqo.txt")
# SDDataCollumns = ['Aluminium Thickness', 'Total Dose', 'Electrons', 'Bremsstrahlung', 'Protons']
plt.plot(SDData[:, 0], (SDData[:, 2] + SDData[:, 3]) / 1000, '-.', label="SHIELDOSE-2Q Electrons")
plt.plot(SDData[:, 0], SDData[:, 4] / 1000, '--', label="SHIELDOSE-2Q Protons")
plt.plot(SDData[:, 0], SDData[:, 1] / 1000, '-', label="SHIELDOSE-2Q Total Dose")

####### Plot 10kRad line #########
CriticalDose = [10] * SDData.shape[0]
plt.plot(SDData[:, 0], CriticalDose, color='k', linewidth=2, label='Critical Dose of 10 krad')

# ----------------------------------------- Plot electrons500kev ---------------------------------------------------------
plt.errorbar(ThickList, Data[:, 0, 2], Data[:, 0, 3], fmt='ks', capsize=7, markersize=6, label="Geant-4 Electrons > 500 keV")

# ----------------------------------------- Plot electronsfull ---------------------------------------------------------
plt.errorbar(ThickList, Data[:, 1, 2], Data[:, 1, 3], fmt='bD', capsize=7, markersize=5, label="Geant-4 Electrons > 40 keV")

# ----------------------------------------- Plot protons10MeV -------------------------------------------------------
plt.errorbar(ThickList, Data[:, 2, 2], Data[:, 2, 3], fmt='rp', capsize=7, markersize=5, label="Geant-4 Protons > 10 MeV")

# ----------------------------------------- Plot protonsfull -------------------------------------------------------
plt.errorbar(ThickList, Data[:, 3, 2], Data[:, 3, 3], fmt='yo', capsize=7, markersize=5, label="Geant-4 Protons > 100 keV")

# plt.xlim([0, 17])
# plt.ylim([0.2, 2e2])
plt.yscale("log")
# plt.xscale("log")
plt.grid(which='major')
plt.title("Trapped particle spectra shielded by aluminium")
plt.xlabel("Aluminium Shield Thickness [mm]")
plt.ylabel("Dose in silicon detector [krad]")
plt.legend()
#plt.show()
plt.savefig(Path + "CompareG4SHieldose", format='eps', bbox_inches="tight")
