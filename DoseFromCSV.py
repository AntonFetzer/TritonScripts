import matplotlib.pyplot as plt
import numpy as np
from natsort import natsorted
import os
import csv
from Dependencies.MeVtokRad_2D import MeVtokRad_2D
from Dependencies.ReadSD2Q import readSDQ2

# -------------- Change these inputs --------------------------
Folder = "MulasCosmicProton1e8Full"
NORM_FACTOR_SPECTRUM = 2.024537E+07
Npart = 1e8
# ------------------------------------------------------------

Path = "/home/anton/Desktop/triton_work/MULASS/"

# Get list of all root files in that folder
CSVFiles = [f for f in os.listdir(Path + Folder + "/CSV/") if f.endswith('.csv')]
CSVFiles = natsorted(CSVFiles)

ThickList = [int(f.split(".")[0]) for f in CSVFiles]

Content = []  # 2D List to store the whole CSV file
i = 0
for file in CSVFiles:
    print(Path + Folder + "/CSV/" + file)
    Content.append([])
    with open(Path + Folder + "/CSV/" + file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            Content[i].append(row)  # Dump the whole CSV in the 2D list
    i += 1

Edep = np.zeros(len(ThickList))    # Edep 1D vector for all Thicknesses
EdepStd = np.zeros(len(ThickList))  # Has to be ndarray becasue otherwise MeVtokRad_2D does not work
Esec = np.zeros(len(ThickList))
EsecStd = np.zeros(len(ThickList))

for i in range(len(ThickList)):
    Edep[i] = Content[i][2][1]  # From ith CSV file read 3rd row and 2nd collumn
    EdepStd[i] = Content[i][2][2]
    Esec[i] = Content[i][3][1]
    EsecStd[i] = Content[i][3][2]

EdepkRads = MeVtokRad_2D(Edep, NORM_FACTOR_SPECTRUM, Npart)
EdepStdkRads = MeVtokRad_2D(EdepStd, NORM_FACTOR_SPECTRUM, Npart)
EseckRads = MeVtokRad_2D(Esec, NORM_FACTOR_SPECTRUM, Npart)
EsecStdkRads = MeVtokRad_2D(EsecStd, NORM_FACTOR_SPECTRUM, Npart)


# ------------ Plot 10kRad line -------------------------
CriticalDose = [10] * len(ThickList)
plt.plot(ThickList, CriticalDose, color='k', linewidth=2, label='Critical Dose of 10 krad')
# ------------- Import and Plot SHIELDOSE Data --------------------------
SDData = readSDQ2("Dependencies/spenvis_sqo.txt")
# SDDataCollumns = ['Aluminium Thickness', 'Total Dose', 'Electrons', 'Bremsstrahlung', 'Protons']
plt.plot(SDData[:, 0], SDData[:, 4] / 1000, label="SHIELDOSE-2Q trapped Protons")
plt.plot(SDData[:, 0], (SDData[:, 2] + SDData[:, 3]) / 1000, label="SHIELDOSE-2Q trapped Electrons + Bremsstrahlung")

# ------------- Plot Edep and Esec with errorbars ----------------------------------------------------------------------
plt.errorbar(ThickList, EdepkRads, EdepStdkRads, fmt='+', capsize=5, markersize=12, label="Total Dose")
plt.errorbar(ThickList, EseckRads, EsecStdkRads, fmt='+', capsize=5, markersize=12, label="Secondary dose")

plt.yscale("log")
plt.xlim([0, 17])
plt.grid(which='both')
plt.title(Folder)
plt.xlabel("Aluminium Absorber Thickness [mm]")
plt.ylabel("Dose in Si [krad]")
plt.legend()

plt.savefig(Path + Folder + "/Img/" + Folder + "TotalDoseCurve.eps", format='eps')

CSVfile = open(Path + Folder + "/" + Folder + "_DoseTable.txt", 'w')
CSVfile.write("Al Shield Thickness, Edep, EdepStd, Esec, EsecStd\n")
for i in range(len(ThickList)):
    CSVfile.write(','.join([str(ThickList[i]), str(EdepkRads[i]), str(EdepStdkRads[i]), str(EseckRads[i]), str(EsecStdkRads[i]), "\n"]))
CSVfile.close()


# plt.show()
# srun --mem=50G --time=00:15:00 python DoseFromCSV.py
