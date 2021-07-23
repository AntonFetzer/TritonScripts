import matplotlib.pyplot as plt
import numpy as np
from natsort import natsorted
import os
import csv
from MeVtokRad_2D import MeVtokRad_2D

Folder = "MulasTrapProton2e9AlFull"
Path = "/home/anton/Desktop/triton_work/MULASS/"

CSVFiles = [f for f in os.listdir(Path + Folder + "/CSV/") if f.endswith('.csv')]
CSVFiles = natsorted(CSVFiles)

ThickList = [int(f.split(".")[0]) for f in CSVFiles]

Content = []

i = 0
for file in CSVFiles:
    print(Path + Folder + "/CSV/" + file)
    Content.append([])
    with open(Path + Folder + "/CSV/" + file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            Content[i].append(row)
    i += 1

Edep = np.zeros(len(ThickList), dtype=np.double)
EdepStd = np.zeros(len(ThickList), dtype=np.double)
Esec = np.zeros(len(ThickList), dtype=np.double)
EsecStd = np.zeros(len(ThickList), dtype=np.double)

for i in range(len(ThickList)):
    Edep[i] = Content[i][2][1]
    EdepStd[i] = Content[i][2][2]
    Esec[i] = Content[i][3][1]
    EsecStd[i] = Content[i][3][2]

# print(Edep)
# print(Esec)

EdepkRads = MeVtokRad_2D(Edep, 8.003046E+14, 2e9)
EdepStdkRads = MeVtokRad_2D(EdepStd, 8.003046E+14, 2e9)
EseckRads = MeVtokRad_2D(Esec, 8.003046E+14, 2e9)
EsecStdkRads = MeVtokRad_2D(EsecStd, 8.003046E+14, 2e9)


####### Plot 10kRad line #########
CriticalDose = [10] * len(ThickList)
plt.plot(ThickList, CriticalDose, color='k', linewidth=2, label='Critical Dose of 10 krad')

plt.errorbar(ThickList, EdepkRads, EdepStdkRads, fmt='+', capsize=5,  label=Folder)
plt.errorbar(ThickList, EseckRads, EsecStdkRads, fmt='+', capsize=5,  label=Folder)


# plt.plot(Data[0], Data[4], 'o', label="G4" + FolderName + " Spectrum Secondary")

plt.yscale("log")
plt.grid(which='both')
plt.title("G4 total dose")
plt.xlabel("Aluminium Absorber Thickness [mm]")
plt.ylabel("Dose in Si [krad]")
plt.legend()
# plt.savefig("/out/" + FolderName + "TotalDoseCurve.eps", format='eps')

# np.savetxt("/out/" + FolderName + "DoseTable.txt", Data)

plt.show()
# srun --mem=50G --time=00:15:00 python DoseFromCSV.py
