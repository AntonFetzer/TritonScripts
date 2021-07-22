import matplotlib.pyplot as plt
from ReadG4root import readG4root
import numpy as np

# ---------------------- Deprecated -------- DONT USE -----------------------

# ---------------------- Deprecated -------- DONT USE -----------------------

FolderName = "MulasSolarProton2e9Full"
NORM_Spec = 3.798797E+11
N = 2000000000

Data = np.zeros((5, 5), dtype=np.double)  # Data[0] = Thickness, Data[1] = Dose in MeV, Data[2] = Dose in krad
Data[0] = [1, 2, 4, 8, 16]

# Normalisation factor to get from the total recorded dose in MeV in the simulation to the real world dose in krad
SiliconDensity = 2.33  # g/cm3
SiliconThick = 0.05  # cm
MassPerArea = SiliconDensity * SiliconThick / 1000  # kg/cm2
# NORM_Spec = 8.380532E+14  # /cm2
NORM_ANG = 0.25
JoulePerMeV = 1.6E-13  # J/MeV
KradPerGy = 0.1  # krad/Gy
A = NORM_Spec * NORM_ANG / (N * MassPerArea)  # 1/kg
print(A)
C = A * JoulePerMeV * KradPerGy
print(C)

i = 0
for Thick in Data[0]:
    RootData = readG4root("../MULASS/" + FolderName + "/out/" + str(int(Thick)) + ".root")
    # print(Dose)
    Data[1, i] = sum(RootData[1])  # SiVol_Edep_MeV
    Data[2, i] = Data[1, i] * C  # SiVol_Edep_krad
    Data[3, i] = sum(RootData[2])  # SiVol_Esec_MeV
    Data[4, i] = Data[2, i] * C  # SiVol_Esec_krad
    print("Edep for", Thick, "mm of Aluminium Shielding is", str(int(Data[1, i])), "MeV")
    print("Edep for", Thick, "mm of Aluminium Shielding is", str(int(Data[2, i])), "krad")
    print("Esec for", Thick, "mm of Aluminium Shielding is", str(int(Data[3, i])), "MeV")
    print("Esec for", Thick, "mm of Aluminium Shielding is", str(int(Data[4, i])), "krad")
    i += 1


plt.plot(Data[0], Data[2], 'o', label=("G4" + FolderName + " Spectrum"))
plt.plot(Data[0], Data[4], 'o', label="G4" + FolderName + " Spectrum Secondary")
####### Plot 10kRad line #########
CriticalDose = [10] * Data[0].shape[0]
plt.plot(Data[0], CriticalDose, color='k', linewidth=2, label='Critical Dose of 10 krad')
plt.yscale("log")
plt.grid(which='both')
plt.title("G4 total dose")
plt.xlabel("Aluminium Absorber Thickness [mm]")
plt.ylabel("Dose in Si [krad]")
plt.legend()
plt.savefig("/out/" + FolderName + "TotalDoseCurve.eps", format='eps')

np.savetxt("/out/" + FolderName + "DoseTable.txt", Data)

# srun --mem=50G --time=00:15:00 python G4TotalDosetoCSV.py
