import numpy as np
import matplotlib.pyplot as plt
from Dependencies.MeVtokRad_2D import MeVtokRad_2D

Path = "/home/anton/Desktop/triton_work/Gradient/2Material1-5gcm2/PE-Pb/csv/"

Shield = "PE-Pb"
MatA = "Polyethylene"
MatB = "Lead"
a = "pe"
b = "pb"
A = "PE"
B = "Pb"

Ymax = 14

ElecFileA = "gradient-" + a + "-" + b + "-2e9electron.txt"
ProtFileA = "gradient-" + a + "-" + b + "-2e7proton.txt"
ElecFileB = "gradient-" + b + "-" + a + "-2e9electron.txt"
ProtFileB = "gradient-" + b + "-" + a + "-2e7proton.txt"

NORM_FACTOR_SPECTRUM_Elec = 5.886798E+14  # Elec500keV
NORM_FACTOR_SPECTRUM_Prot = 3.381390E+11  # Prots10MeV
Npart_Elec = 2e9 / 99
Npart_Prot = 2e7 / 99


ElecAdata = np.genfromtxt(Path + ElecFileA, delimiter=',', dtype=None, encoding='ASCII')
ElecBdata = np.genfromtxt(Path + ElecFileB, delimiter=',', dtype=None, encoding='ASCII')
ProtAdata = np.genfromtxt(Path + ProtFileA, delimiter=',', dtype=None, encoding='ASCII')
ProtBdata = np.genfromtxt(Path + ProtFileB, delimiter=',', dtype=None, encoding='ASCII')
TotalEdepA = np.zeros(99)
TotalEdepB = np.zeros(99)

# ElecB = np.flip(ElecB)
# ProtB = np.flip(ProtB)
# TotalEdepB = np.flip(TotalEdepB)

ElecA = np.zeros(99)
ElecB = np.zeros(99)
ProtA = np.zeros(99)
ProtB = np.zeros(99)

ElecAErr = np.zeros(99)
ElecBErr = np.zeros(99)
ProtAErr = np.zeros(99)
ProtBErr = np.zeros(99)

for i in range(99):

    ElecAErr[i] = MeVtokRad_2D(ElecAdata[i][2], NORM_FACTOR_SPECTRUM_Elec, Npart_Elec)
    ElecBErr[i] = MeVtokRad_2D(ElecBdata[i][2], NORM_FACTOR_SPECTRUM_Elec, Npart_Elec)
    ProtAErr[i] = MeVtokRad_2D(ProtAdata[i][2], NORM_FACTOR_SPECTRUM_Prot, Npart_Prot)
    ProtBErr[i] = MeVtokRad_2D(ProtBdata[i][2], NORM_FACTOR_SPECTRUM_Prot, Npart_Prot)

    ElecA[i] = MeVtokRad_2D(ElecAdata[i][1], NORM_FACTOR_SPECTRUM_Elec, Npart_Elec)
    ElecB[i] = MeVtokRad_2D(ElecBdata[i][1], NORM_FACTOR_SPECTRUM_Elec, Npart_Elec)
    ProtA[i] = MeVtokRad_2D(ProtAdata[i][1], NORM_FACTOR_SPECTRUM_Prot, Npart_Prot)
    ProtB[i] = MeVtokRad_2D(ProtBdata[i][1], NORM_FACTOR_SPECTRUM_Prot, Npart_Prot)

    TotalEdepA[i] = ElecA[i] + ProtA[i]
    TotalEdepB[i] = ElecB[i] + ProtB[i]


x = np.linspace(1, 99, num=99, dtype=int)

plt.figure(1)
plt.plot(x, ElecA, '.', label="Electrons " + A + " on " + B + " Min=" + str(round(np.min(ElecA), 2)) + " krad at " + str(round(np.argmin(ElecA) + 1)) + " % " + A)
plt.plot(x, ElecB, '.', label="Electrons " + B + " on " + A + " Min=" + str(round(np.min(ElecB), 2)) + " krad at " + str(round(np.argmin(ElecB) + 1)) + " % " + A)

plt.plot(x, ProtA, '.', label="Protons " + A + " on " + B + " Min=" + str(round(np.min(ProtA), 2)) + " krad at " + str(round(np.argmin(ProtA) + 1)) + " % " + A)
plt.plot(x, ProtB, '.', label="Protons " + B + " on " + A + " Min=" + str(round(np.min(ProtB), 2)) + " krad at " + str(round(np.argmin(ProtB) + 1)) + " % " + A)

plt.ylim(0, Ymax)
plt.title(
    "Dose deposited by trapped particles in 0.5 mm Si \n behind 1.5g/cm2 of " + MatA + "-" + MatB + " shielding")  # ---------
plt.xlabel("Percentage of shielding mass in " + MatA + " [%]")
plt.ylabel("Deposited ionising dose [krad]")
plt.grid(which='both')
plt.legend()
# plt.show()
plt.savefig(Path + Shield + "-Gradient.eps", format='eps', bbox_inches="tight")


plt.figure(2)
plt.plot(x, TotalEdepA, '.',
         label=MatA + " on top of " + MatB + " Min=" + str(round(np.min(TotalEdepA), 2)) + " krad at " + str(
             round(np.argmin(TotalEdepA) + 1)) + " % " + A)
plt.plot(x, TotalEdepB, '.',
         label=MatB + " on top of " + MatA + " Min=" + str(round(np.min(TotalEdepB), 2)) + " krad at " + str(
             round(np.argmin(TotalEdepB) + 1)) + " % " + A)

plt.ylim(0, Ymax)
plt.title(
    "Total dose deposited by trapped particles in 0.5 mm Si \n behind 1.5g/cm2 of " + MatA + "-" + MatB + " shielding")  # --------
plt.xlabel("Percentage of shielding mass in " + MatA + " [%]")
plt.ylabel("Deposited ionising dose [krad]")
plt.grid(which='both')
plt.legend()
# plt.show()
plt.savefig(Path + Shield + "-GradientSum.eps", format='eps', bbox_inches="tight")
# plt.savefig(Path + Shield + "-GradientSum.png", format='png', dpi=400)




TotalAmin = np.min(TotalEdepA)
TotalBmin = np.min(TotalEdepB)
TotalAminIndex = np.argmin(TotalEdepA)
TotalBminIndex = np.argmin(TotalEdepB)

ElecAmin = ElecA[TotalAminIndex]
ElecBmin = ElecB[TotalBminIndex]

ElecAminErr = ElecAErr[TotalAminIndex]
ElecBminErr = ElecBErr[TotalBminIndex]

ProtAmin = ProtA[TotalAminIndex]
ProtBmin = ProtB[TotalBminIndex]

ProtAminErr = ProtAErr[TotalAminIndex]
ProtBminErr = ProtBErr[TotalBminIndex]

TotalAminErr = np.sqrt(ElecAminErr * ElecAminErr + ProtAminErr * ProtAminErr)
TotalBminErr = np.sqrt(ElecBminErr * ElecBminErr + ProtBminErr * ProtBminErr)

CSVFile = open(Path + "Results.txt", 'w')
CSVFile.writelines("Material A, Material B, % A, % B, Electron Dose, Electron Error, Proton Dose, Proton Error, Total Dose, Total Error, \n \n")

if TotalAmin < TotalBmin:
    List = (A, B, TotalAminIndex+1, 99-TotalAminIndex, ElecAmin, ElecAminErr, ProtAmin, ProtAminErr, TotalAmin, TotalAminErr)
elif TotalAmin > TotalBmin:
    List = (A, B, TotalBminIndex+1, 99-TotalBminIndex, ElecBmin, ElecBminErr, ProtBmin, ProtBminErr, TotalBmin, TotalBminErr)

String = ', '.join(map(str, List))
print(String)
CSVFile.writelines(String + "\n")
CSVFile.close()
