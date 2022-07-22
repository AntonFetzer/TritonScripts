import os
import numpy as np
import matplotlib.pyplot as plt
from ReadGRASCSV import readGrasCsv
from MeVtokRad2DGras import MeVtokRad_2D
from TotalKRadGras import totalkRadGras

Path = "/home/anton/Desktop/triton_work/2MatPhys/DONE"
ShieldingDepth = "1.5"  # g/cm2

MatA = "Aluminium"
MatB = "Lead"
a = "al"
b = "pb"
A = "Al"
B = "Pb"
#Shield = A + "-" + B
Shield = "FTFP_BERT"
NumTiles = 99

Ymax = 15  # Max kRad shown in plots so that every plot has the same scale

Path += "/" + Shield + "/"
ElecA = totalkRadGras(Path + "Res/", "ElectronsA") * NumTiles
ProtA = totalkRadGras(Path + "Res/", "ProtonsA") * NumTiles

ElecB = totalkRadGras(Path + "Res/", "ElectronsB") * NumTiles
ProtB = totalkRadGras(Path + "Res/", "ProtonsB") * NumTiles

x = np.linspace(1, 99, num=99, dtype=int)

plt.figure(1)
plt.errorbar(x, ElecA[0], ElecA[1], fmt=' ', capsize=2,
             label="Electrons " + A + " on " + B + " Min=" + str(round(np.min(ElecA[0]), 2)) + " kRad at " + str(
                 round(np.argmin(ElecA[0]) + 1)) + " % " + A)
plt.errorbar(x, ElecB[0], ElecB[1], fmt=' ', capsize=2,
             label="Electrons " + B + " on " + A + " Min=" + str(round(np.min(ElecB[0]), 2)) + " kRad at " + str(
                 round(np.argmin(ElecB[0]) + 1)) + " % " + A)

plt.errorbar(x, ProtA[0], ProtA[1], fmt=' ', capsize=2,
             label="Protons " + A + " on " + B + " Min=" + str(round(np.min(ProtA[0]), 2)) + " kRad at " + str(
                 round(np.argmin(ProtA[0]) + 1)) + " % " + A)
plt.errorbar(x, ProtB[0], ProtB[1], fmt=' ', capsize=2,
             label="Protons " + B + " on " + A + " Min=" + str(round(np.min(ProtB[0]), 2)) + " kRad at " + str(
                 round(np.argmin(ProtB[0]) + 1)) + " % " + A)

plt.ylim(0, Ymax)
plt.title("Dose deposited by trapped particles in 0.5 mm Si \n behind " + ShieldingDepth + "/cm2 of " + MatA + "-" + MatB + " shielding")
plt.xlabel("Percentage of shielding mass in " + MatA + " [%]")
plt.ylabel("Deposited ionising dose [kRad]")
plt.grid(which='both')
plt.legend()
# plt.show()
plt.savefig(Path + Shield + "-Gradient.eps", format='eps', bbox_inches="tight")

TotalA = ElecA
TotalA[0] = ElecA[0] + ProtA[0]
TotalA[1] = np.sqrt(ElecA[1] * ElecA[1] + ProtA[1] * ProtA[1])

TotalB = ElecB
TotalB[0] = ElecB[0] + ProtB[0]
TotalB[1] = np.sqrt(ElecB[1] * ElecB[1] + ProtB[1] * ProtB[1])

plt.figure(2)
plt.errorbar(x, TotalA[0], TotalA[1], fmt=' ', capsize=2,
             label=MatA + " on top of " + MatB + " Min=" + str(round(np.min(TotalA[0]), 2)) + " krad at " + str(
             round(np.argmin(TotalA[0]) + 1)) + " % " + A)
plt.errorbar(x, TotalB[0], TotalB[1], fmt=' ', capsize=2,
             label=MatB + " on top of " + MatA + " Min=" + str(round(np.min(TotalB[0]), 2)) + " krad at " + str(
             round(np.argmin(TotalB[0]) + 1)) + " % " + A)

plt.ylim(0, Ymax)
plt.title(
    "Total dose deposited by trapped particles in 0.5 mm Si \n behind " + ShieldingDepth + "/cm2 of " + MatA + "-" + MatB + " shielding")  # --------
plt.xlabel("Percentage of shielding mass in " + MatA + " [%]")
plt.ylabel("Deposited ionising dose [krad]")
plt.grid(which='both')
plt.legend()
# plt.show()
plt.savefig(Path + Shield + "-GradientSum.eps", format='eps', bbox_inches="tight")


TotalAmin = np.min(TotalA[0])
TotalBmin = np.min(TotalB[0])
TotalAminIndex = np.argmin(TotalA[0])
TotalBminIndex = np.argmin(TotalB[0])

ElecAmin = ElecA[0][TotalAminIndex]
ElecBmin = ElecB[0][TotalBminIndex]

ElecAminErr = ElecA[1][TotalAminIndex]
ElecBminErr = ElecB[1][TotalBminIndex]

ProtAmin = ProtA[0][TotalAminIndex]
ProtBmin = ProtB[0][TotalBminIndex]

ProtAminErr = ProtA[1][TotalAminIndex]
ProtBminErr = ProtB[1][TotalBminIndex]

TotalAminErr = np.sqrt(ElecAminErr * ElecAminErr + ProtAminErr * ProtAminErr)
TotalBminErr = np.sqrt(ElecBminErr * ElecBminErr + ProtBminErr * ProtBminErr)

CSVFile = open(Path + "Results.txt", 'w')
CSVFile.writelines("Material A, Material B, % A, % B, Electron Dose, Electron Error, Proton Dose, Proton Error, Total Dose, Total Error, \n")

print("TotalAmin:", TotalAmin, "TotalBmin:", TotalBmin)

if TotalAmin < TotalBmin:
    List = (A, B, TotalAminIndex+1, 99-TotalAminIndex, ElecAmin, ElecAminErr, ProtAmin, ProtAminErr, TotalAmin, TotalAminErr)
elif TotalAmin > TotalBmin:
    List = (B, A, 99-TotalBminIndex, TotalBminIndex+1, ElecBmin, ElecBminErr, ProtBmin, ProtBminErr, TotalBmin, TotalBminErr)

String = ', '.join(map(str, List))
print(String)
CSVFile.writelines(String + "\n")
CSVFile.close()
