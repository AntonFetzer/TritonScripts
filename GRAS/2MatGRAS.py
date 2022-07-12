import os
import numpy as np
import matplotlib.pyplot as plt
from ReadGRASCSV import readGrasCsv
from MeVtokRad2DGras import MeVtokRad_2D
from TotalKRadGras import totalkRadGras

Path = "/home/anton/Desktop/triton_work/GRAS-2Mat"
ShieldingDepth = "1.5"  # g/cm2

MatA = "Aluminium"
MatB = "Lead"
a = "al"
b = "pb"
A = "Al"
B = "Pb"
Shield = A + "-" + B
NumTiles = 99

Ymax = 10  # Max kRad shown in plots so that every plot has the same scale

Path += "/" + Shield + "/"
ElecA = totalkRadGras(Path + "/Results/", "Elec") * NumTiles
ProtA = totalkRadGras(Path + "/Results/", "Prot") * NumTiles

x = np.linspace(1, 99, num=99, dtype=int)

plt.figure(1)
plt.errorbar(x, ElecA[0], ElecA[1], fmt=' ', capsize=2,
             label="Electrons " + A + " on " + B + " Min=" + str(round(np.min(ElecA[0]), 2)) + " kRad at " + str(
                 round(np.argmin(ElecA[0]) + 1)) + " % " + A)

plt.errorbar(x, ProtA[0], ProtA[1], fmt=' ', capsize=2,
             label="Protons " + A + " on " + B + " Min=" + str(round(np.min(ProtA[0]), 2)) + " kRad at " + str(
                 round(np.argmin(ProtA[0]) + 1)) + " % " + A)

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

plt.figure(2)
plt.errorbar(x, TotalA[0], TotalA[1], fmt=' ', capsize=2,
             label=MatA + " on top of " + MatB + " Min=" + str(round(np.min(TotalA[0]), 2)) + " krad at " + str(
             round(np.argmin(TotalA[0]) + 1)) + " % " + A)

plt.ylim(0, Ymax)
plt.title(
    "Total dose deposited by trapped particles in 0.5 mm Si \n behind " + ShieldingDepth + "/cm2 of " + MatA + "-" + MatB + " shielding")  # --------
plt.xlabel("Percentage of shielding mass in " + MatA + " [%]")
plt.ylabel("Deposited ionising dose [krad]")
plt.grid(which='both')
plt.legend()
# plt.show()
plt.savefig(Path + Shield + "-GradientSum.eps", format='eps', bbox_inches="tight")