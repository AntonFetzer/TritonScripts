import numpy as np
import matplotlib.pyplot as plt
from MeVtokRad_2D import MeVtokRad_2D

Path = "/home/anton/Desktop/triton_work/Gradient/2MaterialGradient/pe-zn/csv/"

Shield = "Pe-Zn"
MatA = "Polyethylene"
MatB = "Zinc"

ElecFileA = "gradientpezn2e9electron.txt"
ElecFileB = "gradientznpe2e9electron.txt"
ProtFileA = "gradientpezn1e8proton.txt"
ProtFileB = "gradientznpe1e8proton.txt"

NORM_FACTOR_SPECTRUM_Elec = 5.886798E+14  # Elec500keV
NORM_FACTOR_SPECTRUM_Prot = 3.381390E+11  # Prots10MeV
Npart_Elec = 2e9 / 100
Npart_Prot = 1e8 / 100

ElecA = np.loadtxt(Path + ElecFileA)[:99]
ElecB = np.loadtxt(Path + ElecFileB)[:99]
ProtA = np.loadtxt(Path + ProtFileA)[:99]
ProtB = np.loadtxt(Path + ProtFileB)[:99]
TotalEdepA = np.zeros(99)
TotalEdepB = np.zeros(99)

for i in range(99):
    ElecA[i] = MeVtokRad_2D(ElecA[i], NORM_FACTOR_SPECTRUM_Elec, Npart_Elec)
    ElecB[i] = MeVtokRad_2D(ElecB[i], NORM_FACTOR_SPECTRUM_Elec, Npart_Elec)
    ProtA[i] = MeVtokRad_2D(ProtA[i], NORM_FACTOR_SPECTRUM_Prot, Npart_Prot)
    ProtB[i] = MeVtokRad_2D(ProtB[i], NORM_FACTOR_SPECTRUM_Prot, Npart_Prot)
    TotalEdepA[i] = ElecA[i] + ProtA[i]
    TotalEdepB[i] = ElecB[i] + ProtB[i]

x = np.linspace(0, 98, num=99, dtype=int)
print(x)

plt.figure(1)
plt.plot(x + 1, ElecA, '.', label="Electrons " + MatA + " on top of " + MatB)
plt.plot(x + 1, ElecB, '.', label="Electrons " + MatB + " on top of " + MatA)

plt.plot(x + 1, ProtA, '.', label="Protons " + MatA + " on top of " + MatB)
plt.plot(x + 1, ProtB, '.', label="Protons " + MatB + " on top of " + MatA)

plt.title("Dose deposited by trapped particles in 0.5 mm Si \n behind 2.5g/cm2 of " + Shield + " shielding")  # ---------
plt.xlabel("Aluminium mass ratio [%]")
plt.ylabel("Deposited ionising dose [krad]")
plt.grid(which='both')
plt.legend()
# plt.show()
plt.savefig(Path + Shield + "Gradient.eps", format='eps')

plt.figure(2)
plt.plot(x + 1, TotalEdepA, '.', label=MatA + " on top of " + MatB)
plt.plot(x + 1, TotalEdepB, '.', label=MatB + " on top of " + MatA)

plt.title("Total dose deposited by trapped particles in 0.5 mm Si \n behind 2.5g/cm2 of " + Shield + " shielding")  # ---------
plt.xlabel("Aluminium mass ratio [%]")
plt.ylabel("Deposited ionising dose [krad]")
plt.grid(which='both')
plt.legend()
# plt.show()
plt.savefig(Path + Shield + "GradientSum.eps", format='eps')
