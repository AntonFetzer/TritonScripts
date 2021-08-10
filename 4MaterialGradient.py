import numpy as np
import matplotlib.pyplot as plt
from MeVtokRad_2D import MeVtokRad_2D

Path = "/home/anton/Desktop/triton_work/Gradient/4MaterialGradient/"

Shield = "Al-Pe-Al-Pe"
MatA = "Aluminium"
MatB = "Polyethylene"
MatC = "Aluminium"
MatD = "Polyethylene"

ElecFile = "gradient-al-pe-al-pe-1e8electron.txt"
ProtFile = "gradient-al-pe-al-pe-1e8proton.txt"

NORM_FACTOR_SPECTRUM_Elec = 5.886798E+14  # Elec500keV
NORM_FACTOR_SPECTRUM_Prot = 3.381390E+11  # Prots10MeV
Npart_Elec = 1e8 / 100
Npart_Prot = 1e8 / 100

ElecMeV = np.loadtxt(Path + ElecFile)
ProtMeV = np.loadtxt(Path + ProtFile)

EleckRad = np.zeros((10, 10))
ProtkRad = np.zeros((10, 10))

i = 0
for y in range(10):
    for x in range(10):
        EleckRad[x][y] = MeVtokRad_2D(ElecMeV[i], NORM_FACTOR_SPECTRUM_Elec, Npart_Elec)
        ProtkRad[x][y] = MeVtokRad_2D(ProtMeV[i], NORM_FACTOR_SPECTRUM_Prot, Npart_Prot)
        i += 1

fig1 = plt.figure(1, figsize=(10, 10))
plt.imshow(EleckRad, extent=(5, 105, 5, 105))
plt.jet()
plt.title("Electron dose map behind gradient shielding" + Shield)
plt.xlabel("Gradient layer 1 " + MatA + " to layer 2 " + MatB)
plt.ylabel("Gradient layer 1 " + MatA + " to layer 3 " + MatC)
plt.colorbar()
plt.show()

fig2 = plt.figure(2, figsize=(10, 10))
plt.imshow(ProtkRad, extent=(5, 105, 5, 105))
plt.jet()
plt.title("Proton dose map behind gradient shielding" + Shield)
plt.xlabel("Gradient layer 1 " + MatA + " to layer 2 " + MatB)
plt.ylabel("Gradient layer 1 " + MatA + " to layer 3 " + MatC)
plt.colorbar()
plt.show()

TotalkRad = EleckRad + ProtkRad
fig3 = plt.figure(3, figsize=(10, 10))
plt.imshow(TotalkRad, extent=(5, 105, 5, 105))
plt.jet()
plt.title("Total dose map behind gradient shielding" + Shield)
plt.xlabel("Gradient layer 1 " + MatA + " to layer 2 " + MatB)
plt.ylabel("Gradient layer 1 " + MatA + " to layer 3 " + MatC)
plt.colorbar()
plt.show()









'''
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
'''
