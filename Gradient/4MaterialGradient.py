import numpy as np
import matplotlib.pyplot as plt
from Dependencies.MeVtokRad_2D import MeVtokRad_2D
import nmmn.plots

turbo = nmmn.plots.turbocmap()

Path = "/home/anton/Desktop/triton_work/Gradient/3Material1-5gcm2/Al-PE-W/csv/"

Plates = 32

Shield = "Al-PE-W-W"
MatA = "Aluminium"
MatB = "Polyethylene"
MatC = "Tungsten"
MatD = "Tungsten"

ElecFile = "al-pe-w-w-2e8electron.txt"
ProtFile = "al-pe-w-w-2e6proton.txt"

NORM_FACTOR_SPECTRUM_Elec = 5.886798E+14  # Elec500keV
NORM_FACTOR_SPECTRUM_Prot = 3.381390E+11  # Prots10MeV
Npart_Elec = 2e8 / (Plates * Plates)
Npart_Prot = 2e6 / (Plates * Plates)

ElecMeV = np.loadtxt(Path + ElecFile)
ProtMeV = np.loadtxt(Path + ProtFile)

EleckRad = np.zeros((Plates, Plates))
ProtkRad = np.zeros((Plates, Plates))

i = 0
for y in range(Plates):
    for x in range(Plates):
        EleckRad[x][y] = MeVtokRad_2D(ElecMeV[i], NORM_FACTOR_SPECTRUM_Elec, Npart_Elec)
        ProtkRad[x][y] = MeVtokRad_2D(ProtMeV[i], NORM_FACTOR_SPECTRUM_Prot, Npart_Prot)
        i += 1

fig1 = plt.figure(1)  # , figsize=(5, 5)
plt.imshow(EleckRad, cmap=turbo, extent=(0, Plates, Plates, 0))
plt.title("Electrons vs 1.5g/cm2 " + Shield + " Min=" + str(round(np.min(EleckRad), 2)) + " krad \n" + MatA + " --> " + MatC)
plt.xlabel(MatB + " --> " + MatD)
plt.ylabel(MatB + " <-- " + MatA)
cbar = plt.colorbar()
cbar.set_label("Dose in krad")
plt.savefig(Path + Shield + "-Electron-Gradient.eps", format='eps')

fig2 = plt.figure(2)
plt.imshow(ProtkRad, cmap=turbo, extent=(0, Plates, Plates, 0))
plt.title("Proton vs 1.5g/cm2 " + Shield + " Min=" + str(round(np.min(ProtkRad), 2)) + " krad \n" + MatA + " --> " + MatC)
plt.xlabel(MatB + " --> " + MatD)
plt.ylabel(MatB + " <-- " + MatA)
cbar = plt.colorbar()
cbar.set_label("Dose in krad")
plt.savefig(Path + Shield + "-Proton-Gradient.eps", format='eps')

TotalkRad = EleckRad + ProtkRad
fig3 = plt.figure(3)
plt.imshow(TotalkRad, cmap=turbo, extent=(0, Plates, Plates, 0))
plt.title("Trapped Particles vs 1.5g/cm2 " + Shield + " Min=" + str(round(np.min(TotalkRad), 2)) + " krad \n" + MatA + " --> " + MatC)
plt.xlabel(MatB + " --> " + MatD)
plt.ylabel(MatB + " <-- " + MatA)
cbar = plt.colorbar()
cbar.set_label("Dose in krad")
plt.savefig(Path + Shield + "-Total-Gradient.eps", format='eps')

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
