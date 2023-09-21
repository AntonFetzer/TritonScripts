import numpy as np
import matplotlib.pyplot as plt
from instrument_simulation.Dependencies.MeVtokRad_2D import MeVtokRad_2D
import nmmn.plots

Path = "/l/triton_work/Gradient/3Material2-5gcm2/Pe-Pb-W/csv/"

MatA = "Polyethylene"
MatB = "Lead"
MatC = "Tungsten"
MatD = "Tungsten"
a = "pe-"
b = "pb-"
c = "w-"
d = "w-"

A = "PE"
B = "Pb"
C = "W"
D = "W"
Shield = A + "-" + B + "-" + C + "-" + D

Nelec = "2e9"
Nprot = "1e8"

ElecFile = a + b + c + d + Nelec + "electron.txt"
ProtFile = a + b + c + d + Nprot + "proton.txt"

NORM_FACTOR_SPECTRUM_Elec = 5.886798E+14  # Elec500keV
NORM_FACTOR_SPECTRUM_Prot = 3.381390E+11  # Prots10MeV
Plates = 32
Npart_Elec = float(Nelec) / (Plates*Plates)
Npart_Prot = float(Nprot) / (Plates*Plates)

ElecData = np.genfromtxt(Path + ElecFile, delimiter=',', dtype=None, encoding='ASCII')
ProtData = np.genfromtxt(Path + ProtFile, delimiter=',', dtype=None, encoding='ASCII')

EleckRad = np.zeros((Plates, Plates))
ProtkRad = np.zeros((Plates, Plates))
TotalkRad = np.zeros((Plates, Plates))

ElecErr = np.zeros((Plates, Plates))
ProtErr = np.zeros((Plates, Plates))
TotalErr = np.zeros((Plates, Plates))

i = 0
for y in range(Plates):
    for x in range(Plates):
        EleckRad[x][y] = MeVtokRad_2D(ElecData[i][1], NORM_FACTOR_SPECTRUM_Elec, Npart_Elec)
        ProtkRad[x][y] = MeVtokRad_2D(ProtData[i][1], NORM_FACTOR_SPECTRUM_Prot, Npart_Prot)

        ElecErr[x][y] = MeVtokRad_2D(ElecData[i][2], NORM_FACTOR_SPECTRUM_Elec, Npart_Elec)
        ProtErr[x][y] = MeVtokRad_2D(ProtData[i][2], NORM_FACTOR_SPECTRUM_Prot, Npart_Prot)

        TotalkRad[x][y] = EleckRad[x][y] + ProtkRad[x][y]
        TotalErr[x][y] = np.sqrt(ElecErr[x][y] * ElecErr[x][y] + ProtErr[x][y] * ProtErr[x][y])

        i += 1

# Read ------------------------------------------------------------------

turbo = nmmn.plots.turbocmap()

fig1 = plt.figure(1)
plt.imshow(np.rot90(EleckRad), cmap=turbo, extent=(-0.5, Plates-0.5, -0.5, Plates-0.5))
x, y = np.unravel_index(np.argmin(EleckRad), EleckRad.shape)
plt.plot(x, y, 'wx')
plt.title("Electron dose behind 2.5g/cm2 of multilayer shielding\n Layer C: " + MatC + "    Min=" + str(round(np.min(EleckRad), 2)) + " krad    Layer D: " + MatD)
plt.xlabel("     x-index -->")
plt.ylabel("y-index -->")
plt.figtext(0.08, 0.025, "Layer A: " + MatA, size="large")
plt.figtext(0.6, 0.025, "Layer B: " + MatB, size="large")
cbar = plt.colorbar()
cbar.set_label("Dose in krad")
plt.savefig(Path + Shield + "-Electron-Gradient.eps", format='eps', bbox_inches="tight")

fig2 = plt.figure(2)
plt.imshow(np.rot90(ProtkRad), cmap=turbo, extent=(-0.5, Plates-0.5, -0.5, Plates-0.5))
x, y = np.unravel_index(np.argmin(ProtkRad), ProtkRad.shape)
plt.plot(x, y, 'wx')
plt.title("Proton dose behind 2.5g/cm2 of multilayer shielding\n Layer C: " + MatC + "    Min=" + str(round(np.min(ProtkRad), 2)) + " krad    Layer D: " + MatD)
plt.xlabel("     x-index -->")
plt.ylabel("y-index -->")
plt.figtext(0.08, 0.025, "Layer A: " + MatA, size="large")
plt.figtext(0.6, 0.025, "Layer B: " + MatB, size="large")
cbar = plt.colorbar()
cbar.set_label("Dose in krad")
plt.savefig(Path + Shield + "-Proton-Gradient.eps", format='eps', bbox_inches="tight")


fig3 = plt.figure(3)
plt.imshow(np.rot90(TotalkRad), cmap=turbo, extent=(-0.5, Plates-0.5, -0.5, Plates-0.5))
x, y = np.unravel_index(np.argmin(TotalkRad), TotalkRad.shape)
plt.plot(x, y, 'wx')
plt.title("Total dose behind 2.5g/cm2 of multilayer shielding\n Layer C: " + MatC + "    Min=" + str(round(np.min(TotalkRad), 2)) + " krad    Layer D: " + MatD)
plt.xlabel("     x-index -->")
plt.ylabel("y-index -->")
plt.figtext(0.08, 0.025, "Layer A: " + MatA, size="large")
plt.figtext(0.6, 0.025, "Layer B: " + MatB, size="large")
cbar = plt.colorbar()
cbar.set_label("Dose in krad")
plt.savefig(Path + Shield + "-Total-Gradient.eps", format='eps', bbox_inches="tight")


# Numerical Results -------------------------------------------------

TotalMin = np.min(TotalkRad)
# x, y = np.unravel_index(np.argmin(TotalkRad), TotalkRad.shape)

print(x, y)

TotalMinErr = TotalErr[x][y]

Afrac = (31-x) * (31-y) / 9.61
Bfrac = x * (31-y) / 9.61
Cfrac = (31-x) * y / 9.61
Dfrac = x * y / 9.61

CSVFile = open(Path + a + b + c + d + "Results.txt", 'w')
CSVFile.writelines("A%, Material A, B%, Material B, C%, Material C, D%, Material D, Total Dose, Total Error \n")

List = (Afrac, A, Bfrac, B, Cfrac, C, Dfrac, D, TotalMin, TotalMinErr)

String = ', '.join(map(str, List))
print(String)
CSVFile.writelines(String + "\n")
CSVFile.close()

