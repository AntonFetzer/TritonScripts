import numpy as np
import matplotlib.pyplot as plt
from GRAS.Dependencies.TotalDose import totalkRadGras
from matplotlib import cm

Path = "/l/triton_work/2LayerOptimisation"

Folders = ["PE-Pb", "PE-W", "Ti-PE", "PE-FR4", "PE-Al", "Al-W", "FR4-Ti", "Al-FR4", "Kevlar-Al"]

Matrix = np.zeros((len(Folders), 101), dtype=np.float64)
Mins = []

LeftAxis = []
RightAxis = []

for i, folder in enumerate(Folders):
    path = Path + "/" + folder + "/Res/"
    print(path)

    ElecA = totalkRadGras(path, "Electrons2MatA")
    ProtA = totalkRadGras(path, "Protons2MatA")

    ElecB = totalkRadGras(path, "Electrons2MatB")
    ProtB = totalkRadGras(path, "Protons2MatB")

    ElecB = np.flip(ElecB, 1)
    ProtB = np.flip(ProtB, 1)

    TotalA = ElecA[0] + ProtA[0]
    TotalB = ElecB[0] + ProtB[0]

    if(min(TotalA) > min(TotalB)):
        Total = TotalB
        print(folder.split(sep="-")[1], "on top of", folder.split(sep="-")[0])
        LeftAxis.append(folder.split(sep="-")[1])
        RightAxis.append(folder.split(sep="-")[0])
    else:
        Total = np.flip(TotalA)
        print(folder.split(sep="-")[0], "on top of", folder.split(sep="-")[1])
        LeftAxis.append(folder.split(sep="-")[0])
        RightAxis.append(folder.split(sep="-")[1])

    print(folder, min(Total))
    Mins.append(min(Total))
    Matrix[i] = Total

#np.savetxt(Path + "/Matrix", Matrix)
#Matrix = np.loadtxt(Path + "/Matrix")

plt.imshow(Matrix, cmap=cm.viridis, aspect='auto')
cbar = plt.colorbar(orientation="horizontal")
cbar.set_label("Total ionizing dose per month in krad")

ax1 = plt.gca()
ax1.set_yticks(np.arange(9), LeftAxis)

ax2 = ax1.secondary_yaxis("right")
ax2.set_yticks(np.arange(9), RightAxis)

plt.title("Total Ionizing dose for two layer shielding")
plt.xlabel("Percentage of shielding mass in bottom layer [%]")

plt.savefig(Path + "/Summary.pdf", format='pdf', bbox_inches="tight")

#print(Mins)

plt.show()




