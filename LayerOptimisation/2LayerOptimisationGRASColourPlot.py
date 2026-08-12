import numpy as np
import matplotlib.pyplot as plt
from Dependencies.MergeTotalDose import mergeTotalDose
from Dependencies.TotalDose import totalDose
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

    ElecA = totalDose(path, filename_contains="Electrons2MatA")
    ProtA = totalDose(path, filename_contains="Protons2MatA")

    ElecB = totalDose(path, filename_contains="Electrons2MatB")
    ProtB = totalDose(path, filename_contains="Protons2MatB")

    ElecB = {key: np.flip(values) for key, values in ElecB.items()}
    ProtB = {key: np.flip(values) for key, values in ProtB.items()}

    TotalA = mergeTotalDose([ElecA, ProtA])["dose"]
    TotalB = mergeTotalDose([ElecB, ProtB])["dose"]

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


