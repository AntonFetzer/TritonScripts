import numpy as np
import matplotlib.pyplot as plt
from Dependencies.MergeTotalDose import mergeTotalDose
from Dependencies.TotalDose import totalDose

Path = "/l/triton_work/2LayerStackedCurves/"
Shield = "PE-Pb"
DepthsStr = ["-02", "-04", "-08", "-16", "-32"]
#DepthsStr = ["-02", "-04"]
Depths = [0.2, 0.4, 0.8, 1.6, 3.2]

A = "PE"
B = "Pb"
MatA = "Polyethylene"
MatB = "Lead"

x = np.linspace(0, 100, num=101, dtype=int)

plt.figure(1, [6.5, 8])

for i, D in enumerate(DepthsStr):
    Elec = totalDose(Path + Shield + D + "/Res/", filename_contains="Electrons")
    Prot = totalDose(Path + Shield + D + "/Res/", filename_contains="Protons")

    Total = mergeTotalDose([Elec, Prot])

    plt.errorbar(x, Total["dose"], Total["error"], fmt=' ', capsize=2, label= str(Depths[i]) + " g/cm2 Shielding depth")
    #plt.errorbar(x, Elec["dose"], Elec["error"], fmt=' ', capsize=2, label="Shielding depth =" + str(D))
    #plt.errorbar(x, Prot["dose"], Prot["error"], fmt=' ', capsize=2, label="Shielding depth =" + str(D))

    #print("The relative error for D= " + str(Depths[i]) + " is " + str(100 * sum(Total["error"]) / sum(Total["dose"])) + " %")

#plt.ylim(5e-2, 2e2)
#plt.ylim(1e-1, 2e2)
plt.title("Dose deposited by trapped particles in 0.5 mm Si \n behind " + MatA + " on top of " + MatB)
plt.xlabel("Percentage of shielding mass in " + MatA + " [%]")
plt.ylabel("Deposited ionising dose [kRad]")
plt.grid(which='both')
plt.legend(ncol=2, loc='center', bbox_to_anchor=(0.5, 0.6))
plt.yscale("log")
#plt.show()
plt.savefig(Path + "StackedCurves.eps", format='eps', bbox_inches="tight")
