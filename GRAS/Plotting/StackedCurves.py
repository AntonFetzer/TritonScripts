import numpy as np
import matplotlib.pyplot as plt
from GRAS.Dependencies.TotalDose import totalkRadGras
from scipy import signal

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
    Elec = totalkRadGras(Path + Shield + D + "/Res/", "Electrons")
    Prot = totalkRadGras(Path + Shield + D + "/Res/", "Protons")

    Total = Elec
    Total[0] = Elec[0] + Prot[0]
    Total[1] = np.sqrt(Elec[1] * Elec[1] + Prot[1] * Prot[1])

    plt.errorbar(x, Total[0], Total[1], fmt=' ', capsize=2, label= str(Depths[i]) + " g/cm2 Shielding depth")
    #plt.errorbar(x, Elec[0], Elec[1], fmt=' ', capsize=2, label="Shielding depth =" + str(D))
    #plt.errorbar(x, Prot[0], Prot[1], fmt=' ', capsize=2, label="Shielding depth =" + str(D))

    #print("The relative error for D= " + str(Depths[i]) + " is " + str( 100 * sum(Total[1]) / sum(Total[0]) ) + " %")
    TotalFiltered = signal.savgol_filter(Total, 1, 0)

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
