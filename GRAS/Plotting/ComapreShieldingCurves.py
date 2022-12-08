from GRAS.Dependencies.TotalKRadGras import totalkRadGras
import numpy as np
import matplotlib.pyplot as plt
import os
from natsort import natsorted

Paths = ["/home/anton/Desktop/triton_work/ShieldingCurves/SuperGTO/Prot/",
         "/home/anton/Desktop/triton_work/ShieldingCurves/SuperGTO/SolProt/",
         "/home/anton/Desktop/triton_work/ShieldingCurves/SuperGTO/SolHe/",
         "/home/anton/Desktop/triton_work/ShieldingCurves/SuperGTO/SolO/",
         "/home/anton/Desktop/triton_work/ShieldingCurves/SuperGTO/SolFe/"]

Labels = ["Trapped Protons", "Solar Protons", "Solar Helium", "Solar Oxygen", "Solar Iron"]

Data = []
for path in Paths:
    Data.append(totalkRadGras(path, ""))

NumTiles = np.shape(Data[0])[1]

plt.figure(1)

x = np.linspace(0, 9.626, num=NumTiles, endpoint=True)

for i, data in enumerate(Data):
    plt.errorbar(x, data[0], data[1], fmt=' ', markersize=4, capsize=5, label=Labels[i])

plt.xlim(-0.5, 9.626)
plt.grid(which="both")
plt.yscale("log")
plt.title("Shielding curve comparison")
plt.xlabel("Aluminium Shielding Thickness [mm]")
plt.ylabel("Ionising Dose [krad]")
plt.legend()

####### Plot 10kRad line #########
CriticalDose = [10 for i in x]
plt.plot(x, CriticalDose, color='k', linewidth=2, label='10 krad')
CriticalDose = [100 for i in x]
plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='100 krad')

plt.savefig("/home/anton/Desktop/triton_work/ShieldingCurves/SuperGTO/ShieldingCurveComparisonOld.eps", format='eps', bbox_inches="tight")


plt.figure(2)
DataDiff = abs(Data[0] - Data[1])

plt.errorbar(x, DataDiff[0], DataDiff[1], fmt=' ', capsize=5, label="DataDiff")
plt.yscale("log")
plt.title("DataDiff")
plt.xlabel("Entry#")
plt.ylabel("Difference in Ionising Dose [krad]")
plt.legend()

plt.figure(3)
DataRatio = Data[0] / Data[1]

plt.scatter(x, DataRatio[0], label="DataRatio")
plt.yscale("log")
plt.title("DataRatio")
plt.xlabel("Entry#")
plt.ylabel("Difference in Ionising Dose [krad]")
plt.legend()

plt.show()

print(np.average(DataRatio[0]))
