from GRAS.Dependencies.TotalKRadGras import totalkRadGras
import numpy as np
import matplotlib.pyplot as plt
import os
from natsort import natsorted

Paths = ["/home/anton/Desktop/triton_work/ShieldingCurves/Carrington-SEP-Plus2Sigma-Int-With0/Res/",
         "/home/anton/Desktop/triton_work/ShieldingCurves/Carrington-SEP-Expected-Int-With0/Res/",
         "/home/anton/Desktop/triton_work/ShieldingCurves/Carrington-SEP-Minus2Sigma-Int-With0/Res/",
         "/home/anton/Desktop/triton_work/ShieldingCurves/AP9-10MeV-32mm/Res/",
         "/home/anton/Desktop/triton_work/ShieldingCurves/CarringtonElectron-32mm/Res/",
         "/home/anton/Desktop/triton_work/ShieldingCurves/SEP2003-INTEGRAL-FluxBasedOnFluenceDividedBy24h/Res/"]

Labels = ["Carrington SEP +2 Sigma",
          "Carrington SEP Expected",
          "Carrington SEP -2 Sigma",
          "AP9 GT0",
          "Carrington Peak Electron Flux",
          "SEP 2003"]

Colours = ['C1', 'C0', 'C2', 'C3', 'C9', 'C10']

Data = []
for path in Paths:
    Data.append(totalkRadGras(path, "")/ 30 )  # DAILY dose

NumTiles = np.shape(Data[0])[1]

plt.figure(1)

x = np.linspace(0, 32, num=NumTiles, endpoint=True)

for i, data in enumerate(Data):
    plt.errorbar(x, data[0], data[1], fmt='o', markersize=4, capsize=5, label=Labels[i], color=Colours[i])

plt.xlim(-0.5, 20)
plt.grid()
plt.yscale("log")
plt.title("Daily Ionising Dose behind shielding of varying thickness")
plt.xlabel("Aluminium Shielding Thickness [mm]")
plt.ylabel("Ionising Dose per Day [krad]")
plt.legend()

####### Plot 10kRad line #########
CriticalDose = [10 for i in x]
plt.plot(x, CriticalDose, color='k', linewidth=2, label='10 krad')
CriticalDose = [100 for i in x]
plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='100 krad')

plt.savefig("/home/anton/Desktop/triton_work/ShieldingCurves/Plots/CarringtonCurve.pdf", format='pdf', bbox_inches="tight")

'''
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
'''