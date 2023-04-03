from GRAS.Dependencies.TotalKRadGras import totalkRadGras
import numpy as np
import matplotlib.pyplot as plt
import os
from natsort import natsorted

ImagePath = "/home/anton/Desktop/TritonPlots/Luna/"

Paths = ["/home/anton/Desktop/triton_work/ShieldingCurves/LunarSEP10MeVFlux/Res/",
         #"/home/anton/Desktop/triton_work/ShieldingCurves/LunarSEPMissionFluence/Res/",
         "/home/anton/Desktop/triton_work/ShieldingCurves/LunarCosmic-H-Flux/Res/",
         # "/home/anton/Desktop/triton_work/ShieldingCurves/Carrington-SEP-Plus2Sigma-Int-With0/Res/",
         # "/home/anton/Desktop/triton_work/ShieldingCurves/Carrington-SEP-Expected-Int-With0/Res/",
         # "/home/anton/Desktop/triton_work/ShieldingCurves/Carrington-SEP-Minus2Sigma-Int-With0/Res/",
         # "/home/anton/Desktop/triton_work/ShieldingCurves/SEP2003-INTEGRAL-FluxBasedOnFluenceDividedBy24h/Res/",
         "/home/anton/Desktop/triton_work/ShieldingCurves/AP9-10MeV-32mm/Res/",
         "/home/anton/Desktop/triton_work/ShieldingCurves/AE9-500keV-32mm/Res/",
         # "/home/anton/Desktop/triton_work/ShieldingCurves/ISS-LEO-Electron500keV/Res/",
         # "/home/anton/Desktop/triton_work/ShieldingCurves/CarringtonElectron-32mm/Res/"
         ]

Labels = ["LunarSEP",
          #"LunarSEPMissionFluence",
          "Cosmic H",
          # "Carrington SEP +2 Sigma",
          # "Carrington SEP EVT",
          # "Carrington SEP -2 Sigma",
          # "2003 SPE",
          "AP9 GTO trapped protons",
          "AE9 GTO trapped electrons",
          # "Carrington Peak Electron Flux"
          ]

Colours = ['C1', 'C0', 'C2', 'C8', 'C3', 'C9']

Data = []
for i, path in enumerate(Paths):
    Data.append(totalkRadGras(path, "") * 12 )  # Yearly dose

NumTiles = np.shape(Data[0])[1]

plt.figure(1)

x = np.linspace(0, 32, num=NumTiles, endpoint=True)

# plt.fill_between(x, Data[0][0], Data[1][0], color='C1', alpha=0.5)
# plt.fill_between(x, Data[1][0], Data[2][0], color='C2', alpha=0.5)

for i, data in enumerate(Data):
    plt.errorbar(x, data[0], data[1], fmt='.', markersize=5, capsize=5, label=Labels[i], color=Colours[i])

# LEO = totalkRadGras("/home/anton/Desktop/triton_work/ShieldingCurves/ISS-LEO-Electron500keV/Res/", "")
# plt.errorbar(x, LEO[0]*12*10, LEO[1], fmt='.', markersize=5, capsize=5, label="10 Year LEO Electrons")

# GTO = totalkRadGras("/home/anton/Desktop/triton_work/ShieldingCurves/AE9-500keV-32mm/Res/", "")
# plt.errorbar(x, GTO[0]*12, GTO[1], fmt='.', markersize=5, capsize=5, label="1 Year GTO Electrons")

####### Plot 10kRad line #########
CriticalDose = [10 for i in x]
plt.plot(x, CriticalDose, color='k', linewidth=2, label='10 krad')
CriticalDose = [100 for i in x]
plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='100 krad')

# plt.ylim(3e-3, 1e3)
plt.xlim(-0.5, 10)
plt.grid()
plt.yscale("log")
#plt.xscale("log")
plt.title("Yearly Ionising Dose behind shielding of varying thickness")
plt.xlabel("Aluminium Shielding Thickness [mm]")
plt.ylabel("Ionising Dose per Year [krad]")
plt.legend()

# plt.savefig("/home/anton/Desktop/triton_work/ShieldingCurves/Plots/CarringtonCurve.pdf", format='pdf', bbox_inches="tight")
plt.savefig(ImagePath + "/ShieldingCurves.svg", format='svg', bbox_inches="tight")

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
