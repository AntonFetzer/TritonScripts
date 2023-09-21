from GRAS.Dependencies.TotalKRadGras import totalkRadGras
import numpy as np
import matplotlib.pyplot as plt
import os
from natsort import natsorted
from GRAS.Dependencies.TotalKRadGras import totalkRadGras

Path = "/l/triton_work/ShieldingCurves/Carrington/32mm/Res/"

Data = totalkRadGras(Path, "Elec")

NumTiles = np.shape(Data)[1]

fig1 = plt.figure(1)

x = np.linspace(0, 32, num=NumTiles, endpoint=True)

# Dose is in per month.
# To get the total dose after the carrington event with 29.52 hours of flux
ScalingFactor = 29.52/(30*24)

plt.errorbar(x, Data[0] * ScalingFactor, Data[1] * ScalingFactor, fmt=' ', capsize=5, label="Planar Shielding")

####### Plot 10kRad line #########
CriticalDose = [10 for i in x]
plt.plot(x, CriticalDose, color='k', linewidth=2, label='10 krad')
CriticalDose = [100 for i in x]
plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='100 krad')

######################  3D data  ##################################

Path = "/l/triton_work/1U/"

Folders = [f for f in os.listdir(Path) if f.endswith('mm')]
Folders = natsorted(Folders)

print(Folders)

ThickList = []

for folder in Folders:
    ThickList.append(int(folder.split("mm")[0]))

print(ThickList)

Colours = ['C3', 'C1', 'C8', 'C2', 'C9', 'C7']

for i, folder in enumerate(Folders):
    Electrons = totalkRadGras(Path + folder + "/Res/", "Elec")

    x = []
    y = []
    Err = []
    for vol in range(len(Electrons[0])):
        x.append(ThickList[i])
        y.append(Electrons[0][vol] * ScalingFactor)
        Err.append(Electrons[1][vol] * ScalingFactor)
    plt.errorbar(x, y, Err, label="3D " + folder, fmt='+', capsize=8, markersize=10, color=Colours[i])
    print(np.average(x), np.average(y), max(y) - min(y))

plt.ylim(0.5, 1e4)
#plt.ylim(0, 120)
plt.xlim(-0.5, 17.5)
plt.grid(which="major")
plt.yscale("log")
plt.title("Carrington Total Electron Fluence\nTotal Dose deposited behind planar aluminium shielding")
plt.xlabel("Aluminium Shielding Thickness [mm]")
plt.ylabel("Ionising Dose [krad]")
plt.legend()
plt.savefig("/l/TritonPlots/Carrington/CarringtonShieldingCurveWith3D.eps", format='eps', bbox_inches="tight")

#plt.show()
