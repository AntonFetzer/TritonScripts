import matplotlib.pyplot as plt
from GRAS.Dependencies.TotalGRASHistos import totalGRASHistos
import os
from natsort import natsorted
import numpy as np

Num = 7

Path = "/home/anton/Desktop/triton_work/Histograms/StackedCarringtonSEP/"

lowerID = 0
upperID = 1
meanID = 2
valueID = 3
errorID = 4
entriesID = 5

Folders = [f for f in os.listdir(Path) if f.endswith('mm')]
Folders = natsorted(Folders)

print(Folders)

ThickList = []

for folder in Folders:
    ThickList.append(int(folder.split("m")[0]))

print(ThickList)

DoseHists = []
PrimaryHists = []
TotalNumberEntries = []

for i in range(Num):
    Data = totalGRASHistos(Path + Folders[i] + "/Res/", "")
    #print("Shape of Data:", np.shape(Data))
    TotalNumberEntries.append(sum(Data[0][:, entriesID]))
    print("Total number of particles = " + f"{TotalNumberEntries[i]:.3}")
    print("Number of Data points for", str(ThickList[i]) + "mm:", "{:,}".format(len(Data[0])))

    DoseHists.append(Data[0])
    PrimaryHists.append(Data[1])

print("Shape of DoseHists:", np.shape(DoseHists))
# [Thickness][Bin][Variable]

Colours = ['C3', 'C1', 'C8', 'C2', 'C9', 'C0', 'C7']
#          red, orange,yellow,green,cyan,blue, grey
#Colours = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
#           blue,orange,green,red,violet,brown

# C0    blue
# C1    Orange
# C2    Green
# C3    red
# C4    purple
# C5    brown
# C6    pink
# C7    grey
# C8    yellow
# C9    cyan


plt.figure(2)

for i in range(Num):
    plt.bar(PrimaryHists[i][:, lowerID], PrimaryHists[i][:, valueID], width=PrimaryHists[i][:, upperID] - PrimaryHists[i][:, lowerID], align='edge', color=Colours[i], label=str(ThickList[i]) + "mm Al")

plt.yscale("log")
plt.xscale("log")
plt.xlim(5, 55)
plt.grid()
plt.title("Dose deposited VS primary kinetic energy\n" + f"{TotalNumberEntries[0]:.3}" + " particles each")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Dose [krad per Month]")
plt.legend()
plt.savefig(Path + "DoseVSPrimary.pdf", format='pdf', bbox_inches="tight")


plt.figure(3)

for i in range(Num):
    plt.bar(PrimaryHists[i][:, lowerID], PrimaryHists[i][:, entriesID], width=PrimaryHists[i][:, upperID] - PrimaryHists[i][:, lowerID], align='edge', color=Colours[i], label=str(ThickList[i]) + "mm Al")

plt.yscale("log")
plt.xscale("log")
plt.xlim(5, 55)
plt.grid()
plt.title("Particle count VS primary kinetic energy\n" + f"{TotalNumberEntries[0]:.3}" + " particles each")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Number of entries")
plt.legend()
plt.savefig(Path + "PrimaryHistogram.pdf", format='pdf', bbox_inches="tight")
