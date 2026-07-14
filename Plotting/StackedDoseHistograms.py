import matplotlib.pyplot as plt
from Dependencies.TotalDoseHistograms import totalDoseHistograms
import os
from natsort import natsorted
import numpy as np

Path = "/l/triton_work/LunarBackscatter/LunarGCR-WideHist/"

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

for i in range(len(Folders)):
    DoseHist, PrimaryHist = totalDoseHistograms(Path + Folders[i] + "/Res/")
    TotalNumberEntries.append(sum(DoseHist['entries']))
    print("Total number of particles = " + f"{TotalNumberEntries[i]:.3}")
    print("Number of Data points for", str(ThickList[i]) + "mm:", "{:,}".format(len(DoseHist['lower'])))

    DoseHists.append(DoseHist)
    PrimaryHists.append(PrimaryHist)

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

for i in range(len(Folders)):
    plt.bar(PrimaryHists[i]['lower'], PrimaryHists[i]['value'], width=PrimaryHists[i]['upper'] - PrimaryHists[i]['lower'], align='edge', color=Colours[i], label=str(ThickList[i]) + "mm Al")

plt.yscale("log")
plt.xscale("log")
#plt.xlim(8, 1200)
plt.grid()
plt.title("Dose deposited VS primary kinetic energy\n" + f"{TotalNumberEntries[0]:.3}" + " particles each")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Dose [krad per Month]")
plt.legend()
plt.savefig(Path + "DoseVSPrimary.pdf", format='pdf', bbox_inches="tight")


plt.figure(3)

for i in range(len(Folders)):
    plt.bar(PrimaryHists[i]['lower'], PrimaryHists[i]['entries'], width=PrimaryHists[i]['upper'] - PrimaryHists[i]['lower'], align='edge', color=Colours[i], label=str(ThickList[i]) + "mm Al")

plt.yscale("log")
plt.xscale("log")
#plt.xlim(8, 1200)
plt.grid()
plt.title("Particle count VS primary kinetic energy\n" + f"{TotalNumberEntries[0]:.3}" + " particles each")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Number of entries")
plt.legend()
plt.savefig(Path + "PrimaryHistogram.pdf", format='pdf', bbox_inches="tight")
