import matplotlib.pyplot as plt
from ReadG4root import readG4root
import numpy as np
import os
from natsort import natsorted
import sys

# ---------------------- Dont Modify -------- Dont break compatability -----------------------

# ---------------------- Dont Modify -------- Dont break compatability -----------------------

# ---------------------- Dont Modify -------- Dont break compatability -----------------------

# Path = "/home/anton/Desktop/triton_work/MULASS/"
Path = "/scratch/work/fetzera1/MULASS/"
# Folder = "MulasTrapProton2e9AlFull"
Folder = sys.argv[1]

RootFiles = [f for f in os.listdir(Path + Folder + "/root/") if f.endswith('.root')]
RootFiles = natsorted(RootFiles)

ThickList = [int(f.split(".")[0]) for f in RootFiles]

i = 0
Data = []
for Thick in ThickList:
    Data.append(readG4root(Path + Folder + "/root/" + RootFiles[i]))
    print("Number of Data points for", str(int(Thick)) + "mm:", "{:,}".format(len(Data[i][0])))
    i += 1

Colours = ['C3', 'C1', 'C2', 'C9', 'C0', 'C7']
# Colours = ['C4', 'C3', 'C2', 'C1', 'C0']
# ----------------- Primary Histogram ----------------------------
plt.figure(1)
i = 0
for Thick in ThickList:
    plt.hist(Data[i][0], bins=100, label=str(int(Thick)) + "mm Al", color=Colours[i])
    i += 1

plt.yscale("log")
plt.grid(which='major')
plt.title("Histogram of primary energies \n" + Folder)
plt.xlabel("Primary particle energy [MeV]")
plt.ylabel("Number of counts per primary energy bin")
plt.legend()
# plt.show()
plt.savefig(Path + Folder + "/Img/" + Folder + "PrimaryHist.eps", format='eps')

# ----------------- Dose Histogram ----------------------------
plt.figure(2)
i = 0
for Thick in ThickList:
    plt.hist(Data[i][1], bins=100, label=str(int(Thick)) + "mm Al", color=Colours[i])
    i += 1

plt.yscale("log")
plt.grid(which='major')
plt.title("Histogram of dose depositions\n" + Folder)
plt.xlabel("Deposited dose [MeV]")
plt.ylabel("Number of counts per dose bin")
plt.legend()
# plt.show()
plt.savefig(Path + Folder + "/Img/" + Folder + "DoseHist.eps", format='eps')

# ----------------- Primaries weighted for Dose ----------------------------
plt.figure(3)
i = 0
for Thick in ThickList:
    plt.hist(Data[i][0], weights=Data[i][1], bins=100, label=str(int(Thick)) + "mm Al", color=Colours[i])
    i += 1

plt.yscale("log")
plt.grid(which='major')
plt.title("Histogram of dose depositions per primary energy bin\n" + Folder)
plt.xlabel("Primary particle energy [MeV]")
plt.ylabel("Dose deposited per primary energy bin [MeV]")
plt.legend()
# plt.show()
plt.savefig(Path + Folder + "/Img/" + Folder + "DoseWeight.eps", format='eps')

# ----------------- Primaries weighted for Dose ZOOM ----------------------------
plt.figure(4)
i = 0
for Thick in ThickList:
    plt.hist(Data[i][0], weights=Data[i][1], range=(0, 2), bins=100, label=str(int(Thick)) + "mm Al", color=Colours[i])
    i += 1

plt.yscale("log")
plt.grid(which='major')
plt.title("Histogram of dose depositions per primary energy bin\n" + Folder)
plt.xlabel("Primary particle energy [MeV]")
plt.ylabel("Dose deposited per primary energy bin [MeV]")
plt.legend()
# plt.show()
plt.savefig(Path + Folder + "/Img/" + Folder + "DoseWeightZOOM.eps", format='eps')

# ----------------- Secondaries Histogram ----------------------------
plt.figure(5)
i = 0
for Thick in ThickList:
    plt.hist(Data[i][0], weights=Data[i][2], bins=100, label=str(int(Thick)) + "mm Al", color=Colours[i])
    i += 1

plt.yscale("log")
plt.grid(which='major')
plt.title("Histogram of dose depositions pue to secondary particles\n" + Folder)
plt.xlabel("Primary particle energy [MeV]")
plt.ylabel("Dose from secondaries deposited per primary energy bin [MeV]")
plt.legend()
# plt.show()
plt.savefig(Path + Folder + "/Img/" + Folder + "SecondaryWeight.eps", format='eps')

# # srun --mem=50G --time=00:15:00 python ShieldedHistograms.py
