import matplotlib.pyplot as plt
from Dependencies.ReadG4root import readG4root
import os
from natsort import natsorted
import sys

# ---------------------- Dont Modify -------- Dont break compatability -----------------------

# ---------------------- Dont Modify -------- Dont break compatability -----------------------

# ---------------------- Dont Modify -------- Dont break compatability -----------------------

# Path = "/l/triton_work/MulasTest2/"
# Folder = "ProtonsFull"

Path = "/scratch/work/fetzera1/MulasTest2/"
Folder = sys.argv[1]
Max = int(sys.argv[2])

RootFiles = [f for f in os.listdir(Path + Folder + "/root/") if f.endswith('.root')]
RootFiles = natsorted(RootFiles)

ThickList = [0, 1, 2, 4, 8, 16]

i = 0
Data = []
for Thick in ThickList:
    Data.append(readG4root(Path + Folder + "/root/" + RootFiles[i]))
    print("Number of Data points for", str(int(Thick)) + "mm:", "{:,}".format(len(Data[i][0])))
    i += 1

Colours = ['C3', 'C1', 'C8', 'C2', 'C9', 'C0', 'C5']
#          red, orange, yellow,green, cyan, blue, grey
#Colours = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
#         blue,orange,green,red, violet, brown
#Colours = ['C6', 'C7', 'C8', 'C9', 'C10', 'C11']
#         pink,grey,yellow,cyan, blue, orange
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
plt.savefig(Path + Folder + "/img/" + Folder + "PrimaryHist.eps", format='eps', bbox_inches="tight")

# ----------------- Primary Histogram ZOOM ----------------------------
plt.figure(2)

i = 0
for Thick in ThickList:
    plt.hist(Data[i][0], bins=100, range=(0, Max), label=str(int(Thick)) + "mm Al", color=Colours[i])
    i += 1

plt.yscale("log")
plt.grid(which='major')
plt.title("Histogram of primary energies \n" + Folder)
plt.xlabel("Primary particle energy [MeV]")
plt.ylabel("Number of counts per primary energy bin")
plt.legend()
# plt.show()
plt.savefig(Path + Folder + "/img/" + Folder + "PrimaryHistZoom.eps", format='eps', bbox_inches="tight")

# ----------------- Dose Histogram ----------------------------
plt.figure(3)
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
plt.savefig(Path + Folder + "/img/" + Folder + "DoseHist.eps", format='eps', bbox_inches="tight")

# ----------------- Primaries weighted for Dose ----------------------------
plt.figure(4)
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
plt.savefig(Path + Folder + "/img/" + Folder + "DoseWeight.eps", format='eps', bbox_inches="tight")

# ----------------- Primaries weighted for Dose ZOOM ----------------------------
plt.figure(5)
i = 0
for Thick in ThickList:
    plt.hist(Data[i][0], weights=Data[i][1], range=(0, Max), bins=100, label=str(int(Thick)) + "mm Al", color=Colours[i])
    i += 1

plt.yscale("log")
plt.grid(which='major')
plt.title("Histogram of dose depositions per primary energy bin\n" + Folder)
plt.xlabel("Primary particle energy [MeV]")
plt.ylabel("Dose deposited per primary energy bin [MeV]")
plt.legend()
# plt.show()
plt.savefig(Path + Folder + "/img/" + Folder + "DoseWeightZOOM.eps", format='eps', bbox_inches="tight")

# ----------------- Secondaries Histogram ----------------------------
plt.figure(6)
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
plt.savefig(Path + Folder + "/img/" + Folder + "SecondaryWeight.eps", format='eps', bbox_inches="tight")

# srun --mem=50G --time=00:15:00 python ShieldedHistograms.py