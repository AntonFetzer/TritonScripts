import matplotlib.pyplot as plt
from GRAS.Dependencies.TotalGRASLETHistos import totalGRASLETHistos
import os
from natsort import natsorted
import numpy as np
from uncertainties import ufloat


absolute_path = os.path.abspath(__file__)
print("Full path: " + absolute_path)
print("Directory Path: " + os.path.dirname(absolute_path))

Path = "/l/triton_work/LET/A9-LEO/AP9Mission/"
Title = "AP9Mission LEO"

lowerID = 0
upperID = 1
meanID = 2
valueID = 3
errorID = 4
entriesID = 5

Folders = [f for f in os.listdir(Path) if f.endswith('mm')]
Folders = natsorted(Folders)

print(Folders)
Num = len(Folders)

ThickList = []

for folder in Folders:
    ThickList.append(int(folder.split("m")[0]))

print(ThickList)

LETHist = []
EffHist = []

for i in range(Num):
    Data = totalGRASLETHistos(Path + Folders[i] + "/Res/", "")
    # print("Shape of Data:", np.shape(Data))
    TotalNumberEntries = sum(Data[0][:, entriesID])
    print("Total number of particles = " + f"{TotalNumberEntries:.5}")
    print("Number of Data points for", str(ThickList[i]) + "mm:", "{:,}".format(len(Data[0])))

    LETHist.append(Data[0])
    EffHist.append(Data[1])

print("Shape of LETHist:", np.shape(LETHist))
# [Thickness][Bin][Variable]

Colours = ['C3', 'C1', 'C8', 'C2', 'C9', 'C0', 'C7']
#          red, orange,yellow,green,cyan,blue, grey
# Colours = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
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

#NumberEntriesLETHist = sum(LETHist[0][:, entriesID])

C = 2330  # to convert from MeV/cm to Mev cm2 mg-1
for i in range(Num):
    LETHist[i][:, lowerID] = LETHist[i][:, lowerID] / C
    LETHist[i][:, upperID] = LETHist[i][:, upperID] / C
    LETHist[i][:, meanID] = LETHist[i][:, meanID] / C


plt.figure(0)
for i in range(Num):
    plt.bar(LETHist[i][:, lowerID], LETHist[i][:, entriesID], width=LETHist[i][:, upperID] - LETHist[i][:, lowerID],
            align='edge', alpha=0.5, color=Colours[i])
    plt.step(LETHist[i][:, lowerID], LETHist[i][:, entriesID], where='post', label=str(ThickList[i]) + "mm Al", color=Colours[i])

plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title(Title)
plt.xlabel("LET [MeV cm2 mg-1]")
plt.ylabel("Number of entries per LET bin")
plt.legend()
plt.savefig(Path + "StackedPlots/LETentries.pdf", format='pdf', bbox_inches="tight")

#TotalLETbyValues = []
#TotalLETERRORbyValues = []
#  " + str(ufloat(TotalLETbyValues[i], TotalLETERRORbyValues[i])) + " MeV/cm total LET"

plt.figure(1)
for i in range(Num):
    #TotalLETbyValues.append(sum(LETHist[i][:, meanID] * LETHist[i][:, valueID]))
    #TotalLETERRORbyValues.append(sum(LETHist[i][:, meanID] * LETHist[i][:, errorID]))
    plt.bar(LETHist[i][:, lowerID], LETHist[i][:, valueID], width=LETHist[i][:, upperID] - LETHist[i][:, lowerID],
            align='edge', alpha=0.5, color=Colours[i])
    plt.step(LETHist[i][:, lowerID], LETHist[i][:, valueID], where='post', label=str(ThickList[i]) + "mm Al", color=Colours[i])
    
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title(Title)
plt.xlabel("LET [MeV cm2 mg-1]")
plt.ylabel("Rate per LET bin [cm-2 s-1]")
plt.legend()
plt.savefig(Path + "StackedPlots/LETvalues.pdf", format='pdf', bbox_inches="tight")

#plt.show()

'''
NumberEntriesEffHist = sum(EffHist[0][:, entriesID])

plt.figure(2)
for i in range(Num):
    plt.bar(EffHist[i][:, lowerID], EffHist[i][:, entriesID], width=EffHist[i][:, upperID] - EffHist[i][:, lowerID],
            align='edge', alpha=0.5, color=Colours[i])
    plt.step(EffHist[i][:, lowerID], EffHist[i][:, valueID], where='post', label=str(ThickList[i]) + "mm Al", color=Colours[i])
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("EffLET Histogram" + f"{NumberEntriesEffHist:.2}" + " events each")
plt.xlabel("EffLET [MeV/cm]")
plt.ylabel("Number of entries per EffLET bin")
plt.legend()
plt.savefig(Path + "StackedPlots/EFFentries.pdf", format='pdf', bbox_inches="tight")

plt.figure(3)
for i in range(Num):
    plt.bar(EffHist[i][:, lowerID], EffHist[i][:, valueID], width=EffHist[i][:, upperID] - EffHist[i][:, lowerID],
            align='edge', label=str(ThickList[i]) + "mm Al", alpha=0.5, color=Colours[i])
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("EffLET Histogram " + f"{NumberEntriesEffHist:.2}" + " events each")
plt.xlabel("EffLET [MeV/cm]")
plt.ylabel("Rate per LET bin [cm-2 s-1]")
plt.legend()
#plt.xlim(10, 1e6)
plt.savefig(Path + "StackedPlots/EFFvalues.pdf", format='pdf', bbox_inches="tight")

#plt.show()
'''