import matplotlib.pyplot as plt
from GRAS.Dependencies.TotalGRASLETHistos import totalGRASLETHistos
import os
from natsort import natsorted
import numpy as np

Paths = ["/home/anton/Desktop/triton_work/LET/Carrington-SEP-Plus2Sigma-Int-With0/16mm/Res/",
         "/home/anton/Desktop/triton_work/LET/Carrington-SEP-Expected-Int-With0/16mm/Res/",
         "/home/anton/Desktop/triton_work/LET/Carrington-SEP-Minus2Sigma-Int-With0/16mm/Res/"]

Labels = ["Carington SEP +2 Sigma SEP",
          "Carington SEP Expected SEP",
          "Carington SEP -2 Sigma SEP"]

#Colours = ['C3', 'C1', 'C8', 'C2', 'C9', 'C0', 'C7']
Colours = ['C1', 'C0', 'C2']

lowerID = 0
upperID = 1
meanID = 2
valueID = 3
errorID = 4
entriesID = 5

LETHist = []
EffHist = []

for path in Paths:
    Temp = totalGRASLETHistos(path, "")
    LETHist.append(Temp[0])
    EffHist.append(Temp[1])

Num = len(Paths)

plt.figure(0)
for i in range(Num):
    plt.bar(LETHist[i][:, lowerID], LETHist[i][:, entriesID], width=LETHist[i][:, upperID] - LETHist[i][:, lowerID],
            align='edge', alpha=0.5, color=Colours[i])
    plt.step(LETHist[i][:, lowerID], LETHist[i][:, entriesID], where='post', label=Labels[i], color=Colours[i])
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("LET Histogram Carrington SEP vs. 16mm Aluminium")
plt.xlabel("LET [MeV/cm]")
plt.ylabel("Number of entries per LET bin")
plt.legend()
plt.savefig("/home/anton/Desktop/triton_work/LET/Plots/LETentries.pdf", format='pdf', bbox_inches="tight")

plt.figure(1)
for i in range(Num):
    plt.bar(LETHist[i][:, lowerID], LETHist[i][:, valueID], width=LETHist[i][:, upperID] - LETHist[i][:, lowerID],
            align='edge', alpha=0.5, color=Colours[i])
    plt.step(LETHist[i][:, lowerID], LETHist[i][:, valueID], where='post', label=Labels[i], color=Colours[i])
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("LET Histogram Carrington SEP vs. 16mm Aluminium")
plt.xlabel("LET [MeV/cm]")
plt.ylabel("Rate per LET bin [s-1]")
plt.legend()
plt.savefig("/home/anton/Desktop/triton_work/LET/Plots/LETvalues.pdf", format='pdf', bbox_inches="tight")
'''
NumberEntriesEffHist = sum(EffHist[0][:, entriesID])

plt.figure(2)
for i in range(Num):
    plt.bar(EffHist[i][:, lowerID], EffHist[i][:, entriesID], width=EffHist[i][:, upperID] - EffHist[i][:, lowerID],
            align='edge', label=Labels[i], alpha=0.3)
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("EffLET Histogram by entries")
plt.xlabel("EffLET [MeV/cm]")
plt.ylabel("Number of entries per EffLET bin")
plt.legend()
#plt.savefig("/home/anton/Desktop/triton_work/LET/EFFentries.eps", format='eps', bbox_inches="tight")

plt.figure(3)
for i in range(Num):
    plt.bar(EffHist[i][:, lowerID], EffHist[i][:, valueID], width=EffHist[i][:, upperID] - EffHist[i][:, lowerID],
            align='edge', label=Labels[i], alpha=0.3)
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("EffLET Histogram by values")
plt.xlabel("EffLET [MeV/cm]")
plt.ylabel("Rate per LET bin [s-1]")
plt.legend()
#plt.savefig("/home/anton/Desktop/triton_work/LET/EFFvalues.eps", format='eps', bbox_inches="tight")
'''
# plt.show()
