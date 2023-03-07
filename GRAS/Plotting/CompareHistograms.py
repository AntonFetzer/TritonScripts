import matplotlib.pyplot as plt
from GRAS.Dependencies.TotalGRASHistos import totalGRASHistos
import os
from natsort import natsorted
import numpy as np

Paths = [#"/home/anton/Desktop/triton_work/Histograms/AP910MeV/Res/",
#"/home/anton/Desktop/triton_work/Histograms/AP910MeV-DIffwith0/Res/",
#"/home/anton/Desktop/triton_work/Histograms/AP910MeV-Integral/Res/",
#"/home/anton/Desktop/triton_work/Histograms/AP910MeV-Integral-With0/Res/"
         #"/home/anton/Desktop/triton_work/Histograms/Carrington-SEP-Expected-Diff/Res/",
         #"/home/anton/Desktop/triton_work/Histograms/Carrington-SEP-Expected-Int/Res/",
         "/home/anton/Desktop/triton_work/Histograms/Carrington-SEP-Plus2Sigma-Int-With0/Res/",
         "/home/anton/Desktop/triton_work/Histograms/Carrington-SEP-Expected-Int-With0/Res/",
         "/home/anton/Desktop/triton_work/Histograms/Carrington-SEP-Minus2Sigma-Int-With0/Res/",
         "/home/anton/Desktop/triton_work/Histograms/AP910MeV/Res/",
         "/home/anton/Desktop/triton_work/Histograms/SEP2003-INTEGRAL-FluxBasedOnFluenceDividedBy24h/Res/"
         ]


Labels = [#"AP910MeV",
#"AP910MeV-DIffwith0",
#"AP910MeV-Integral",
#"AP910MeV-Integral-With0"
          #"Car SEP Expected Diff",
          #"Car SEP Expected Int",
          "Carrington SEP +2 Sigma",
          "Carrington SEP Expected",
          "Carrington SEP -2 Sigma",
          "AP9 GTO trapped protons",
          "2003 SEP"
          ]

Colours = ['C1', 'C0', 'C2', 'C3', 'C8']

lowerID = 0
upperID = 1
meanID = 2
valueID = 3
errorID = 4
entriesID = 5

DoseHists = []
PrimaryHists = []
for path in Paths:
    Temp = totalGRASHistos(path, "")
    DoseHists.append(Temp[0])
    PrimaryHists.append(Temp[1])

plt.figure(1)
for i, PrimaryHist in enumerate(PrimaryHists):
    plt.bar(PrimaryHist[:, lowerID], PrimaryHist[:, valueID], width=PrimaryHist[:, upperID] - PrimaryHist[:, lowerID],
            align='edge', alpha=0.5, color=Colours[i])
    plt.step(PrimaryHist[:, lowerID], PrimaryHist[:, valueID], where='post', label=Labels[i], color=Colours[i])

#plt.ylim(1e-6, 10)
plt.xlim(8, 300)
plt.ylim(5e-4, 2e2)
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("Ionising dose VS primary kinetic energy in unshielded Si")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Relative ionising Dose [a.u.]")
plt.legend()
plt.savefig("/home/anton/Desktop/triton_work/Histograms/HistogramComparisonDose.pdf", format='pdf', bbox_inches="tight")

plt.figure(2)

for i, PrimaryHist in enumerate(PrimaryHists):
    plt.bar(PrimaryHist[:, lowerID], PrimaryHist[:, entriesID], width=PrimaryHist[:, upperID] - PrimaryHist[:, lowerID],
            align='edge', alpha=0.5, color=Colours[i])
    plt.step(PrimaryHist[:, lowerID], PrimaryHist[:, entriesID], where='post', label=Labels[i], color=Colours[i])

plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("Particle count VS primary kinetic energy")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Number of entries")
plt.legend()

plt.savefig("/home/anton/Desktop/triton_work/Histograms/HistogramComparisonCounts.pdf", format='pdf', bbox_inches="tight")

'''
plt.figure(3)
plt.bar(DoseHist1[:, lowerID], DoseHist1[:, entriesID], width=DoseHist1[:, upperID] - DoseHist1[:, lowerID], align='edge', label=Data1, alpha=0.5)
plt.bar(DoseHist2[:, lowerID], DoseHist2[:, entriesID], width=DoseHist2[:, upperID] - DoseHist2[:, lowerID], align='edge', label=Data2, alpha=0.5)
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("Dose Depositions Histogram")
plt.xlabel("Dose [krad per Month]")
plt.ylabel("Number of entries per dose bin")
plt.legend()

plt.figure(4)
plt.bar(DoseHist1[:, lowerID], DoseHist1[:, valueID], width=DoseHist1[:, upperID] - DoseHist1[:, lowerID], align='edge', label=Data1, alpha=0.5)
plt.bar(DoseHist2[:, lowerID], DoseHist2[:, valueID], width=DoseHist2[:, upperID] - DoseHist2[:, lowerID], align='edge', label=Data2, alpha=0.5)
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("Dose Depositions Histogram")
plt.xlabel("Dose [krad per Month]")
plt.ylabel("???")
plt.legend()
'''
#plt.show()
