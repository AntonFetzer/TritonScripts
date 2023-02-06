import matplotlib.pyplot as plt
from GRAS.Dependencies.TotalGRASHistos import totalGRASHistos
import os
from natsort import natsorted
import numpy as np

Paths = [
    "/home/anton/Desktop/triton_work/Histograms/CarringtonSEPTwoPointPowerLaw/Res/",
    "/home/anton/Desktop/triton_work/Histograms/CarringtonSEP-AP9lowEndextrapolation/Res/",
    "/home/anton/Desktop/triton_work/Histograms/CarringtonSEP-AP9highEndextrapolation/Res/"]

Labels = ["CaringtonSEP",
          "Low-End extrapolation",
          "High-End extrapolation"]

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
            align='edge', alpha=0.5)
    plt.step(PrimaryHist[:, lowerID], PrimaryHist[:, valueID], where='post', label=Labels[i])

plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("Ionising dose VS primary kinetic energy")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Relative ionising Dose [a.u.]")
plt.legend()
plt.savefig("/home/anton/Desktop/triton_work/Histograms/HistogramComparison.pdf", format='pdf', bbox_inches="tight")
'''
plt.figure(2)

plt.bar(PrimaryHist1[:, lowerID], PrimaryHist1[:, entriesID], width=PrimaryHist1[:, upperID] - PrimaryHist1[:, lowerID], align='edge', label=Data1, alpha=0.5)
plt.bar(PrimaryHist2[:, lowerID], PrimaryHist2[:, entriesID], width=PrimaryHist2[:, upperID] - PrimaryHist2[:, lowerID], align='edge', label=Data2, alpha=0.5)

plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("Particle count VS primary kinetic energy")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Number of entries")
plt.legend()

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
plt.show()
