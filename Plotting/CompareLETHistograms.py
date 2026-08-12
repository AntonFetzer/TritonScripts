import matplotlib.pyplot as plt
from Dependencies.TotalLETHistos import totalLETHistos
import os
from natsort import natsorted
import numpy as np

Thick = "16mm"

Paths = ["/l/triton_work/LET/Carrington/Carrington-SEP-Expected-Int-With0/" + Thick + "/Res/",
         "/l/triton_work/LET/Carrington/Carrington-SEP-Minus2Sigma-Int-With0/" + Thick + "/Res/",
         "/l/triton_work/LET/Carrington/Carrington-SEP-Plus2Sigma-Int-With0/" + Thick + "/Res/",
         "/l/triton_work/LET/Carrington/SEP2003-INTEGRAL-FluxBasedOnFluenceDividedBy24h/" + Thick + "/Res/",
         "/l/triton_work/LET/A9-GTO/AE9Mission/" + Thick + "/Res/",
         "/l/triton_work/LET/A9-GTO/AP9Mission/" + Thick + "/Res/",
         "/l/triton_work/LET/A9-LEO/AE9Mission/" + Thick + "/Res/",
         "/l/triton_work/LET/A9-LEO/AP9Mission/" + Thick + "/Res/"]

Labels = ["SEP-Expected", "Minus2Sigma", "Plus2Sigma", "SEP2003", "AE9 GTO", "AP9 GTO", "AE9 LEO", "AP9 LEO"]

Colours = ['C1', 'C0', 'C2', 'C8', 'C3', 'C7', 'C9', 'b', 'k']

LETHist = []
EffHist = []

for p, path in enumerate(Paths):
    LET, Eff = totalLETHistos(path)
    LETHist.append(LET)
    EffHist.append(Eff)

Num = len(Paths)
'''
plt.figure(0)
for i in range(Num):
    plt.bar(LETHist[i]['lower'], LETHist[i]['entries'], width=LETHist[i]['upper'] - LETHist[i]['lower'],
            align='edge', alpha=0.5, color=Colours[i])
    plt.step(LETHist[i]['lower'], LETHist[i]['entries'], where='post', label=Labels[i], color=Colours[i])
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("LET Histogram Carrington SEP vs. 16mm Aluminium")
plt.xlabel("LET [MeV cm2 mg-1]")
plt.ylabel("Number of entries per LET bin")
plt.legend()
plt.savefig("/l/triton_work/LET/Plots/LETentries.pdf", format='pdf', bbox_inches="tight")
'''
plt.figure(1)
for i in range(Num):
    if Labels[i].startswith("A"):
        LETHist[i]['value'] /= 30 * 24 * 60 * 60
        LETHist[i]['error'] /= 30 * 24 * 60 * 60

    plt.bar(LETHist[i]['lower'], LETHist[i]['value'], width=LETHist[i]['upper'] - LETHist[i]['lower'],
            align='edge', alpha=0.5, color=Colours[i])
    plt.step(LETHist[i]['lower'], LETHist[i]['value'], where='post', label=Labels[i], color=Colours[i])
plt.yscale("log")
plt.xscale("log")
# plt.ylim(1e-7, 1e5)
plt.grid()
plt.title("LET Histogram Carrington SEP vs. " + Thick +" Aluminium")
plt.xlabel("LET [MeV/cm]")
plt.ylabel("Rate per LET bin [s-1]")
plt.legend()
# plt.savefig("/l/triton_work/LET/Plots/LETvalues.pdf", format='pdf', bbox_inches="tight")
# plt.savefig("/l/TritonPlots/Luna/LETHistogramComparison.svg", format='svg', bbox_inches="tight")
'''
NumberEntriesEffHist = sum(EffHist[0]['entries'])

plt.figure(2)
for i in range(Num):
    plt.bar(EffHist[i]['lower'], EffHist[i]['entries'], width=EffHist[i]['upper'] - EffHist[i]['lower'],
            align='edge', label=Labels[i], alpha=0.3)
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("EffLET Histogram by entries")
plt.xlabel("EffLET [MeV cm2 mg-1]")
plt.ylabel("Number of entries per EffLET bin")
plt.legend()
#plt.savefig("/l/triton_work/LET/EFFentries.eps", format='eps', bbox_inches="tight")

plt.figure(3)
for i in range(Num):
    plt.bar(EffHist[i]['lower'], EffHist[i]['value'], width=EffHist[i]['upper'] - EffHist[i]['lower'],
            align='edge', label=Labels[i], alpha=0.3)
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("EffLET Histogram by values")
plt.xlabel("EffLET [MeV cm2 mg-1]")
plt.ylabel("Rate per LET bin [s-1]")
plt.legend()
#plt.savefig("/l/triton_work/LET/EFFvalues.eps", format='eps', bbox_inches="tight")
'''
plt.show()
