import matplotlib.pyplot as plt
from Dependencies.TotalDoseHistograms import totalDoseHistograms
import numpy as np
from uncertainties import ufloat

Paths = ["/l/triton_work/LunarBackscatter/LunarSEP/0mm/Res/",
         "/l/triton_work/Histograms/LunarSEP/Res/"]


Labels = ["SEP Backscatter",
          "SEP"]

Colours = ['C1', 'C0', 'C2', 'C8', 'C3', 'C7']

DoseHists = []
PrimaryHists = []
for path in Paths:
    DoseHist, PrimaryHist = totalDoseHistograms(path)
    DoseHists.append(DoseHist)
    PrimaryHists.append(PrimaryHist)

plt.figure(1)
for i, PrimaryHist in enumerate(PrimaryHists):
    TotalDose = np.sum(PrimaryHist['value'])
    # Bins of the same run share the same primaries, so their errors are
    # correlated rather than independent; quadrature over bins is the
    # conventional estimate here and at worst understates the total error.
    TotalError = np.sqrt(np.sum(PrimaryHist['error'] ** 2))
    plt.bar(PrimaryHist['lower'], PrimaryHist['value'], width=PrimaryHist['upper'] - PrimaryHist['lower'],
            align='edge', alpha=0.5, color=Colours[i])
    plt.step(PrimaryHist['lower'], PrimaryHist['value'], where='post', label=Labels[i] + " " + str(ufloat(TotalDose, TotalError)) + " krad/month", color=Colours[i])

plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("Ionising dose VS primary kinetic energy in unshielded Si")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Relative ionising Dose [a.u.]")
plt.legend()
plt.savefig("/l/triton_work/LunarBackscatter/ComparisonPlots/HistogramComparisonDose.pdf", format='pdf', bbox_inches="tight")

plt.figure(2)

for i, PrimaryHist in enumerate(PrimaryHists):
    plt.bar(PrimaryHist['lower'], PrimaryHist['entries'], width=PrimaryHist['upper'] - PrimaryHist['lower'],
            align='edge', alpha=0.5, color=Colours[i])
    plt.step(PrimaryHist['lower'], PrimaryHist['entries'], where='post', label=Labels[i], color=Colours[i])

plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("Particle count VS primary kinetic energy")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Number of entries")
plt.legend()

plt.savefig("/l/triton_work/LunarBackscatter/ComparisonPlots/HistogramComparisonCounts.pdf", format='pdf', bbox_inches="tight")

#plt.show()
