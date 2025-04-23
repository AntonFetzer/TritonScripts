from GRAS.Dependencies.TotalDose import totalDose
import numpy as np
import matplotlib.pyplot as plt
from os import path
from GRAS.Dependencies.MergeTotalDose import mergeTotalDose
import csv


Path = "/l/triton_work/Shielding_Curves/Mono/"
res_suffix = "/Res/"

Labels = [
    "32MeVProtons",
    "64MeVProtons",
    "128MeVProtons"
]

Colours = ['C0', 'C1', 'C2', 'C8', 'C3', 'C9', 'C7', 'k', 'C4', 'C5', 'C6']

# Initialise the Shielding Curves dictionary to store the shielding curves
SC = {}
for label in Labels:
    SC[label] = totalDose(Path + label + res_suffix)


NumTiles = np.shape(SC[Labels[0]]['dose'])[0]

## Plot shielding curves with error bars
plt.figure(1)

x = np.linspace(0, 10, num=NumTiles, endpoint=True)

# print("Shielding curves: ", SC)

for i, label in enumerate(Labels):
   plt.errorbar(x, SC[label]['dose'], SC[label]['error'], fmt='', capsize=5, label=label, color=Colours[i], linestyle='')

####### Plot 10kRad line #########
CriticalDose = [1 for i in x]
plt.plot(x, CriticalDose, color='k', linewidth=2, label='1 krad')
CriticalDose = [10 for i in x]
plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='10 krad')

plt.grid()
plt.yscale("log")
# plt.xscale("log")
plt.title("Ionising Dose behind shielding of varying thickness")
plt.xlabel("Aluminium Thickness [mm]")
plt.ylabel("Ionising Dose [krad/1E11 protons]")
plt.legend()
# plt.ylim(0.1, 100)
# plt.xlim(0, 2.5)

plt.savefig(Path + "/ShieldingCurves.pdf", format='pdf', bbox_inches="tight")