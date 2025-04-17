from GRAS.Dependencies.TotalDose import totalDose
import numpy as np
import matplotlib.pyplot as plt
from os import path
from GRAS.Dependencies.MergeTotalDose import mergeTotalDose
import csv


Path = "/l/triton_work/Shielding_Curves/CompareMaterials/"
res_suffix = "/Res/"

Labels = [
    "GTO-AE9-mission-W",
    "GTO-AP9-mission-W",
    "GTO-AE9-mission-Al",
    "GTO-AP9-mission-Al",
    "GTO-AE9-mission-PE",
    "GTO-AP9-mission-PE"
]

ShieldingCurves = {}
for label in Labels:
    ShieldingCurves[label] = totalDose(Path + label + res_suffix)

Colours = ['C0', 'C1', 'C2', 'C8', 'C3', 'C9', 'C7', 'k', 'C4', 'C5', 'C6']

NumTiles = np.shape(ShieldingCurves[Labels[0]]['dose'])[0]

# Merge Electron and Proton Data
ShieldingCurves["W"] = mergeTotalDose([ShieldingCurves["GTO-AE9-mission-W"], ShieldingCurves["GTO-AP9-mission-W"]])
ShieldingCurves["Al"] = mergeTotalDose([ShieldingCurves["GTO-AE9-mission-Al"], ShieldingCurves["GTO-AP9-mission-Al"]])
ShieldingCurves["PE"] = mergeTotalDose([ShieldingCurves["GTO-AE9-mission-PE"], ShieldingCurves["GTO-AP9-mission-PE"]])


## Plot shielding curves with error bars
plt.figure(1)

x = np.linspace(0, 2.5, num=NumTiles, endpoint=True)

#for i, Curve in enumerate(ShieldingCurves):
#    plt.errorbar(x, Curve['dose'], Curve['error'], fmt='', capsize=5, label=Labels[i], color=Colours[i], linestyle='')

# plt.errorbar(x, ShieldingCurves["W"]['dose'], ShieldingCurves["W"]['error'], capsize=2, label="W", color=Colours[0], linestyle='')
# plt.errorbar(x, ShieldingCurves["Al"]['dose'], ShieldingCurves["Al"]['error'], capsize=2, label="Al", color=Colours[1], linestyle='')
# plt.errorbar(x, ShieldingCurves["PE"]['dose'], ShieldingCurves["PE"]['error'], capsize=2, label="PE", color=Colours[2], linestyle='')

# plt.errorbar(x, ShieldingCurves["GTO-AE9-mission-W"]['dose'], ShieldingCurves["GTO-AE9-mission-W"]['error'], capsize=2, label="W AE9", color=Colours[0], linestyle='')
# plt.errorbar(x, ShieldingCurves["GTO-AE9-mission-Al"]['dose'], ShieldingCurves["GTO-AE9-mission-Al"]['error'], capsize=2, label="Al AE9", color=Colours[1], linestyle='')
# plt.errorbar(x, ShieldingCurves["GTO-AE9-mission-PE"]['dose'], ShieldingCurves["GTO-AE9-mission-PE"]['error'], capsize=2, label="PE AE9", color=Colours[2], linestyle='')

plt.errorbar(x, ShieldingCurves["GTO-AP9-mission-W"]['dose'], ShieldingCurves["GTO-AP9-mission-W"]['error'], capsize=2, label="W AP9", color=Colours[0], linestyle='')
plt.errorbar(x, ShieldingCurves["GTO-AP9-mission-Al"]['dose'], ShieldingCurves["GTO-AP9-mission-Al"]['error'], capsize=2, label="Al AP9", color=Colours[1], linestyle='')
plt.errorbar(x, ShieldingCurves["GTO-AP9-mission-PE"]['dose'], ShieldingCurves["GTO-AP9-mission-PE"]['error'], capsize=2, label="PE AP9", color=Colours[2], linestyle='')


####### Plot 10kRad line #########
CriticalDose = [1 for i in x]
plt.plot(x, CriticalDose, color='k', linewidth=2, label='1 krad')
CriticalDose = [10 for i in x]
plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='10 krad')

plt.grid()
plt.yscale("log")
# plt.xscale("log")
plt.title("Ionising Dose behind shielding of varying thickness")
plt.xlabel("Shielding depth [g/cm2]")
plt.ylabel("Ionising Dose [krad/month]")
plt.legend()
plt.ylim(0.1, 100)
plt.xlim(0, 2.5)

plt.savefig(Path + "/ShieldingCurves.pdf", format='pdf', bbox_inches="tight")

'''
## Plot absolute deviation from the first dataset
plt.figure(2)

for i, Curve in enumerate(ShieldingCurves[1:]):  # Skip the first dataset as it's the reference
    DataDiff = ShieldingCurves[0]['dose'] - Curve['dose']  # Calculate the difference from the first dataset and plot it
    plt.plot(x, DataDiff, '.', label=Labels[i+1], color=Colours[i+1])

plt.title("Difference from the first dataset")
plt.xlabel("Shielding depth [g/cm2]")
plt.ylabel("Difference in Ionising Dose [krad/month]")
plt.legend()
plt.grid()

plt.savefig(Path + "/Difference.pdf", format='pdf', bbox_inches="tight")

## Plot relative deviation from the first dataset
plt.figure(3)

for i, Curve in enumerate(ShieldingCurves[1:]):  # Skip the first dataset as it's the reference
    DataRatio = ( ( Curve['dose'] / ShieldingCurves[0]['dose'] ) - 1 ) * 100  # Calculate the relative deviation from the first dataset and plot it
    plt.plot(x, DataRatio, '.', label=Labels[i+1], color=Colours[i+1])

# plt.yscale("log")
plt.title("Relative Deviation from the first dataset [%]")
plt.xlabel("Shielding depth [g/cm2]")
plt.ylabel("Relative Deviation")
plt.legend()
plt.grid()

plt.savefig(Path + "/RelativeDeviation.pdf", format='pdf', bbox_inches="tight")
#plt.show()
'''

for label in Labels:
    csv_file = Path + label + "_ShieldingCurve.csv"
    with open(csv_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        # Write the header
        writer.writerow(["Shielding depth [g/cm2]", "dose [krad/month]", "error [krad/month]"])
        # Write the data
        for i in range(NumTiles):
            writer.writerow([x[i], ShieldingCurves[label]['dose'][i], ShieldingCurves[label]['error'][i]])