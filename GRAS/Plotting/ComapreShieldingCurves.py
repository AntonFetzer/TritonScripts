from GRAS.Dependencies.TotalDose import totalDose
import numpy as np
import matplotlib.pyplot as plt
from os import path


Path = "/l/triton_work/ShieldingCurves/Carrington/"
res_suffix = "/Res/"

Labels = [
'CarringtonElectronDiffPowTabelated-10mm',
'CarringtonElectronINTEGRALPowTabelated-10mm'
]


Paths = [Path + label + res_suffix for label in Labels]

Colours = ['C0', 'C1', 'C2', 'C8', 'C3', 'C9', 'C7', 'k', 'C4', 'C5', 'C6']

ShieldingCurves = []
for i, path in enumerate(Paths):
   ShieldingCurves.append(totalDose(path))

NumTiles = np.shape(ShieldingCurves[0]['dose'])[0]


## Plot shielding curves with error bars
plt.figure(1)

x = np.linspace(0, 10, num=NumTiles, endpoint=True)

for i, Curve in enumerate(ShieldingCurves):
    plt.errorbar(x, Curve['dose'], Curve['error'], fmt='', markersize=5, capsize=5, label=Labels[i], color=Colours[i], linestyle='')

####### Plot 10kRad line #########
CriticalDose = [10 for i in x]
plt.plot(x, CriticalDose, color='k', linewidth=2, label='10 krad')
CriticalDose = [100 for i in x]
plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='100 krad')

plt.grid()
plt.yscale("log")
# plt.xscale("log")
plt.title("Ionising Dose behind shielding of varying thickness")
plt.xlabel("Aluminium Shielding thickness [mm]")
plt.ylabel("Ionising Dose [krad/s]")
plt.legend()

# plt.savefig(Path + "/ShieldingCurves.pdf", format='pdf', bbox_inches="tight")


## Plot absolute deviation from the first dataset
plt.figure(2)

for i, Curve in enumerate(ShieldingCurves[1:]):  # Skip the first dataset as it's the reference
    DataDiff = ShieldingCurves[0]['dose'] - Curve['dose']  # Calculate the difference from the first dataset and plot it
    plt.plot(x, DataDiff, '.', label=Labels[i+1], color=Colours[i+1])

plt.title("Difference from the first dataset")
plt.xlabel("Aluminium Shielding thickness [mm]")
plt.ylabel("Difference in Ionising Dose [krad/s]")
plt.legend()
plt.grid()

#plt.savefig(Path + "/Difference.pdf", format='pdf', bbox_inches="tight")

## Plot relative deviation from the first dataset
plt.figure(3)

for i, Curve in enumerate(ShieldingCurves[1:]):  # Skip the first dataset as it's the reference
    DataRatio = ( ( Curve['dose'] / ShieldingCurves[0]['dose'] ) - 1 ) * 100  # Calculate the relative deviation from the first dataset and plot it
    plt.plot(x, DataRatio, '.', label=Labels[i+1], color=Colours[i+1])

# plt.yscale("log")
plt.title("Relative Deviation from the first dataset [%]")
plt.xlabel("Aluminium Shielding Depth [g/cm2]")
plt.ylabel("Relative Deviation")
plt.legend()
plt.grid()

# plt.savefig(Path + "/RelativeDeviation.pdf", format='pdf', bbox_inches="tight")
plt.show()

