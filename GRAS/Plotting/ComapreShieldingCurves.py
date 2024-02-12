from GRAS.Dependencies.TotalKRadGras import totalkRadGras
import numpy as np
import matplotlib.pyplot as plt
from os import path


Path = "/l/triton_work/ShieldingCurves/Cobalt-60/"
res_suffix = "/Res/"

Labels = [
'Co-60 10mm Al',
'Co-60 10mm Al + 1m Air',
'Co-60 10mm Al + 1m Air + 1-5mm Pb + 0-7mm Al',
# 'Co-60 10mm Al + 1mAir + 5mm Plex',
'Co-60 10mm Al + 1mmPlex + small Gap',
'Co-60 10mm Al + 2mmPlex + small Gap',
'Co-60 10mm Al + 2mmPlex + small Gap + 2mmPlex behind',
# 'Co-60 10mm Al + 5mmPlex + large Gap',
'Co-60 10mm Al + 5mmPlex + small Gap',
]


Paths = [Path + label + res_suffix for label in Labels]

Colours = ['C0', 'C1', 'C2', 'C8', 'C3', 'C9', 'C7', 'k', 'C4', 'C5', 'C6',
           'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'white',
           '#800000', '#808000', '#800080', '#008080',  # some HTML hex colors
           (0.1, 0.2, 0.3), (0.4, 0.5, 0.6), (0.7, 0.8, 0.9),  # some RGB colors
           (0.1, 0.2, 0.3, 0.4), (0.5, 0.6, 0.7, 0.8),  # some RGBA colors
            'C0', 'C1', 'C2', 'C8', 'C3', 'C9', 'C7', 'k', 'C4', 'C5', 'C6',    
           'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'white']

Data = []
for i, path in enumerate(Paths):
   Data.append(totalkRadGras(path, "")*6e10)

NumTiles = np.shape(Data[0])[1]

plt.figure(1)

x = np.linspace(0, 10, num=NumTiles, endpoint=True)

for i, data in enumerate(Data):
    plt.errorbar(x, data[0], data[1], fmt='', markersize=5, capsize=5, label=Labels[i], color=Colours[i], linestyle='')

####### Plot 10kRad line #########
#CriticalDose = [10 for i in x]
#plt.plot(x, CriticalDose, color='k', linewidth=2, label='10 krad')
# CriticalDose = [100 for i in x]
# plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='100 krad')

plt.ylim(0.6, 1.05)
plt.xlim(-0.1, 5)
plt.grid()
# plt.yscale("log")
# plt.xscale("log")
plt.title("Ionising Dose behind shielding of varying thickness")
plt.xlabel("Aluminium Shielding thickness [mm]")
plt.ylabel("Ionising Dose [a.u.]")
plt.legend()


plt.savefig(Path + "/ShieldingCurves.pdf", format='pdf', bbox_inches="tight")
# plt.show()
'''
plt.figure(2)

# Calculate the difference with the first dataset and plot it
for i, data in enumerate(Data[1:]):  # Skip the first dataset as it's the reference
    DataDiff = Data[0][0] - data[0]
    plt.plot(x, DataDiff, '-', label=Labels[i+1], color=Colours[i+1])

plt.yscale("log")
plt.title("Difference from the first dataset")
plt.xlabel("Aluminium Shielding Depth [g/cm2]")
plt.ylabel("Difference in Ionising Dose [a.u.]")
plt.legend()
plt.grid()
plt.savefig(Path + "/Difference.pdf", format='pdf', bbox_inches="tight")

plt.figure(3)

# Calculate the relative deviation with the first dataset and plot it
for i, data in enumerate(Data[1:]):  # Skip the first dataset as it's the reference
    DataRatio = ( ( data[0] / Data[0][0] ) - 1 ) * 100  # Dividing by the first dataset
    plt.plot(x, DataRatio, '-', label=Labels[i+1], color=Colours[i+1])

# plt.yscale("log")
plt.title("Relative Deviation from the first dataset [%]")
plt.xlabel("Aluminium Shielding Depth [g/cm2]")
plt.ylabel("Relative Deviation")
plt.legend()
plt.grid()
plt.savefig(Path + "/RelativeDeviation.pdf", format='pdf', bbox_inches="tight")

# plt.show()
'''
