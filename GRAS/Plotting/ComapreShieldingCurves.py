from GRAS.Dependencies.TotalKRadGras import totalkRadGras
import numpy as np
import matplotlib.pyplot as plt
from os import path

Path = "/l/triton_work/ShieldingCurves/MultilayerPaper/"
res_suffix = "/Res/"

Labels = [
    "AE9-GTO",
    "AP9-GTO",
    "SolarGTO-Protons",
    "SolarGTO-He",
    "SolarGTO-C",
    "SolarGTO-O",
    "SolarGTO-Ne",
    "SolarGTO-Mg",
    "SolarGTO-Si",
    "SolarGTO-Fe",
    # "ISO-GTO-Protons",
    # "ISO-GTO-He",
    # "ISO-GTO-C",
    # "ISO-GTO-N",
    # "ISO-GTO-O",
    # "ISO-GTO-Mg",
    # "ISO-GTO-Si",
    # "ISO-GTO-Fe",
    "ISO-GTO-Protons-mission",
    "ISO-GTO-He-mission",
    "ISO-GTO-C-mission",
    "ISO-GTO-N-mission",
    "ISO-GTO-O-mission",
    "ISO-GTO-Mg-mission",
    "ISO-GTO-Si-mission",
    "ISO-GTO-Fe-mission",
]

Paths = [Path + label + res_suffix for label in Labels]

Colours = ['C0', 'C1', 'C2', 'C8', 'C3', 'C9', 'C7', 'k', 'C4', 'C5', 'C6',
           'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'white',
           '#800000', '#808000', '#800080', '#008080',  # some HTML hex colors
           (0.1, 0.2, 0.3), (0.4, 0.5, 0.6), (0.7, 0.8, 0.9),  # some RGB colors
           (0.1, 0.2, 0.3, 0.4), (0.5, 0.6, 0.7, 0.8)]  # some RGBA colors

Data = []
for i, path in enumerate(Paths):
    Data.append(totalkRadGras(path, ""))

NumTiles = np.shape(Data[0])[1]

plt.figure(1)

x = np.linspace(0, 2.5, num=NumTiles, endpoint=True)

# plt.fill_between(x, Data[0][0], Data[1][0], color='C1', alpha=0.5)
# plt.fill_between(x, Data[1][0], Data[2][0], color='C2', alpha=0.5)

seconds_in_a_month = 60 * 60 * 24 * 30.44  # number of seconds in a month


for i, data in enumerate(Data):
    if "mission" in Labels[i]:
        plt.errorbar(x, data[0] / seconds_in_a_month, data[1] / seconds_in_a_month, fmt='', markersize=5, capsize=5, label=Labels[i], color=Colours[i], linestyle='')
    else:
        plt.errorbar(x, data[0], data[1], fmt='', markersize=5, capsize=5, label=Labels[i], color=Colours[i], linestyle='')

####### Plot 10kRad line #########
CriticalDose = [10 for i in x]
#plt.plot(x, CriticalDose, color='k', linewidth=2, label='10 krad')
# CriticalDose = [100 for i in x]
# plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='100 krad')

# plt.ylim(3e-3, 1e3)
plt.xlim(0.25, 2.5)
plt.grid()
plt.yscale("log")
# plt.xscale("log")
plt.title("Ionising Dose per month behind shielding of varying thickness")
plt.xlabel("Aluminium Shielding Depth [g/cm2]")
plt.ylabel("Ionising Dose per month [krad]")
plt.legend()

plt.savefig(Path + "/ShieldingCurves.pdf", format='pdf', bbox_inches="tight")

'''
plt.figure(2)
DataDiff = abs(Data[0] - Data[1])

plt.errorbar(x, DataDiff[0], DataDiff[1], fmt=' ', capsize=5, label="DataDiff")
plt.yscale("log")
plt.title("DataDiff")
plt.xlabel("Entry#")
plt.ylabel("Difference in Ionising Dose [krad]")
plt.legend()

plt.figure(3)
DataRatio = Data[0] / Data[1]

plt.scatter(x, DataRatio[0], label="DataRatio")
plt.yscale("log")
plt.title("DataRatio")
plt.xlabel("Entry#")
plt.ylabel("Difference in Ionising Dose [krad]")
plt.legend()

plt.show()

print(np.average(DataRatio[0]))
'''
