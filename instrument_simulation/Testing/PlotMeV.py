import numpy as np
import matplotlib.pyplot as plt

Path = "/l/triton_work/Gradient/2Material1-5gcm2/TEST/root/gradient-test-2e9electron.txt"

# Data = np.loadtxt(Path)[:99]
Data = np.genfromtxt(Path, delimiter=',', dtype=None, encoding='ASCII')

Keys = []
Edep = []
Err = []

for row in Data:
    Keys.append(row[0])
    Edep.append(row[1])
    Err.append(row[2])

x = list(range(1, 100))
print(x, Edep, Err)

plt.errorbar(x, Edep, Err, fmt=".", capsize=10)

plt.grid(which='both')
plt.title("Total Dose behind 5mm of aluminium shielding")
plt.xlabel("Number of layers")
plt.ylabel("Total Ionising Dose [MeV]")
plt.savefig(Path.split(".")[0] + ".eps", format='eps')
plt.show()
