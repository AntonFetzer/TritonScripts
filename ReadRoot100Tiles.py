import uproot
import numpy as np
import matplotlib.pyplot as plt

FileName = "output.root"

file = uproot.open(FileName)
print("Read in: " + FileName)

tree = file["Detector Data 0"]

Edep = []

for x in range(100):
	DoseArray = tree["Sivol_" + str(x) + "_Edep_MeV"].array()
	Edep.append(np.sum(DoseArray))

plt.plot(Edep)
plt.title("Energy absorbed in 0.5 mm Si behind Al shielding")
plt.xlabel("Aluminium Shielding thickness [mm]")
plt.ylabel("Absorbed energy [MeV]")
plt.yscale("log")
plt.grid(which='both')
#plt.show()
plt.savefig("100Tiles.eps", format='eps')

   # srun --mem=50G --time=00:15:00 python ReadRoot100Tiles.py
