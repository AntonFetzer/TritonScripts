import uproot
import numpy as np
import matplotlib.pyplot as plt
from MeVtokRad_2D import MeVtokRad_2D

# Path = "/scratch/work/fetzera1/2DGradient/2MaterialGradient/"
Path = "/home/anton/Desktop/triton_work/Gradient/2MaterialGradient/"
File = "2materialsgradient1e8protonalta.root"  # ---------------------------------------------------------------------
# NORM_FACTOR_SPECTRUM = 5.886798E+14  # Electron500keV
NORM_FACTOR_SPECTRUM = 3.381390E+11  # Protons10MeV
Npart = 1e8 / 100

file = uproot.open(Path + File)
print("Read in: " + Path + File)

tree = file["Detector Data 0"]

Edep = np.zeros(99)
x = np.linspace(0, 98, num=99, dtype=int)

print(x)

for i in x:
    Edep[i] = np.sum(tree["Sivol_" + str(i) + "_Edep_MeV"].array(library="np"))
    Edep[i] = MeVtokRad_2D(Edep[i], NORM_FACTOR_SPECTRUM, Npart)

plt.plot(x+1, Edep, '.')
plt.title("Dose deposited in 0.5 mm Si behind 2.5g/cm2 of Ta-Al shielding \n Tantalum on top of Aluminium")  # ---------
plt.xlabel("Aluminium mass ratio [%]")
plt.ylabel("Deposited ionising dose [krad]")
# plt.yscale("log")
plt.grid(which='both')
plt.savefig(Path + "/Img/protonAlTa.eps", format='eps')  # -------------------------------------------------------------

# srun --mem=50G --time=00:10:00 python 2DGradient.py
