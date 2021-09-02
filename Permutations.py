import numpy as np
import matplotlib.pyplot as plt
from MeVtokRad_2D import MeVtokRad_2D

# Materials = ["Polyethylene", "Teflon", "Aluminium", "Titanium", "Iron", "Copper", "Zink", "Tantalum", "Tungsten", "Lead"]


# for i1 in range(10):
# for i2 in range(10):
# for i3 in range(10):
# print(Materials[i3], Materials[i2], Materials[i1])

# print(Shields[999])

Prot = np.loadtxt("/home/anton/Desktop/triton_work/Permutations/csv/3layer2e7proton.txt")
Elec = np.loadtxt("/home/anton/Desktop/triton_work/Permutations/csv/3layer2e8electron.txt")


NORM_FACTOR_SPECTRUM_Elec = 5.886798E+14  # Elec500keV
NORM_FACTOR_SPECTRUM_Prot = 3.381390E+11  # Prots10MeV
Npart_Elec = 2e8 / 1000
Npart_Prot = 2e7 / 1000

for i in range(1000):
    print(MeVtokRad_2D(Elec[i], NORM_FACTOR_SPECTRUM_Elec, Npart_Elec))
    #print(MeVtokRad_2D(Prot[i], NORM_FACTOR_SPECTRUM_Prot, Npart_Prot))
