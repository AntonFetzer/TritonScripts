from GRAS.Dependencies.TotalKRadGras import totalkRadGras
import numpy as np
from uncertainties import ufloat

Materials = ["Al-7075", "G4_POLYETHYLENE", "G4_W", "FR4"]


print("Number of Names:", len(Materials))

Path = "/home/anton/Desktop/triton_work/Permutations/5Layer/Res/"

Electrons = totalkRadGras(Path, "Elec")
Protons = totalkRadGras(Path, "Prot")

print("Electrons Shape:", np.shape(Electrons))

Total = Electrons + Protons
Total[1] = np.sqrt(Electrons[1] * Electrons[1] + Protons[1] * Protons[1])
# print(Total)

NumMat = len(Materials)

for i1 in range(NumMat):
    for i2 in range(NumMat):
        for i3 in range(NumMat):
            for i4 in range(NumMat):
                for i5 in range(NumMat):
                    i = i1 * NumMat * NumMat * NumMat * NumMat + i2 * NumMat * NumMat * NumMat + i3 * NumMat * NumMat + i4 * NumMat + i5
                    print(i+1, Materials[i1], Materials[i2], Materials[i3], Materials[i4], Materials[i5], ufloat(Electrons[0][i], Electrons[1][i]), ufloat(Protons[0][i], Protons[1][i]), ufloat(Total[0][i], Total[1][i]))