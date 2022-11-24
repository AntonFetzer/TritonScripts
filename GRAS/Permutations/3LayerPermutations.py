from GRAS.Dependencies.TotalKRadGras import totalkRadGras
import numpy as np
from uncertainties import ufloat

Materials = ["G4_Al", "Al-7075", "G4_POLYETHYLENE", "G4_KEVLAR", "G4_Pb", "G4_W", "G4_STAINLESS-STEEL", "CarbonFibre", "FR4", "G4_Ta", "G4_TEFLON", "Ti-6AL-4V"]

#Densities = [2, 0.534, 1.55, 2.33, 1.74, 2.37, 2.699, 4.54, 11.35, 19.3]

print("Number of Names:", len(Materials))
#print("Number of Densities:", len(Densities))

#VolumesStr = ''
#for i, Name in enumerate(Names):
#    VolumesStr += '        <volume name ="ShieldVol_' + str(i) + '">\n            <materialref ref="' + Name + '"/>\n            <solidref ref="Shield_' + str(i) + '"/>\n        </volume>\n\n'
#print(VolumesStr)

#for x in Densities:
#    print(x)

Path = "/home/anton/Desktop/triton_work/Permutations/3Layer/Res/"

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
            i = i1 * NumMat * NumMat + i2 * NumMat + i3
            print(i, Materials[i1], Materials[i2], Materials[i3], ufloat(Electrons[0][i], Electrons[1][i]), ufloat(Protons[0][i], Protons[1][i]), ufloat(Total[0][i], Total[1][i]))
