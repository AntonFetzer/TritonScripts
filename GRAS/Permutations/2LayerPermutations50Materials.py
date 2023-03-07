import nmmn.plots
from GRAS.Dependencies.TotalKRadGras import totalkRadGras
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from uncertainties import ufloat


Names = ["G4_lH2", "G4_He", "G4_Li", "G4_Be", "G4_B", "G4_C", "G4_lN2", "G4_lO2", "G4_F", "G4_Ne", "G4_Na", "G4_Mg",
         "G4_Al", "G4_Si", "G4_P", "G4_S", "G4_Cl", "G4_lAr", "G4_K", "G4_Ca", "G4_Sc", "G4_Ti", "G4_V", "G4_Cr",
         "G4_Mn", "G4_Fe", "G4_Co", "G4_Ni", "G4_Cu", "G4_Zn", "G4_Ga", "G4_Ge", "G4_As", "G4_Se", "G4_lBr", "G4_lKr",
         "G4_Rb", "G4_Sr", "G4_Y", "G4_Zr", "G4_Nb", "G4_Mo", "G4_Tc", "G4_Ru", "G4_Rh", "G4_Pd", "G4_Ag", "G4_Cd",
         "G4_In", "G4_Sn"]

Densities = [0.0708, 0.000166322, 0.534, 1.848, 2.37, 2, 0.807, 1.141, 0.00158029, 0.000838505, 0.971, 1.74, 2.699,
             2.33, 2.2, 2, 0.00299473, 1.396, 0.862, 1.55, 2.989, 4.54, 6.11, 7.18, 7.44, 7.874, 8.9, 8.902, 8.96,
             7.133, 5.904, 5.323, 5.73, 4.5, 3.1028, 2.418, 1.532, 2.54, 4.469, 6.506, 8.57, 10.22, 11.5, 12.41, 12.41,
             12.02, 10.5, 8.65, 7.31, 7.31]

print("Number of Names:", len(Names))
print("Number of Densities:", len(Densities))

Cmap = cm.viridis

#VolumesStr = ''
#for i, Name in enumerate(Names):
#   VolumesStr += '        <volume name ="ShieldVol_' + str(i) + '">\n            <materialref ref="' + Name + '"/>\n            <solidref ref="Shield_' + str(i) + '"/>\n        </volume>\n\n'
#print(VolumesStr)

# for x in Densities:
#    print(x)

Path = "/home/anton/Desktop/triton_work/Permutations/2Layer50/Res/"

Electrons = totalkRadGras(Path, "Elec")

Protons = totalkRadGras(Path, "Prot")

print("Electrons Shape:", np.shape(Electrons))

Total = Electrons + Protons
Total[1] = np.sqrt(Electrons[1] * Electrons[1] + Protons[1] * Protons[1])
# print(Total)


NumMat = 50

for x in range(NumMat):
    for y in range(NumMat):
        i = x * NumMat + y
        print(i+1, x+1, y+1, Names[x], Names[y],  ufloat(Electrons[0][x], Electrons[1][x]), ufloat(Protons[0][x], Protons[1][x]), ufloat(Total[0][x], Total[1][x]))

'''
NumTiles = np.shape(Protons)[1]
print("NumTiles:", NumTiles)

ProtonMap = np.zeros((NumMat, NumMat), dtype=float)
ElectronMap = np.zeros((NumMat, NumMat), dtype=float)

for x in range(NumMat):
    for y in range(NumMat):
        ProtonMap[x][y] = Protons[0][x * NumMat + y]
        ElectronMap[x][y] = Electrons[0][x * NumMat + y]


fig1 = plt.figure(1)
image = plt.imshow(np.rot90(np.transpose(ElectronMap)), cmap=Cmap, extent=(0.5, len(Names)+0.5, 0.5, len(Names)+0.5))
ax = plt.gca()
ax.set_xticks(np.arange(5, 55, 5))
ax.set_yticks(np.arange(5, 55, 5))
ax.set_xticklabels(np.arange(5, 55, 5))
ax.set_yticklabels(np.arange(5, 55, 5))
cbar = plt.colorbar()
cbar.set_label("Ionizing dose per month [kRad]")
plt.title("Electron dose for two layer shielding")
plt.xlabel("Z-number of bottom-layer")
plt.ylabel("Z-number of top-layer")
plt.savefig(Path + "../ElectronMap.pdf", format='pdf', bbox_inches="tight")

fig2 = plt.figure(2)
plt.imshow(np.rot90(np.transpose(ProtonMap)), cmap=Cmap, extent=(0.5, len(Names)+0.5, 0.5, len(Names)+0.5))
ax = plt.gca()
ax.set_xticks(np.arange(5, 55, 5))
ax.set_yticks(np.arange(5, 55, 5))
ax.set_xticklabels(np.arange(5, 55, 5))
ax.set_yticklabels(np.arange(5, 55, 5))
cbar = plt.colorbar()
cbar.set_label("Ionizing dose per month [kRad]")
plt.title("Proton dose for two layer shielding")
plt.xlabel("Z-number of bottom-layer")
plt.ylabel("Z-number of top-layer")
plt.savefig(Path + "../ProtonMap.pdf", format='pdf', bbox_inches="tight")

TotalMap = ElectronMap + ProtonMap

fig3 = plt.figure(3)
plt.imshow(np.rot90(np.transpose(TotalMap)), cmap=Cmap, extent=(0.5, len(Names)+0.5, 0.5, len(Names)+0.5))
ax = plt.gca()
ax.set_xticks(np.arange(5, 55, 5))
ax.set_yticks(np.arange(5, 55, 5))
ax.set_xticklabels(np.arange(5, 55, 5))
ax.set_yticklabels(np.arange(5, 55, 5))
cbar = plt.colorbar()
cbar.set_label("Ionizing dose per month [kRad]")
plt.title("Total Ionizing dose for two layer shielding")
plt.xlabel("Z-number of bottom-layer")
plt.ylabel("Z-number of top-layer")
plt.savefig(Path + "../TotalMap.pdf", format='pdf', bbox_inches="tight")
'''