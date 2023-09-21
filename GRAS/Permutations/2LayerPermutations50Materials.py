from GRAS.Dependencies.TotalKRadGras import totalkRadGras
import numpy as np
from uncertainties import ufloat

Materials = ["G4_lH2", "G4_He", "G4_Li", "G4_Be", "G4_B", "G4_C", "G4_lN2", "G4_lO2", "G4_F", "G4_Ne", "G4_Na", "G4_Mg",
         "G4_Al", "G4_Si", "G4_P", "G4_S", "G4_Cl", "G4_lAr", "G4_K", "G4_Ca", "G4_Sc", "G4_Ti", "G4_V", "G4_Cr",
         "G4_Mn", "G4_Fe", "G4_Co", "G4_Ni", "G4_Cu", "G4_Zn", "G4_Ga", "G4_Ge", "G4_As", "G4_Se", "G4_lBr", "G4_lKr",
         "G4_Rb", "G4_Sr", "G4_Y", "G4_Zr", "G4_Nb", "G4_Mo", "G4_Tc", "G4_Ru", "G4_Rh", "G4_Pd", "G4_Ag", "G4_Cd",
         "G4_In", "G4_Sn"]

Densities = [0.0708, 0.000166322, 0.534, 1.848, 2.37, 2, 0.807, 1.141, 0.00158029, 0.000838505, 0.971, 1.74, 2.699,
             2.33, 2.2, 2, 0.00299473, 1.396, 0.862, 1.55, 2.989, 4.54, 6.11, 7.18, 7.44, 7.874, 8.9, 8.902, 8.96,
             7.133, 5.904, 5.323, 5.73, 4.5, 3.1028, 2.418, 1.532, 2.54, 4.469, 6.506, 8.57, 10.22, 11.5, 12.41, 12.41,
             12.02, 10.5, 8.65, 7.31, 7.31]

print("Number of Names:", len(Materials))
print("Number of Densities:", len(Densities))



#VolumesStr = ''
#for i, Name in enumerate(Names):
#   VolumesStr += '        <volume name ="ShieldVol_' + str(i) + '">\n            <materialref ref="' + Name + '"/>\n            <solidref ref="Shield_' + str(i) + '"/>\n        </volume>\n\n'
#print(VolumesStr)

# for x in Densities:
#    print(x)

Names = Materials

Path = "/l/triton_work/Permutations/2Layer50/Res/"
file_name = Path + "../Analysis/2Layer50-Raw.csv"

Electrons = totalkRadGras(Path, "Elec")
Protons = totalkRadGras(Path, "Prot")

print("Electrons Shape:", np.shape(Electrons))

Total = Electrons + Protons
Total[1] = np.sqrt(Electrons[1] * Electrons[1] + Protons[1] * Protons[1])
# print(Total)

NumMat = 50
List = []

#for i1 in range(NumMat):
#    for i2 in range(NumMat):
#            i = i1 * NumMat + i2
#            print(i+1, i1, i2, Materials[i1], Materials[i2], ufloat(Electrons[0][i], Electrons[1][i]), ufloat(Protons[0][i], Protons[1][i]), ufloat(Total[0][i], Total[1][i]))

with open(file_name, 'w') as file:
    file.write("Combination #,Material 1 Z-Number,Material 2 Z-Number,Material 1,Material 2,Electron Dose [krad/Month],Electron Err [krad/Month],Proton Dose [krad/Month],Proton Err [krad/Month],Total Dose [krad/Month],Total Err [krad/Month]\n")
    for i1 in range(NumMat):
        for i2 in range(NumMat):
            i = i1 * NumMat + i2
            line = f"{i+1},{i1+1},{i2+1},{Materials[i1]},{Materials[i2]},{ufloat(Electrons[0][i], Electrons[1][i])},{ufloat(Protons[0][i], Protons[1][i])},{ufloat(Total[0][i], Total[1][i])}\n"
            line = line.replace("+/-", ",")
            file.write(line)


'''
NumTiles = np.shape(Protons)[1]
print("NumTiles:", NumTiles)

ProtonMap = np.zeros((NumMat, NumMat), dtype=float)
ElectronMap = np.zeros((NumMat, NumMat), dtype=float)
TotalMap = np.zeros((NumMat, NumMat), dtype=float)

for x in range(NumMat):
    for y in range(NumMat):
        ProtonMap[x][y] = Protons[0][x * NumMat + y]
        ElectronMap[x][y] = Electrons[0][x * NumMat + y]
        TotalMap[x][y] = Total[0][x * NumMat + y]


cmapE = cm.viridis
cmapP = cm.plasma
#cmapE = cm.jet
#cmapP = cm.turbo
cmapT = create_average_colormap(cmapE, cmapP)


fig1 = plt.figure(1)
plt.imshow(np.rot90(np.transpose(np.log(ElectronMap))), cmap=cmapE, extent=(0.5, len(Names)+0.5, 0.5, len(Names)+0.5))
ax = plt.gca()
ax.set_xticks(np.arange(5, 55, 5))
ax.set_yticks(np.arange(5, 55, 5))
ax.set_xticklabels(np.arange(5, 55, 5))
ax.set_yticklabels(np.arange(5, 55, 5))

cbar_ticks = np.geomspace(ElectronMap.min()*1.05, ElectronMap.max()*0.95, num=8)  # Adjust the number of ticks as needed
cbar_ticks = [sigfig.round(num, sigfigs=2) for num in cbar_ticks]
cbar = plt.colorbar(ticks=np.log(cbar_ticks))
cbar.set_ticklabels(cbar_ticks)

cbar.set_label("Ionizing dose per month [kRad]")
plt.title("Electron dose for two layer shielding")
plt.xlabel("Z-number of bottom-layer")
plt.ylabel("Z-number of top-layer")
plt.savefig(Path + "../Plots/ElectronMap.pdf", format='pdf', bbox_inches="tight")


fig2 = plt.figure(2)
plt.imshow(np.rot90(np.transpose(np.log(ProtonMap))), cmap=cmapP, extent=(0.5, len(Names)+0.5, 0.5, len(Names)+0.5))
ax = plt.gca()
ax.set_xticks(np.arange(5, 55, 5))
ax.set_yticks(np.arange(5, 55, 5))
ax.set_xticklabels(np.arange(5, 55, 5))
ax.set_yticklabels(np.arange(5, 55, 5))

cbar_ticks = np.geomspace(ProtonMap.min()*1.05, ProtonMap.max()*0.95, num=8)  # Adjust the number of ticks as needed
cbar_ticks = [sigfig.round(num, sigfigs=2) for num in cbar_ticks]
cbar = plt.colorbar(ticks=np.log(cbar_ticks))
cbar.set_ticklabels(cbar_ticks)

cbar.set_label("Ionizing dose per month [kRad]")
plt.title("Proton dose for two layer shielding")
plt.xlabel("Z-number of bottom-layer")
plt.ylabel("Z-number of top-layer")
plt.savefig(Path + "../Plots/ProtonMap.pdf", format='pdf', bbox_inches="tight")




fig3 = plt.figure(3)
plt.imshow(np.rot90(np.transpose(np.log(TotalMap))), cmap=cmapT, extent=(0.5, len(Names)+0.5, 0.5, len(Names)+0.5))
ax = plt.gca()
ax.set_xticks(np.arange(5, 55, 5))
ax.set_yticks(np.arange(5, 55, 5))
ax.set_xticklabels(np.arange(5, 55, 5))
ax.set_yticklabels(np.arange(5, 55, 5))

cbar_ticks = np.geomspace(TotalMap.min()*1.05, TotalMap.max()*0.95, num=8)  # Adjust the number of ticks as needed
cbar_ticks = [sigfig.round(num, sigfigs=2) for num in cbar_ticks]
cbar = plt.colorbar(ticks=np.log(cbar_ticks))
cbar.set_ticklabels(cbar_ticks)

cbar.set_label("Ionizing dose per month [kRad]")
plt.title("Total Ionizing dose for two layer shielding")
plt.xlabel("Z-number of bottom-layer")
plt.ylabel("Z-number of top-layer")
plt.savefig(Path + "../Plots/TotalMap.pdf", format='pdf', bbox_inches="tight")
'''
