from GRAS.Dependencies.TotalKRadGras import totalkRadGras
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from uncertainties import ufloat
import sigfig
from GRAS.GPT4Experiments.CombinedColormaps import create_average_colormap

Names = ["G4_Li", "G4_Be", "G4_B", "G4_C", "G4_Mg", "G4_Al", "G4_Si", "G4_S", "G4_Ca", "G4_Sc", "G4_Ti", "G4_V",
         "G4_Cr", "G4_Mn", "G4_Fe", "G4_Co", "G4_Ni", "G4_Cu", "G4_Zn", "G4_Ge", "G4_As", "G4_Se", "G4_Sr", "G4_Y",
         "G4_Zr", "G4_Nb", "G4_Mo", "G4_Ru", "G4_Rh", "G4_Pd", "G4_Ag", "G4_Cd", "G4_In", "G4_Sn", "G4_Sb", "G4_Te",
         "G4_I", "G4_Ba", "G4_La", "G4_Ce", "G4_Pr", "G4_Nd", "G4_Sm", "G4_Eu", "G4_Gd", "G4_Tb", "G4_Dy", "G4_Ho",
         "G4_Er", "G4_Tm", "G4_Yb", "G4_Lu", "G4_Hf", "G4_Ta", "G4_W", "G4_Re", "G4_Os", "G4_Ir", "G4_Pt", "G4_Au",
         "G4_Tl", "G4_Pb"]

Densities = [0.534, 1.848, 2.37, 2, 1.74, 2.699, 2.33, 1.62, 2.989, 4.54, 6.11, 7.18, 7.44, 7.874, 8.9, 8.902, 8.96,
             7.133, 5.323, 5.73, 4.5, 2.54, 4.469, 6.506, 8.57, 10.22, 12.41, 12.02, 10.5, 8.65, 7.31, 6.691, 6.24,
             4.93, 3.5, 6.154, 6.657, 6.71, 6.9, 7.46, 5.243, 7.9004, 8.229, 8.62, 8.795, 9.066, 9.321, 6.73, 9.84,
             13.31, 16.654, 19.3, 21.02, 22.57, 22.42, 21.45, 19.32, 11.72, 11.35]

'''
VolumesStr = ''
for i, Name in enumerate(Names):
    VolumesStr += '        <volume name ="ShieldVol_' + str(i) + '">\n            <materialref ref="' + Name + '"/>\n            <solidref ref="Shield_' + str(i) + '"/>\n        </volume>\n\n'
print(VolumesStr)



file = "/l/triton_work/Permutations/NIST.txt"
Text = []

with open(file, 'r') as f:
    reader = csv.reader(f)
    for line in reader:
        if 'G4' in line[0] and len(line) == 1:
            Text.append(line[0].split())

Data = {}

for i, line in enumerate(Text):
    if 'G4' in line[1]: # and i < 97:
        Data[line[1]] = float(line[2])

for i, mat in enumerate(Data):
    print(i, mat, Data[mat])



Data["FR4"] = 1.85
Data["G4_POLYETHYLENE"] = 0.95
Data["G4_LITHIUM_HYDRIDE"] = 0.78

for x, top in enumerate(Data):
    for y, bottom in enumerate(Data):
        print(x*100+y, x, y, top, bottom)
        #print(Data[x])


MatrixStr = '<matrix name="MaterialDens" coldim="1" values="'

for key in Data:
    MatrixStr += str(Data[key]) + ' '

MatrixStr += '" />'

print(MatrixStr)

VolumesStr = ''

for i, key in enumerate(Data):
    VolumesStr += '        <volume name ="ShieldVol_' + str(i) + '">\n            <materialref ref="' + key + '"/>\n            <solidref ref="Shield_' + str(i) + '"/>\n        </volume>\n\n'

print(VolumesStr)
'''

Materials = Names

Path = "/l/triton_work/Permutations/2Layer62/Res/"
file_name = Path + "../Analysis/2Layer62-Raw.csv"

Electrons = totalkRadGras(Path, "Elec")
Protons = totalkRadGras(Path, "Prot")

print("Electrons Shape:", np.shape(Electrons))

Total = Electrons + Protons
Total[1] = np.sqrt(Electrons[1] * Electrons[1] + Protons[1] * Protons[1])
# print(Total)

NumMat = 62
List = []

#for i1 in range(NumMat):
#    for i2 in range(NumMat):
#            i = i1 * NumMat + i2
#            print(i+1, i1, i2, Materials[i1], Materials[i2], ufloat(Electrons[0][i], Electrons[1][i]), ufloat(Protons[0][i], Protons[1][i]), ufloat(Total[0][i], Total[1][i]))

with open(file_name, 'w') as file:
    file.write("Combination #,Mat 1 Z-Num,Mat 2 Z-Num,Material 1,Material 2,Electron Dose [krad/Month],Electron Err [krad/Month],Proton Dose [krad/Month],Proton Err [krad/Month],Total Dose [krad/Month],Total Err [krad/Month]\n")
    for i1 in range(NumMat):
        for i2 in range(NumMat):
            i = i1 * NumMat + i2
            line = f"{i+1},{i1+1},{i2+1},{Materials[i1]},{Materials[i2]},{ufloat(Electrons[0][i], Electrons[1][i])},{ufloat(Protons[0][i], Protons[1][i])},{ufloat(Total[0][i], Total[1][i])}\n"
            line = line.replace("+/-", ",")
            file.write(line)



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
ax.set_xticks(np.arange(5, 62, 5))
ax.set_yticks(np.arange(5, 62, 5))
ax.set_xticklabels(np.arange(5, 62, 5))
ax.set_yticklabels(np.arange(5, 62, 5))

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
ax.set_xticks(np.arange(5, 62, 5))
ax.set_yticks(np.arange(5, 62, 5))
ax.set_xticklabels(np.arange(5, 62, 5))
ax.set_yticklabels(np.arange(5, 62, 5))

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
ax.set_xticks(np.arange(5, 62, 5))
ax.set_yticks(np.arange(5, 62, 5))
ax.set_xticklabels(np.arange(5, 62, 5))
ax.set_yticklabels(np.arange(5, 62, 5))

cbar_ticks = np.geomspace(TotalMap.min()*1.05, TotalMap.max()*0.95, num=8)  # Adjust the number of ticks as needed
cbar_ticks = [sigfig.round(num, sigfigs=2) for num in cbar_ticks]
cbar = plt.colorbar(ticks=np.log(cbar_ticks))
cbar.set_ticklabels(cbar_ticks)

cbar.set_label("Ionizing dose per month [kRad]")
plt.title("Total Ionizing dose for two layer shielding")
plt.xlabel("Z-number of bottom-layer")
plt.ylabel("Z-number of top-layer")
plt.savefig(Path + "../Plots/TotalMap.pdf", format='pdf', bbox_inches="tight")
