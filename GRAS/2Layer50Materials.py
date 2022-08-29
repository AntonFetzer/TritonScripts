import csv
import nmmn.plots
from TotalKRadGras import totalkRadGras
import numpy as np
import matplotlib.pyplot as plt

Names = ["G4_lH2", "G4_He", "G4_Li", "G4_Be", "G4_B", "G4_C", "G4_lN2", "G4_lO2", "G4_F", "G4_Ne", "G4_Na", "G4_Mg",
         "G4_Al", "G4_Si", "G4_P", "G4_S", "G4_Cl", "G4_lAr", "G4_K", "G4_Ca", "G4_Sc", "G4_Ti", "G4_V", "G4_Cr",
         "G4_Mn", "G4_Fe", "G4_Co", "G4_Ni", "G4_Cu", "G4_Zn", "G4_Ga", "G4_Ge", "G4_As", "G4_Se", "G4_lBr", "G4_lKr",
         "G4_Rb", "G4_Sr", "G4_Y", "G4_Zr", "G4_Nb", "G4_Mo", "G4_Tc", "G4_Ru", "G4_Rh", "G4_Ta", "G4_W", "G4_Os",
         "G4_Pb", "G4_U"]

Densities = [0.07, 0.00, 0.53, 1.85, 2.37, 2.00, 0.81, 1.14, 0.97, 1.74, 2.70, 2.33, 2.20, 1.40, 0.86, 1.55, 2.99, 4.54,
             6.11, 7.18, 7.44, 7.87, 8.90, 8.96, 7.13, 5.90, 5.32, 5.73, 4.50, 3.10, 2.42, 1.53, 2.54, 4.47, 6.51, 8.57,
             10.22, 11.50, 12.41, 12.02, 16.65, 19.30, 22.57, 11.35, 18.95]

print("Number of Names:", len(Names))
print("Number of Densities:", len(Densities))

'''
VolumesStr = ''
for i, Name in enumerate(Names):
    VolumesStr += '        <volume name ="ShieldVol_' + str(i) + '">\n            <materialref ref="' + Name + '"/>\n            <solidref ref="Shield_' + str(i) + '"/>\n        </volume>\n\n'
print(VolumesStr)



file = "/home/anton/Desktop/triton_work/Permutations/NIST.txt"
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


Path = "/home/anton/Desktop/triton_work/Permutations/2Layer50/Res/"

Electrons = totalkRadGras(Path, "Elec")

Protons = totalkRadGras(Path, "Prot")

print("Electrons Shape:", np.shape(Electrons))

Total = Electrons + Protons
Total[1] = np.sqrt(Electrons[1] * Electrons[1] + Protons[1] * Protons[1])
# print(Total)


NumMat = 62

for x in range(62):
    for y in range(62):
        i = x * NumMat + y
        print(i, x, y, Names[x], Names[y], Electrons[0][i], Electrons[1][i], Protons[0][i], Protons[1][i])

NumTiles = np.shape(Protons)[1]
print("NumTiles:", NumTiles)

ProtonMap = np.zeros((NumMat, NumMat), dtype=float)
ElectronMap = np.zeros((NumMat, NumMat), dtype=float)

for x in range(NumMat):
    for y in range(NumMat):
        ProtonMap[x][y] = Protons[0][x * NumMat + y]
        ElectronMap[x][y] = Electrons[0][x * NumMat + y]

turbo = nmmn.plots.turbocmap()
fig1 = plt.figure(1)
plt.imshow(ElectronMap, cmap=turbo)
cbar = plt.colorbar()
cbar.set_label("Electron Dose in krad")
# plt.show()
plt.savefig(Path + "../ElectronMap.eps", format='eps', bbox_inches="tight")

fig2 = plt.figure(2)
plt.imshow(ProtonMap, cmap=turbo)
cbar = plt.colorbar()
cbar.set_label("Proton Dose in krad")
# plt.show()
plt.savefig(Path + "../ProtonMap.eps", format='eps', bbox_inches="tight")

TotalMap = ElectronMap + ProtonMap

fig3 = plt.figure(3)
plt.imshow(TotalMap, cmap=turbo)
cbar = plt.colorbar()
cbar.set_label("Total Dose in krad")
# plt.show()
plt.savefig(Path + "../TotalMap.eps", format='eps', bbox_inches="tight")
'''
