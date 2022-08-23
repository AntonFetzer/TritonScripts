import csv
import nmmn.plots
from TotalKRadGras import totalkRadGras
import numpy as np
import matplotlib.pyplot as plt

Names = ["G4_Li", "G4_Be", "G4_B", "G4_C", "G4_Mg", "G4_Al", "G4_Si", "G4_S", "G4_Ca", "G4_Sc", "G4_Ti", "G4_V",
         "G4_Cr", "G4_Mn", "G4_Fe", "G4_Co", "G4_Ni", "G4_Cu", "G4_Zn", "G4_Ge", "G4_As", "G4_Se", "G4_Sr", "G4_Y",
         "G4_Zr", "G4_Nb", "G4_Mo", "G4_Ru", "G4_Rh", "G4_Pd", "G4_Ag", "G4_Cd", "G4_In", "G4_Sn", "G4_Sb", "G4_Te",
         "G4_I", "G4_Ba", "G4_La", "G4_Ce", "G4_Pr", "G4_Nd", "G4_Sm", "G4_Eu", "G4_Gd", "G4_Tb", "G4_Dy", "G4_Ho",
         "G4_Er", "G4_Tm", "G4_Yb", "G4_Lu", "G4_Hf", "G4_Ta", "G4_W", "G4_Re", "G4_Os", "G4_Ir", "G4_Pt", "G4_Au",
         "G4_Tl", "G4_Pb"]

Densities = [0.534, 1.848, 2.37, 2, 1.74, 2.699, 2.33, 1.55, 2.989, 4.54, 6.11, 7.18, 7.44, 7.874, 8.9, 8.902, 8.96,
             7.133, 5.323, 5.73, 4.5, 2.54, 4.469, 6.506, 8.57, 10.22, 12.41, 12.02, 10.5, 8.65, 7.31, 6.691, 6.24,
             4.93, 3.5, 6.154, 6.657, 6.71, 6.9, 7.46, 5.243, 7.9004, 8.229, 8.55, 8.795, 9.066, 9.321, 6.73, 9.84,
             13.31, 16.654, 19.3, 21.02, 22.57, 22.42, 21.45, 19.32, 11.72, 11.35]

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
'''

Path = "/home/anton/Desktop/triton_work/Permutations/2Layer/Res/"

Electrons = totalkRadGras(Path, "Elec")

Protons = totalkRadGras(Path, "Prot")

print("Electrons Shape:", np.shape(Electrons))

Total = Electrons + Protons
Total[1] = np.sqrt(Electrons[1] * Electrons[1] + Protons[1] * Protons[1])
# print(Total)


NumMat = 62

for x in range(62):
    for y in range(62):
        i = x*NumMat + y
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
#plt.imshow(ProtonMap, cmap=turbo)
plt.imshow(ElectronMap, cmap=turbo)

cbar = plt.colorbar()
cbar.set_label("Dose in krad")
# plt.yscale("log")
# plt.grid(which='both')
plt.show()