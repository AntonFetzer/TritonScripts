import csv
from TotalKRadGras import totalkRadGras
import numpy as np
import matplotlib.pyplot as plt

file = "/home/anton/Desktop/triton_work/Permutations/NIST.txt"
Text = []

with open(file, 'r') as f:
    reader = csv.reader(f)
    for line in reader:
        if 'G4' in line[0] and len(line) == 1:
            Text.append(line[0].split())

Data = {}

for i, line in enumerate(Text):
    if 'G4' in line[1] and float(line[2]) > 0.5:
        Data[line[1]] = float(line[2])

'''
MatrixStr = '<matrix name="MaterialThicknesses" coldim="1" values="'

for key in Data:
    MatrixStr += str(Data[key]) + ' '

MatrixStr += '" />'

print(MatrixStr)

VolumesStr = ''

for i, key in enumerate(Data):
    VolumesStr += '        <volume name ="ShieldVol_' + str(i) + '">\n            <materialref ref="' + key + '"/>\n            <solidref ref="Shield_' + str(i) + '"/>\n        </volume>\n\n'

print(VolumesStr)
'''

Path = "/home/anton/Desktop/triton_work/Permutations/1Layer/Res/"

Electrons = totalkRadGras(Path, "Elec") * 299

Protons = totalkRadGras(Path, "Prot") * 299

print(np.shape(Electrons))

x = np.linspace(1, 300, num=299, dtype=int)
plt.errorbar(x, Electrons[0], Electrons[1], fmt=' ', capsize=2)
plt.errorbar(x, Protons[0], Protons[1], fmt=' ', capsize=2)
plt.show()

Total = Electrons + Protons

Total[1] = np.sqrt(Electrons[1] * Electrons[1] + Protons[1] * Protons[1])

print(Total)

ResultsStr = ''

for i, key in enumerate(Data):
    ResultsStr += key + " " + str(Data[key]) + " " + str(Electrons[0][i]) + " " + str(Electrons[1][i]) + " " + str(Protons[0][i]) + " " + str(Protons[1][i]) + " " + str(Total[0][i]) + " " + str(Total[1][i]) + "\n"

CSVFile = open(Path + "../Results.txt", 'w')
CSVFile.writelines("Material Density Electrons ElecErr Protons ProtErr TotalDose TotalErr\n")
CSVFile.writelines(ResultsStr + "\n")
CSVFile.close()
