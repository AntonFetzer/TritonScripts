import csv

file = "/home/anton/Desktop/triton_work/Materials/NIST.txt"
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

MatrixStr = '<matrix name="MaterialThicknesses" coldim="1" values="'

for key in Data:
    MatrixStr += str(Data[key]) + ' '

MatrixStr += '" />'

print(MatrixStr)

VolumesStr = ''

for i, key in enumerate(Data):
    VolumesStr += '        <volume name ="ShieldVol_' + str(i) + '">\n            <materialref ref="' + key + '"/>\n            <solidref ref="Shield_' + str(i) + '"/>\n        </volume>\n\n'

print(VolumesStr)

