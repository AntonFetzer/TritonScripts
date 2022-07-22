import numpy as np
import matplotlib.pyplot as plt
import csv

file = "/home/anton/Desktop/triton_work/1Tile/PhysicsComparison/Results/FTFP_BERT-Prot_757679_117218.csv"

DoseHist = []
PrimaryHist = []

ReadFlag = 0
print("Reading in File: " + file)
with open(file, 'r') as f:
    reader = csv.reader(f)
    for line in reader:
        if ReadFlag == 0:
            if "'Bin entries'" in line:
                ReadFlag = 1
        elif ReadFlag == 1:
            if "'End of Block'" in line:
                ReadFlag = 2
            else:
                DoseHist.append([float(x) for x in line])
        elif ReadFlag == 2:
            if "'Bin entries'" in line:
                ReadFlag = 3
        elif ReadFlag == 3:
            if "'End of Block'" in line:
                ReadFlag = 4
            else:
                PrimaryHist.append([float(x) for x in line])

DoseHist = np.asarray(DoseHist)

lowerID = 0
upperID = 1
meanID = 2
valueID = 3
errorID = 4
entriesID = 5

Bins = int(len(DoseHist)/10)
NewDoseHist = np.zeros((Bins, 6), dtype=float)

for x in range(len(NewDoseHist)):
    NewDoseHist[x][lowerID] = DoseHist[x*10][lowerID]
    NewDoseHist[x][upperID] = DoseHist[x*10+9][upperID]
    for i in range(10):
        NewDoseHist[x][valueID] += DoseHist[x*10 + i][valueID]
        NewDoseHist[x][errorID] += DoseHist[x*10 + i][errorID]

# plt.pyplot.bar(x, height, width=0.8, bottom=None, *, align='center', DoseHist=None, **kwargs)

plt.figure(1)
plt.bar(NewDoseHist[:, lowerID], NewDoseHist[:, valueID], width=NewDoseHist[:, upperID] - NewDoseHist[:, lowerID], align='edge')
plt.yscale("log")
plt.xscale("log")

PrimaryHist = np.asarray(PrimaryHist)

Bins = int(len(PrimaryHist)/10)
NewPrimaryHist = np.zeros((Bins, 6), dtype=float)

for x in range(len(NewPrimaryHist)):
    NewPrimaryHist[x][lowerID] = PrimaryHist[x*10][lowerID]
    NewPrimaryHist[x][upperID] = PrimaryHist[x*10+9][upperID]
    for i in range(10):
        NewPrimaryHist[x][valueID] += PrimaryHist[x*10 + i][valueID]
        NewPrimaryHist[x][errorID] += PrimaryHist[x*10 + i][errorID]

plt.figure(2)
plt.bar(NewPrimaryHist[:, lowerID], NewPrimaryHist[:, valueID], width=NewPrimaryHist[:, upperID] - NewPrimaryHist[:, lowerID], align='edge')
plt.yscale("log")
plt.xscale("log")


plt.show()
