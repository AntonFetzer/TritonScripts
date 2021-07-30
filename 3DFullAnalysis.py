import numpy as np
from ReadMultipleRoot import readMultipleRoot
from MeVtokRad_3D import MeVtokRad_3D
import matplotlib.pyplot as plt

Path = "/home/anton/Desktop/triton_work/3D/MultiChipTest/root/Electron/"
File = "1-mmal3dmultichip2e9electrons500kev.root"
NORM_FACTOR_SPECTRUM = 5.886798E+14
Npart = 2e9
Radius = 10  # cm

Data = readMultipleRoot(
    Path + File)  # ['Sivol_0_Edep_MeV', ... , 'Sivol_0_Esec_MeV', ... , 'Gun_energy_MeV', 'Gun_angle_deg']

NumDataSets = len(Data)

TotalMeV = []

for i in range(NumDataSets):
    TotalMeV.append(sum(Data[i]))

print("Edep1TotalMeV:", TotalMeV[0])

# -------------------------- Standard Deviation ---------------------------
NumPoints = len(Data[0])
SampleLen = int(NumPoints / 100)

Samples = np.zeros(100)
StdMeV = []

for j in range(NumDataSets):
    for i in range(100):
        Samples[i] = sum(Data[j][SampleLen * i:SampleLen * (i + 1) - 1])
    # print("From:", SampleLen*i)
    # print("To:", SampleLen*(i+1)-1)
    # print("Samples:", Samples)
    StdMeV.append(np.std(Samples))

print("Edep1StdMeV:", StdMeV[0])

# --------- Conversion ##############
TotalKRAD = []
StdKRAD = []

for i in range(NumDataSets):
    TotalKRAD.append(MeVtokRad_3D(TotalMeV[i], NORM_FACTOR_SPECTRUM, Npart, Radius))
    StdKRAD.append(MeVtokRad_3D(StdMeV[i], NORM_FACTOR_SPECTRUM, Npart, Radius))

print("Edep1TotalKRAD:", TotalKRAD[0])
print("Edep1StdKRAD:", StdKRAD[0])

print("Edep1StdMeV %:", 100 * TotalKRAD[0] / StdKRAD[0])
print("Esec1StdMeV %:", 100 * TotalKRAD[4] / StdKRAD[4])

print("GunEnergyTotalKRAD:", TotalKRAD[-2])
for i in range(4):
    print("EdepTotalKRAD", str(i), TotalKRAD[i], "+-", TotalKRAD[i])
    print("EsecTotalKRAD", str(i), TotalKRAD[i + 4], "+-", TotalKRAD[i + 4])
