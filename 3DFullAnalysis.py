import numpy as np
from ReadMultipleRoot import readMultipleRoot
from MeVtokRad_3D import MeVtokRad_3D

Path = "/home/anton/Desktop/triton_work/3D/SiChipTestModular/root/"
File = "1mmal3dmultichip2e9electrons500kev.root"
NORM_FACTOR_SPECTRUM = 5.886798E+14
Npart = 2e9
Radius = 10  # cm

GunEnergy, Edep, Esec = readMultipleRoot(Path + File)

GunEnergyTotalMeV = sum(GunEnergy)

NumberVolumes = len(Edep)

EdepTotalMeV = []
EsecTotalMeV = []

for i in range(NumberVolumes):
    EdepTotalMeV.append(sum(Edep[i]))
    EsecTotalMeV.append(sum(Esec[i]))

print("GunEnergyTotalMeV:", GunEnergyTotalMeV)
print("EdepTotalMeV:", EdepTotalMeV[0])
print("EsecTotalMeV:", EsecTotalMeV[0])

NumPoints = len(GunEnergy)
SampleLen = int(NumPoints / 100)

GunEnergySamples = np.zeros(100)
EdepSamples = np.zeros((NumberVolumes, 100))
EsecSamples = np.zeros((NumberVolumes, 100))

for i in range(100):
    GunEnergySamples[i] = sum(GunEnergy[SampleLen * i:SampleLen * (i + 1) - 1])

    for j in range(NumberVolumes):
        EdepSamples[j][i] = sum(Edep[j][SampleLen * i:SampleLen * (i + 1) - 1])
        EsecSamples[j][i] = sum(Esec[j][SampleLen * i:SampleLen * (i + 1) - 1])
    # print("From:", SampleLen*i)
    # print("To:", SampleLen*(i+1)-1)

GunEnergyStdMeV = np.std(GunEnergySamples)
EdepStdMeV = []
EsecStdMeV = []

for j in range(NumberVolumes):
    EdepStdMeV.append(np.std(EdepSamples[j][:]))
    EsecStdMeV.append(np.std(EsecSamples[j][:]))

print("GunEnergyStdMeV %:", 100 * GunEnergyStdMeV / GunEnergyTotalMeV)
print("EdepStdMeV %:", 100 * EdepStdMeV[0] / EdepTotalMeV[0])
print("EsecStdMeV %:", 100 * EsecStdMeV[0] / EsecTotalMeV[0])

GunEnergyTotalKRAD = MeVtokRad_3D(GunEnergyTotalMeV, NORM_FACTOR_SPECTRUM, Npart, Radius)
GunEnergyStdKRAD = MeVtokRad_3D(GunEnergyStdMeV, NORM_FACTOR_SPECTRUM, Npart, Radius)

EdepTotalKRAD = []
EsecTotalKRAD = []
EdepStdKRAD = []
EsecStdKRAD = []

for i in range(NumberVolumes):
    EdepTotalKRAD.append(MeVtokRad_3D(EdepTotalMeV[i], NORM_FACTOR_SPECTRUM, Npart, Radius))
    EsecTotalKRAD.append(MeVtokRad_3D(EsecTotalMeV[i], NORM_FACTOR_SPECTRUM, Npart, Radius))
    EdepStdKRAD.append(MeVtokRad_3D(EdepStdMeV[i], NORM_FACTOR_SPECTRUM, Npart, Radius))
    EsecStdKRAD.append(MeVtokRad_3D(EsecStdMeV[i], NORM_FACTOR_SPECTRUM, Npart, Radius))


print("GunEnergyTotalKRAD:", GunEnergyTotalKRAD)
for i in range(NumberVolumes):
    print("EdepTotalKRAD", str(i), EdepTotalKRAD[i], "+-", EdepStdKRAD[i])
