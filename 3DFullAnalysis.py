import numpy as np
from ReadMultipleRoot import readMultipleRoot
from MeVtokRad_3D import MeVtokRad_3D

Path = "/home/anton/Desktop/triton_work/3D/SiChipTestModular/root/"
File = "1mmal3dmultichip1e8solarproton10mev.root"

GunEnergy, Edep, Esec = readMultipleRoot(Path + File)

GunEnergyMeV = sum(GunEnergy)

NumberVolumes = len(Edep)

EdepTotalMeV = []
EsecTotalMeV = []

for i in range(NumberVolumes):
    EdepTotalMeV.append(sum(Edep[i]))
    EsecTotalMeV.append(sum(Esec[i]))

print("GunEnergyMeV:", GunEnergyMeV)
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

print(GunEnergy / EdepTotalMeV[0])
