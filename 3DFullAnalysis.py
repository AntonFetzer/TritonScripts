import numpy as np
from ReadMultipleRoot import readMultipleRoot
from MeVtokRad_3D import MeVtokRad_3D

Path = "/home/anton/Desktop/triton_work/3D/SiChipTestModular/root/"
File = "1mmal3dmultichip1e8solarproton10mev.root"

Data = readMultipleRoot(Path + File)

NumberOfDataSets = len(Data)

TotalMeV = []

for i in range(NumberOfDataSets):
    TotalMeV.append(sum(Data[i]))

print(TotalMeV)
