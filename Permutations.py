import numpy as np
from Dependencies.MeVtokRad_2D import MeVtokRad_2D

Path = "/home/anton/Desktop/triton_work/Permutations/5Layer1-5gcm2/Materials3/"

ElecFile = "5layer3-2e9electron.txt"
ProtFile = "5layer3-2e7proton.txt"

NumberOfLayers = 5
Materials = ["PE", "Al", "FR4", "Sn"]
NumMaterials = len(Materials)

Npart_Elec = 2e9 / NumMaterials**NumberOfLayers
Npart_Prot = 2e7 / NumMaterials**NumberOfLayers
NORM_FACTOR_SPECTRUM_Elec = 5.886798E+14  # Elec500keV
NORM_FACTOR_SPECTRUM_Prot = 3.381390E+11  # Prots10MeV

Thicknesses = [3.19148936170213, 1.11152278621712, 1.62162162162162, 0.410396716826265]

Prot = np.genfromtxt(Path + ProtFile, delimiter=',', dtype=None, encoding='ASCII')
Elec = np.genfromtxt(Path + ElecFile, delimiter=',', dtype=None, encoding='ASCII')

NumPoints = len(Elec)
ElecEdepMeV = np.zeros(NumPoints)
ElecErrMeV = np.zeros(NumPoints)
ElecEdepRad = np.zeros(NumPoints)
ElecErrRad = np.zeros(NumPoints)
ProtEdepMeV = np.zeros(NumPoints)
ProtErrMeV = np.zeros(NumPoints)
ProtEdepRad = np.zeros(NumPoints)
ProtErrRad = np.zeros(NumPoints)

for i in range(NumPoints):
    ElecEdepMeV[i] = Elec[i][1]
    ElecErrMeV[i] = Elec[i][2]
    ProtEdepMeV[i] = Prot[i][1]
    ProtErrMeV[i] = Prot[i][2]

ElecEdepRad = MeVtokRad_2D(ElecEdepMeV, NORM_FACTOR_SPECTRUM_Elec, Npart_Elec)
ElecErrRad = MeVtokRad_2D(ElecErrMeV, NORM_FACTOR_SPECTRUM_Elec, Npart_Elec)
ProtEdepRad = MeVtokRad_2D(ProtEdepMeV, NORM_FACTOR_SPECTRUM_Prot, Npart_Prot)
ProtErrRad = MeVtokRad_2D(ProtErrMeV, NORM_FACTOR_SPECTRUM_Prot, Npart_Prot)

TotalEdep = ElecEdepRad + ProtEdepRad
TotalErr = np.sqrt(ElecErrRad*ElecErrRad + ProtErrRad*ProtErrRad)

CSVFile = open(Path + "Results.txt", 'w')
CSVFile.writelines("Number, Material A, Material B, Material C, Material D, Material E, Electron Dose, Electron Error, Proton Dose, Proton Error, Total Dose, Total Error, Total Thick \n \n")

for i5 in range(NumMaterials):
    for i4 in range(NumMaterials):
        for i3 in range(NumMaterials):
            for i2 in range(NumMaterials):
                for i1 in range(NumMaterials):
                    # print(i3, i2, i1)
                    i = i5*NumMaterials**4 + i4*NumMaterials**3 + i3*NumMaterials**2 + i2*NumMaterials**1 + i1

                    TotalThick = Thicknesses[i1] + Thicknesses[i2] + Thicknesses[i3] + Thicknesses[i4] + Thicknesses[i5]

                    List = (i, Materials[i1], Materials[i2], Materials[i3], Materials[i4], Materials[i5], ElecEdepRad[i], ElecErrRad[i], ProtEdepRad[i],
                            ProtErrRad[i], TotalEdep[i], TotalErr[i], TotalThick)
                    String = ', '.join(map(str, List))
                    print(String)
                    CSVFile.writelines(String + "\n")

CSVFile.close()
