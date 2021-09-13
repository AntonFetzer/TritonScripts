import os

import numpy as np
from natsort import natsorted

from Dependencies.EstimateError import EstimateError
from Dependencies.MeVtokRad_3D import MeVtokRad_3D
from Dependencies.ReadMultipleRoot import readMultipleRoot

Path = "/home/anton/Desktop/triton_work/3D/AlFrame/AlFramePE-W/root/"
SiChipLen = 0.8  # cm --------------------------------------------------------------------------------------------------

#                       Electron500keV, Protons10MeV, ElectronsFull, ProtonFUll, solarproton, cosmicproton
NORM_FACTOR_SPECTRUM_List = [5.886798E+14, 3.381390E+11, 6.159454E+15, 8.003046E+14, 1.109681E+10, 2.024537E+07]
SpectrumTypes = ["electrons500kev.root", "protons10mev.root", "electronsfull.root", "protonsfull.root", "solarproton10mev.root", "cosmicproton.root"]
NpartList = [2e9, 2e9, 2e9, 2e9, 2e9, 1e8]

for Type in range(6):

    # Get list of all root files in that folder
    Files = [f for f in os.listdir(Path) if f.endswith(SpectrumTypes[Type])]
    if not Files:
        continue
    Files = natsorted(Files)

    Spectrum = SpectrumTypes[Type]
    NORM_FACTOR_SPECTRUM = NORM_FACTOR_SPECTRUM_List[Type]
    Npart = NpartList[Type]
    Radius = 10  # cm

    TotalMeV = []
    StdMeV = []
    TotalKRAD = []
    StdKRAD = []
    lol = 0

    CSVfile = open(Path + Files[0].split(".")[0] + "DoseTable.csv", 'w')

    for File in Files:

        Data = readMultipleRoot(Path + File)  # ['Sivol_0_Edep_MeV', ... , 'Sivol_0_Esec_MeV', ... , 'Gun_energy_MeV', 'Gun_angle_deg']

        NumDataSets = len(Data)
        NumSivols = int(NumDataSets / 2) - 1

        TotalMeV = []

        for i in range(NumDataSets):
            TotalMeV.append(sum(Data[i]))

        # -------------------------- Standard Deviation ---------------------------
        NumPoints = len(Data[0])
        tesz = 1
        Samples = np.zeros(100)
        StdMeV = []

        for j in range(NumDataSets):
            StdMeV.append(EstimateError(Data[j], 100))

        # --------- Conversion ##############
        TotalKRAD = []
        StdKRAD = []

        for i in range(NumDataSets):
            TotalKRAD.append(MeVtokRad_3D(TotalMeV[i], NORM_FACTOR_SPECTRUM, Npart, Radius, SiChipLen))
            StdKRAD.append(MeVtokRad_3D(StdMeV[i], NORM_FACTOR_SPECTRUM, Npart, Radius, SiChipLen))

        print("GunEnergyTotalKRAD:", TotalKRAD[-2], "+-", 100 * StdKRAD[-2] / TotalKRAD[-2], "%")
        for i in range(NumSivols):
            print("EdepTotalKRAD", str(i), TotalKRAD[i], "+-", 100 * StdKRAD[i] / TotalKRAD[i], "%")
            print("EsecTotalKRAD", str(i), TotalKRAD[i + 4], "+-", 100 * StdKRAD[i + 4] / TotalKRAD[i + 4], "%")

        # ---------- Write to CSV File ----------------------------------

        CSVfile.write(File + ", Dose Table in kRad, Number of Points: " + str(NumPoints) + "\n")
        CSVfile.write(File.split("-")[0] + "mm Al Shielding Thickness, Edep, EdepStd, Esec, EsecStd\n")
        for i in range(NumSivols):
            CSVfile.write(','.join(["Sivol_" + str(i), str(TotalKRAD[i]), str(StdKRAD[i]), str(TotalKRAD[i + NumSivols]), str(StdKRAD[i + NumSivols]) + "\n"]))

        CSVfile.write("\n")
    CSVfile.close()
