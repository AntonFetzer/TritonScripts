import numpy as np
from ReadMultipleRoot import readMultipleRoot
from MeVtokRad_3D import MeVtokRad_3D
from natsort import natsorted
import os

Path = "/home/anton/Desktop/triton_work/3D/MultiChipTest/root/CosmicProton/"
# NORM_FACTOR_SPECTRUM = 5.886798E+14  # Electron500keV
# NORM_FACTOR_SPECTRUM = 3.381390E+11  # Protons10MeV
# NORM_FACTOR_SPECTRUM = 1.109681E+10  # SolarProton10MeV
NORM_FACTOR_SPECTRUM = 2.024537E+07  # CosmicProtonsFull

Npart = 2e9
Radius = 10  # cm

# Get list of all root files in that folder
Files = [f for f in os.listdir(Path) if f.endswith('.root')]
Files = natsorted(Files)

NumDataSets = 0
TotalMeV = []
StdMeV = []
TotalKRAD = []
StdKRAD = []

CSVfile = open(Path + Files[0].split(".")[0].split("-")[1] + "DoseTable.csv", 'w')

for File in Files:

    Data = readMultipleRoot(
        Path + File)  # ['Sivol_0_Edep_MeV', ... , 'Sivol_0_Esec_MeV', ... , 'Gun_energy_MeV', 'Gun_angle_deg']

    NumDataSets = len(Data)
    NumSivols = int(NumDataSets / 2) - 1

    TotalMeV = []

    for i in range(NumDataSets):
        TotalMeV.append(sum(Data[i]))

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

    # --------- Conversion ##############
    TotalKRAD = []
    StdKRAD = []

    for i in range(NumDataSets):
        TotalKRAD.append(MeVtokRad_3D(TotalMeV[i], NORM_FACTOR_SPECTRUM, Npart, Radius))
        StdKRAD.append(MeVtokRad_3D(StdMeV[i], NORM_FACTOR_SPECTRUM, Npart, Radius))

    print("GunEnergyTotalKRAD:", TotalKRAD[-2], "+-", 100 * StdKRAD[-2] / TotalKRAD[-2])
    for i in range(NumSivols):
        print("EdepTotalKRAD", str(i), TotalKRAD[i], "+-", 100 * StdKRAD[i] / TotalKRAD[i])
        print("EsecTotalKRAD", str(i), TotalKRAD[i + 4], "+-", 100 * StdKRAD[i + 4] / TotalKRAD[i + 4])

    # ---------- Write to CSV File ----------------------------------

    CSVfile.write(File + ", Dose Table in kRad, Number of Points: " + str(NumPoints) + "\n")
    CSVfile.write(File.split("-")[0] + "mm Al Shielding Thickness, Edep, EdepStd, Esec, EsecStd\n")
    for i in range(NumSivols):
        CSVfile.write(','.join(["Sivol_" + str(i), str(TotalKRAD[i]), str(StdKRAD[i]), str(TotalKRAD[i + NumSivols]), str(StdKRAD[i + NumSivols]) + "\n"]))

    CSVfile.write("\n")
CSVfile.close()
