import os
from ReadGRASCSV1Tile import readGrasCsv1Tile
from MeVtokRad2DGras import MeVtokRad_2D
import matplotlib.pyplot as plt

path = "/home/anton/Desktop/triton_work/GRAS-1Mat/PhysicsComparison/Results/"
'''
ElecFiles = [f for f in os.listdir(path) if "Elec" in f]
ProtFiles = [f for f in os.listdir(path) if "Prot" in f]

Elec = dict.fromkeys(ElecFiles, {})
Prot = dict.fromkeys(ProtFiles, {})

NORM_FACTOR_SPECTRUM_Elec = 7.891281E+14
NORM_FACTOR_SPECTRUM_Prot = 3.389664E+11

for key in Elec:
    Elec[key] = readGrasCsv1Tile(path + key)
    Elec[key]["Dose"] = MeVtokRad_2D(Elec[key]["Dose"], NORM_FACTOR_SPECTRUM_Elec)
    Elec[key]["Error"] = MeVtokRad_2D(Elec[key]["Error"], NORM_FACTOR_SPECTRUM_Elec)

for key in Prot:
    Prot[key] = readGrasCsv1Tile(path + key)
    Prot[key]["Dose"] = MeVtokRad_2D(Prot[key]["Dose"], NORM_FACTOR_SPECTRUM_Prot)
    Prot[key]["Error"] = MeVtokRad_2D(Prot[key]["Error"], NORM_FACTOR_SPECTRUM_Prot)

# print(Elec["QGSP_BERT_HP-Elec_463330_452375.csv"]["Dose"])

# print(Prot["QGSP_BIC-Prot_9412_1456.csv"])

plt.figure(1)
for key in Elec:
    plt.barh(key, Elec[key]["Dose"])
    print(Elec[key]["Dose"], key)

plt.figure(2)
for key in Prot:
    plt.barh(key, Prot[key]["Dose"])
    print(Prot[key]["Dose"], key)

plt.show()
'''
PhysList = [
"binary_had.mac",
"em_standard_opt4.mac",
"binary_ion.mac",
"em_standard_remizovich.mac",
"QBBC.mac",
"decay.mac",
"em_standard_space.mac",
"QGSP_BERT_HP.mac",
"Default.mac",
"em_standardSS.mac",
"QGSP_BERT.mac",
"elastic.mac",
"em_standardWVI.mac",
"QGSP_BIC_HP.mac",
"firsov.mac",
"QGSP_BIC.mac",
"em_livermore.mac",
"FTFP_BERT.mac",
"QGSP.mac",
"em_penelope.mac",
"gamma_nuc.mac",
"QGSP_QMD_HP.mac",
"em_standard.mac",
"jqmd_ion.mac",
"raddecay.mac",
"em_standardNR.mac",
"LHEP.mac",
"em_standard_opt1.mac",
"em_standard_opt2.mac",
"Shielding.mac",
"em_standard_opt3.mac",
"stopping.mac"]

print(PhysList)

