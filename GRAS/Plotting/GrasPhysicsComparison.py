import os
import csv
import numpy as np
from ReadGRASCSV1Tile import readGrasCsv1Tile
from GRAS.Dependencies.MeVtokRad2DGras import MeVtokRad_2D
import matplotlib.pyplot as plt


path = "/l/triton_work/1Tile/PhysicsComparison"

ElecFiles = [f for f in os.listdir(path + "/Results/") if "Elec" in f]
ProtFiles = [f for f in os.listdir(path + "/Results/") if "Prot" in f]

ElecFiles = sorted(ElecFiles, key=str.upper)
ProtFiles = sorted(ProtFiles, key=str.upper)

Elec = dict.fromkeys(ElecFiles, {})
Prot = dict.fromkeys(ProtFiles, {})

NORM_FACTOR_SPECTRUM_Elec = 7.891281E+14
NORM_FACTOR_SPECTRUM_Prot = 3.389664E+11

for key in Elec:
    Elec[key] = readGrasCsv1Tile(path + "/Results/" + key)
    Elec[key]["Dose"] = MeVtokRad_2D(Elec[key]["Dose"], NORM_FACTOR_SPECTRUM_Elec)
    Elec[key]["Error"] = MeVtokRad_2D(Elec[key]["Error"], NORM_FACTOR_SPECTRUM_Elec)

for key in Prot:
    Prot[key] = readGrasCsv1Tile(path + "/Results/" + key)
    Prot[key]["Dose"] = MeVtokRad_2D(Prot[key]["Dose"], NORM_FACTOR_SPECTRUM_Prot)
    Prot[key]["Error"] = MeVtokRad_2D(Prot[key]["Error"], NORM_FACTOR_SPECTRUM_Prot)

# print(Elec["QGSP_BERT_HP-Elec_463330_452375.csv"]["Dose"])

# print(Prot["QGSP_BIC-Prot_9412_1456.csv"])

print("Electron Dose")
plt.figure(1)
for key in Elec:
    plt.barh(key, Elec[key]["Dose"])
    plt.title("Electron Dose")
    print(Elec[key]["Dose"])

print("Proton Dose")
plt.figure(2)
for key in Prot:
    plt.barh(key, Prot[key]["Dose"])
    plt.title("Proton Dose")
    print(Prot[key]["Dose"])

print("Electron Entries")
plt.figure(3)
for key in Elec:
    plt.barh(key, Elec[key]["Entries"])
    plt.title("Electron Entries")
    print(Elec[key]["Entries"])

print("Proton Entries")
plt.figure(4)
for key in Prot:
    plt.barh(key, Prot[key]["Entries"])
    plt.title("Proton Entries")
    print(Prot[key]["Entries"])

#plt.show()

print("Number of Proton Datasets:", len(Prot))
print("Number of Electron Datasets:", len(Elec))

PhysList = [
'bertini',
'bertini_hp',
'bertini_preco',
'binary',
'binary_had',
'binary_hp',
'binary_ion',
'decay',
'Default',
'elastic',
'elasticCHIPS',
'elasticHP',
'em_livermore',
'em_penelope',
'em_standard',
'em_standardNR',
'em_standardSS',
'em_standardWVI',
'em_standard_opt1',
'em_standard_opt2',
'em_standard_opt3',
'em_standard_opt4',
'em_standard_remizovich',
'em_standard_space',
'firsov',
'FTFP_BERT',
'gamma_nuc',
'incl_ion',
'jqmd_ion',
'LHEP',
'QBBC',
'QGSP',
'QGSP_BERT',
'QGSP_BERT_HP',
'QGSP_BIC',
'QGSP_BIC_HP',
'QGSP_QMD_HP',
'qmd_ion',
'raddecay',
'secondary_generator',
'Shielding',
'stopping']

print(len(PhysList))

PhysList = sorted(PhysList, key=str.upper)
#print(PhysList)

for i, Str in enumerate(PhysList):
    print("  " + str(i) + ")  MAC=" + Str + ".mac ;;")


#Path = "/l/Test/"
#
#for Str in PhysList:
#    with open(Path + Str + '.mac', 'w') as f:
#        f.write('/control/alias ThisPhysList ' + Str + '\n')
#        f.write('\n')
#        f.write("/control/execute Physcomp.mac")


LogFiles = [f for f in os.listdir(path + "/log/") if ".log" in f]

#LogFiles = sorted(LogFiles)

Times = np.zeros(len(LogFiles), dtype=float)
Phys = []

for i, file in enumerate(LogFiles):
    with open(path + "/log/" + file, 'r') as f:
        reader = csv.reader(f)
        for line in reader:
            if line:
                if 'Total CPU time used: ' in line[0] and Times[i] == 0:
                    Times[i] = - 1
                elif 'Total CPU time used: ' in line[0] and Times[i] < 0:
                    Times[i] = float(line[0][52:59])
                elif '/gras/physics/addPhysics ' in line[0] and '###' not in line[0]:
                    Phys.append(line[0][25:])

plt.figure(5)
for i in range(len(Times)):
    print(Times[i], Phys[i])
    plt.barh(Phys[i], Times[i])
    plt.title("Times")

plt.show()
