from GRAS.Dependencies.TotalDose import totalkRadGras
import numpy as np
from uncertainties import ufloat

Materials = ["Al-7075", "G4_POLYETHYLENE", "G4_KEVLAR", "G4_W", "G4_STAINLESS-STEEL", "FR4"]

#Densities = [2, 0.534, 1.55, 2.699, 11.35, 19.3]

print("Number of Names:", len(Materials))
#print("Number of Densities:", len(Densities))

#VolumesStr = ''
#for i, Name in enumerate(Names):
#    VolumesStr += '        <volume name ="ShieldVol_' + str(i) + '">\n            <materialref ref="' + Name + '"/>\n            <solidref ref="Shield_' + str(i) + '"/>\n        </volume>\n\n'
#print(VolumesStr)

#for x in Densities:
#    print(x)

Path = "/l/triton_work/Permutations/4Layer/Res/"
file_name = Path + "../Analysis/4Layer-Raw.csv"

Electrons = totalkRadGras(Path, "Elec")
Protons = totalkRadGras(Path, "Prot")

print("Electrons Shape:", np.shape(Electrons))

Total = Electrons + Protons
Total[1] = np.sqrt(Electrons[1] * Electrons[1] + Protons[1] * Protons[1])
# print(Total)

NumMat = len(Materials)

#for i1 in range(NumMat):
#    for i2 in range(NumMat):
#        for i3 in range(NumMat):
#            for i4 in range(NumMat):
#                i = i1 * NumMat * NumMat * NumMat + i2 * NumMat * NumMat + i3 * NumMat + i4
#                print(i+1, Materials[i1], Materials[i2], Materials[i3], Materials[i4], ufloat(Electrons[0][i], Electrons[1][i]), ufloat(Protons[0][i], Protons[1][i]), ufloat(Total[0][i], Total[1][i]))

with open(file_name, 'w') as file:
    file.write("Combination #,Material 1 Z-Number,Material 2 Z-Number,Material 3 Z-Number,Material 4 Z-Number,Material 1,Material 2,Material 3,Material 4,Electron Dose [krad/Month],Electron Err [krad/Month],Proton Dose [krad/Month],Proton Err [krad/Month],Total Dose [krad/Month],Total Err [krad/Month]\n")
    for i1 in range(NumMat):
        for i2 in range(NumMat):
            for i3 in range(NumMat):
                for i4 in range(NumMat):
                    i = i1 * NumMat**3 + i2 * NumMat**2 + i3 * NumMat + i4
                    line = f"{i+1},{i1+1},{i2+1},{i3+1},{i4+1},{Materials[i1]},{Materials[i2]},{Materials[i3]},{Materials[i4]},{ufloat(Electrons[0][i], Electrons[1][i])},{ufloat(Protons[0][i], Protons[1][i])},{ufloat(Total[0][i], Total[1][i])}\n"
                    line = line.replace("+/-", ",")
                    file.write(line)
