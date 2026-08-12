from Dependencies.MergeTotalDose import mergeTotalDose
from Dependencies.TotalDose import totalDose
import numpy as np
from uncertainties import ufloat

Materials = ["G4_Al", "Al_7075", "G4_POLYETHYLENE", "G4_KEVLAR", "G4_Pb", "G4_W", "G4_STAINLESS_STEEL", "CarbonFibre", "FR4", "G4_Ta", "G4_TEFLON", "Ti_6AL_4V"]

#Densities = [2, 0.534, 1.55, 2.33, 1.74, 2.37, 2.699, 4.54, 11.35, 19.3]

print("Number of Names:", len(Materials))
#print("Number of Densities:", len(Densities))

#VolumesStr = ''
#for i, Name in enumerate(Names):
#    VolumesStr += '        <volume name ="ShieldVol_' + str(i) + '">\n            <materialref ref="' + Name + '"/>\n            <solidref ref="Shield_' + str(i) + '"/>\n        </volume>\n\n'
#print(VolumesStr)

#for x in Densities:
#    print(x)

Path = "/l/triton_work/Permutations/3Layer/Res/"
file_name = Path + "../Analysis/3Layer-Raw.csv"

Electrons = totalDose(Path, filename_contains="Elec")
Protons = totalDose(Path, filename_contains="Prot")

print("Electrons Shape:", np.shape(Electrons["dose"]))

Total = mergeTotalDose([Electrons, Protons])
# print(Total)

NumMat = len(Materials)

#for i1 in range(NumMat):
#    for i2 in range(NumMat):
#        for i3 in range(NumMat):
#            i = i1 * NumMat * NumMat + i2 * NumMat + i3
#            print(i+1, Materials[i1], Materials[i2], Materials[i3], ufloat(Electrons["dose"][i], Electrons["error"][i]), ufloat(Protons["dose"][i], Protons["error"][i]), ufloat(Total["dose"][i], Total["error"][i]))

with open(file_name, 'w') as file:
    file.write("Combination #,Material 1 Z-Number,Material 2 Z-Number,Material 3 Z-Number,Material 1,Material 2,Material 3,Electron Dose [krad/Month],Electron Err [krad/Month],Proton Dose [krad/Month],Proton Err [krad/Month],Total Dose [krad/Month],Total Err [krad/Month]\n")
    for i1 in range(NumMat):
        for i2 in range(NumMat):
            for i3 in range(NumMat):
                i = i1 * NumMat**2 + i2 * NumMat + i3
                line = f"{i+1},{i1+1},{i2+1},{i3+1},{Materials[i1]},{Materials[i2]},{Materials[i3]},{ufloat(Electrons['dose'][i], Electrons['error'][i])},{ufloat(Protons['dose'][i], Protons['error'][i])},{ufloat(Total['dose'][i], Total['error'][i])}\n"
                line = line.replace("+/-", ",")
                file.write(line)
