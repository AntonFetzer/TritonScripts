from Dependencies.MergeTotalDose import mergeTotalDose
from Dependencies.TotalDose import totalDose
import numpy as np
from uncertainties import ufloat

Materials = ["Al-7075", "G4_POLYETHYLENE", "G4_W", "FR4"]


print("Number of Names:", len(Materials))

Path = "/l/triton_work/Permutations/5Layer/Res/"
file_name = Path + "../Analysis/5Layer-Raw.csv"

Electrons = totalDose(Path, filename_contains="Elec")
Protons = totalDose(Path, filename_contains="Prot")

print("Electrons Shape:", np.shape(Electrons["dose"]))

Total = mergeTotalDose([Electrons, Protons])
# print(Total)

NumMat = len(Materials)

#for i1 in range(NumMat):
#    for i2 in range(NumMat):
#        for i3 in range(NumMat):
#            for i4 in range(NumMat):
#                for i5 in range(NumMat):
#                    i = i1 * NumMat * NumMat * NumMat * NumMat + i2 * NumMat * NumMat * NumMat + i3 * NumMat * NumMat + i4 * NumMat + i5
#                    print(i+1, Materials[i1], Materials[i2], Materials[i3], Materials[i4], Materials[i5], ufloat(Electrons["dose"][i], Electrons["error"][i]), ufloat(Protons["dose"][i], Protons["error"][i]), ufloat(Total["dose"][i], Total["error"][i]))

with open(file_name, 'w') as file:
    file.write("Combination #,Material 1 Z-Number,Material 2 Z-Number,Material 3 Z-Number,Material 4 Z-Number,Material 5 Z-Number,Material 1,Material 2,Material 3,Material 4,Material 5,Electron Dose [krad/Month],Electron Err [krad/Month],Proton Dose [krad/Month],Proton Err [krad/Month],Total Dose [krad/Month],Total Err [krad/Month]\n")
    for i1 in range(NumMat):
        for i2 in range(NumMat):
            for i3 in range(NumMat):
                for i4 in range(NumMat):
                    for i5 in range(NumMat):
                        i = i1 * NumMat**4 + i2 * NumMat**3 + i3 * NumMat**2 + i4 * NumMat + i5
                        line = f"{i+1},{i1+1},{i2+1},{i3+1},{i4+1},{i5+1},{Materials[i1]},{Materials[i2]},{Materials[i3]},{Materials[i4]},{Materials[i5]},{ufloat(Electrons['dose'][i], Electrons['error'][i])},{ufloat(Protons['dose'][i], Protons['error'][i])},{ufloat(Total['dose'][i], Total['error'][i])}\n"
                        line = line.replace("+/-", ",")
                        file.write(line)
