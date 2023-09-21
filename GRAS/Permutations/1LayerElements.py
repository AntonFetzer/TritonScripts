import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat

from GRAS.Dependencies.TotalKRadGras import totalkRadGras

Names = ["G4_lH2", "G4_He", "G4_Li", "G4_Be", "G4_B", "G4_C", "G4_lN2", "G4_lO2", "G4_F", "G4_Ne", "G4_Na", "G4_Mg",
         "G4_Al", "G4_Si", "G4_P", "G4_S", "G4_Cl", "G4_lAr", "G4_K", "G4_Ca", "G4_Sc", "G4_Ti", "G4_V", "G4_Cr",
         "G4_Mn", "G4_Fe", "G4_Co", "G4_Ni", "G4_Cu", "G4_Zn", "G4_Ga", "G4_Ge", "G4_As", "G4_Se", "G4_lBr", "G4_lKr",
         "G4_Rb", "G4_Sr", "G4_Y", "G4_Zr", "G4_Nb", "G4_Mo", "G4_Tc", "G4_Ru", "G4_Rh", "G4_Pd", "G4_Ag", "G4_Cd",
         "G4_In", "G4_Sn", "G4_Sb", "G4_Te", "G4_I", "G4_lXe", "G4_Cs", "G4_Ba", "G4_La", "G4_Ce", "G4_Pr", "G4_Nd",
         "G4_Pm", "G4_Sm", "G4_Eu", "G4_Gd", "G4_Tb", "G4_Dy", "G4_Ho", "G4_Er", "G4_Tm", "G4_Yb", "G4_Lu", "G4_Hf",
         "G4_Ta", "G4_W", "G4_Re", "G4_Os", "G4_Ir", "G4_Pt", "G4_Au", "G4_Hg", "G4_Tl", "G4_Pb", "G4_Bi", "G4_Po",
         "G4_At", "G4_Rn", "G4_Fr", "G4_Ra", "G4_Ac", "G4_Th", "G4_Pa", "G4_U", "G4_Np", "G4_Pu", "G4_Am", "G4_Cm",
         "G4_Bk", "G4_Cf"]

Densities = [0.0708, 0.000166322, 0.534, 1.848, 2.37, 2, 0.807, 1.141, 0.00158029, 0.000838505, 0.971, 1.74, 2.699,
             2.33, 2.2, 2, 0.00299473, 1.396, 0.862, 1.55, 2.989, 4.54, 6.11, 7.18, 7.44, 7.874, 8.9, 8.902, 8.96,
             7.133, 5.904, 5.323, 5.73, 4.5, 3.1028, 2.418, 1.532, 2.54, 4.469, 6.506, 8.57, 10.22, 11.5, 12.41, 12.41,
             12.02, 10.5, 8.65, 7.31, 7.31, 6.691, 6.24, 4.93, 2.953, 1.873, 3.5, 6.154, 6.657, 6.71, 6.9, 7.22, 7.46,
             5.243, 7.9004, 8.229, 8.55, 8.795, 9.066, 9.321, 6.73, 9.84, 13.31, 16.654, 19.3, 21.02, 22.57, 22.42,
             21.45, 19.32, 13.546, 11.72, 11.35, 9.747, 9.32, 9.32, 0.00900662, 1, 5, 10.07, 11.72, 15.37, 18.95, 20.25,
             19.84, 13.67, 13.51, 14, 10]

print("Number of Names:", len(Names))
print("Number of Densities:", len(Densities))

#VolumesStr = ''
#for i, Name in enumerate(Names):
#    VolumesStr = '        <volume name ="ShieldVol_' + str(i) + '">\n            <materialref ref="' + Name + '"/>\n            <solidref ref="Shield_' + str(i) + '"/>\n        </volume>\n\n'
#    print(VolumesStr)

# for x in Densities:
#    print(x)

Path = "/l/triton_work/Permutations/1Layer/Res/"
file_name = Path + "../Analysis/1Layer-Raw.csv"

Electrons = totalkRadGras(Path, "Elec")
Protons = totalkRadGras(Path, "Prot")

Total = Electrons + Protons
Total[1] = np.sqrt(Electrons[1] * Electrons[1] + Protons[1] * Protons[1])

NumTiles = np.shape(Electrons)[1]
print("NumTiles:", NumTiles)


with open(file_name, 'w') as file:
    file.write("Material 1 Z-Number,Material 1,Electron Dose [krad/Month],Electron Err [krad/Month],Proton Dose [krad/Month],Proton Err [krad/Month],Total Dose [krad/Month],Total Err [krad/Month]\n")
    for i in range(NumTiles):
        line = f"{i+1},{Names[i]},{ufloat(Electrons[0][i], Electrons[1][i])},{ufloat(Protons[0][i], Protons[1][i])},{ufloat(Total[0][i], Total[1][i])}\n"
        line = line.replace("+/-", ",")
        file.write(line)



x = np.linspace(1, NumTiles, num=NumTiles, dtype=int)

#for i in x:
#    print(i, Electrons[0][i-1])
'''
fig1 = plt.figure(1)
plt.errorbar(x, Electrons[0], Electrons[1], fmt='C0', label="AE9 Electrons", linewidth=0.75, capsize=2) #, elinewidth=0.5, capthick=0.5)
plt.errorbar(x, Protons[0], Protons[1], fmt='C1', label="AP9 Protons", linewidth=0.75, capsize=2) #, elinewidth=0.5, capthick=0.5)
plt.errorbar(x, Total[0], Total[1], fmt='C2', label="Total dose", linewidth=0.75, capsize=2) #, elinewidth=0.5, capthick=0.5)
plt.legend()
plt.yscale("log")
plt.yticks([0.25, 0.5, 1, 2], [0.25, 0.5, 1, 2])
plt.grid(which='both')
plt.title("Ionizing dose behind 1.5 g/cm2 of shielding")
plt.ylabel("Total ionizing dose per month in silicon [krad]")
plt.xlabel("Z-Number of shielding material")
plt.show()
#plt.savefig(Path + "../1LayerMaterials.eps", format='eps', bbox_inches="tight")
'''


fig1, ax1 = plt.subplots(figsize=(8, 6))

ax1.errorbar(x, Electrons[0], Electrons[1], fmt='.-', color='C0', label="AE9 Electrons", linewidth=1, capsize=5, elinewidth=1.5, capthick=1.5, alpha=0.8)
ax1.errorbar(x, Protons[0], Protons[1], fmt='.-', color='C1', label="AP9 Protons", linewidth=1, capsize=5, elinewidth=1.5, capthick=1.5, alpha=0.8)
ax1.errorbar(x, Total[0], Total[1], fmt='.-', color='C2', label="Total dose", linewidth=1, capsize=5, elinewidth=1.5, capthick=1.5, alpha=0.8)

ax1.legend()
ax1.set_yscale("log")
ax1.set_yticks([0.25, 0.5, 1, 2])
ax1.set_yticklabels([0.25, 0.5, 1, 2], fontsize=12)
ax1.grid(which='both', linestyle='--', linewidth=0.5)
ax1.set_title("Ionizing dose behind 1.5 g/cm2 of shielding", fontsize=16)
ax1.set_ylabel("Total ionizing dose per month in silicon [krad]", fontsize=14)
ax1.set_xlabel("Z-Number of shielding material", fontsize=14)
ax1.tick_params(axis='both', which='major', labelsize=12)

plt.tight_layout()
plt.show()

