import numpy as np
import matplotlib.pyplot as plt
from GRAS.Dependencies.TotalDose import totalkRadGras
from uncertainties import ufloat

Path = "/l/triton_work/2LayerOptimisation"
ShieldingDepth = "1.5"  # g/cm2

MatA = "Polyethylene"
MatB = "Lead"
A = "PE"
B = "Pb"
Shield = A + "-" + B

Ymax = 1  # Max kRad shown in plots so that every plot has the same scale

Path += "/" + Shield + "/"
ResultsFolder = "Res/"
ElecA = totalkRadGras(Path + ResultsFolder, "Electrons2MatA")
ProtA = totalkRadGras(Path + ResultsFolder, "Protons2MatA")

ElecB = totalkRadGras(Path + ResultsFolder, "Electrons2MatB")
ProtB = totalkRadGras(Path + ResultsFolder, "Protons2MatB")

x = np.linspace(0, 100, num=101, dtype=int)

print(np.shape(ElecA))
print(np.shape(ElecB))
print(np.shape(ProtA))
print(np.shape(ProtB))

ElecB = np.flip(ElecB, 1)
ProtB = np.flip(ProtB, 1)


plt.figure(1)
plt.errorbar(x, ElecA[0], ElecA[1], fmt=' ', capsize=2, label="Electrons " + A + " on " + B) # + " Min=" + str(round(np.min(ElecA[0]), 2)) + " kRad at " + str(round(np.argmin(ElecA[0]))) + " % " + A)
plt.errorbar(x, ElecB[0], ElecB[1], fmt=' ', capsize=2, label="Electrons " + B + " on " + A) # + " Min=" + str(round(np.min(ElecB[0]), 2)) + " kRad at " + str(round(np.argmin(ElecB[0]))) + " % " + A)

plt.errorbar(x, ProtA[0], ProtA[1], fmt=' ', capsize=2, label="Protons " + A + " on " + B) # + " Min=" + str(round(np.min(ProtA[0]), 2)) + " kRad at " + str(round(np.argmin(ProtA[0]))) + " % " + A)
plt.errorbar(x, ProtB[0], ProtB[1], fmt=' ', capsize=2, label="Protons " + B + " on " + A) # + " Min=" + str(round(np.min(ProtB[0]), 2)) + " kRad at " + str(round(np.argmin(ProtB[0]))) + " % " + A)

#plt.ylim(0, Ymax)
plt.title("Dose deposited by trapped particles in 0.5 mm Si \n behind " + ShieldingDepth + " g/cm2 of " + MatA + "-" + MatB + " shielding")
plt.xlabel("Percentage of shielding mass in " + MatA + " [%]")
plt.ylabel("Deposited ionising dose [kRad]")
plt.grid(which='both')
plt.legend()
# plt.yscale("log")
plt.show()
# plt.savefig(Path + Shield + "-Gradient.eps", format='eps', bbox_inches="tight")

TotalA = np.zeros_like(ElecA)
TotalA[0] = ElecA[0] + ProtA[0]
TotalA[1] = np.sqrt(ElecA[1] * ElecA[1] + ProtA[1] * ProtA[1])

TotalB = np.zeros_like(ElecB)
TotalB[0] = ElecB[0] + ProtB[0]
TotalB[1] = np.sqrt(ElecB[1] * ElecB[1] + ProtB[1] * ProtB[1])

plt.figure(2)
plt.errorbar(x, TotalA[0], TotalA[1], fmt=' ', capsize=2, label=MatA + " on top of " + MatB) # + " Min=" + str(round(np.min(TotalA[0]), 2)) + " krad at " + str(round(np.argmin(TotalA[0]))) + " % " + A)
plt.errorbar(x, TotalB[0], TotalB[1], fmt=' ', capsize=2, label=MatB + " on top of " + MatA) # + " Min=" + str(round(np.min(TotalB[0]), 2)) + " krad at " + str(round(np.argmin(TotalB[0]))) + " % " + A)

#plt.ylim(0, Ymax)
plt.title("Total dose deposited by trapped particles in 0.5 mm Si \n behind " + ShieldingDepth + " g/cm2 of " + MatA + "-" + MatB + " shielding")  # --------
plt.xlabel("Percentage of shielding mass in " + MatA + " [%]")
plt.ylabel("Deposited ionising dose [krad]")
plt.grid(which='both')
plt.legend()
# plt.yscale("log")
plt.show()
# plt.savefig(Path + Shield + "-GradientSum.eps", format='eps', bbox_inches="tight")


TotalAmin = np.min(TotalA[0])
TotalBmin = np.min(TotalB[0])
TotalAminIndex = np.argmin(TotalA[0])
TotalBminIndex = np.argmin(TotalB[0])

ElecAmin = ElecA[0][TotalAminIndex]
ElecBmin = ElecB[0][TotalBminIndex]

ElecAminErr = ElecA[1][TotalAminIndex]
ElecBminErr = ElecB[1][TotalBminIndex]

ProtAmin = ProtA[0][TotalAminIndex]
ProtBmin = ProtB[0][TotalBminIndex]

ProtAminErr = ProtA[1][TotalAminIndex]
ProtBminErr = ProtB[1][TotalBminIndex]

TotalAminErr = np.sqrt(ElecAminErr * ElecAminErr + ProtAminErr * ProtAminErr)
TotalBminErr = np.sqrt(ElecBminErr * ElecBminErr + ProtBminErr * ProtBminErr)
'''
CSVFile = open(Path + "Results.txt", 'w')
CSVFile.writelines("Material A, Material B, % A, % B, Electron Dose, Electron Error, Proton Dose, Proton Error, Total Dose, Total Error, \n")

print("TotalAmin:", TotalAmin, "TotalBmin:", TotalBmin)

if TotalAmin < TotalBmin:
    List = (A, B, TotalAminIndex, 100-TotalAminIndex, ElecAmin, ElecAminErr, ProtAmin, ProtAminErr, TotalAmin, TotalAminErr)
elif TotalAmin > TotalBmin:
    List = (B, A, 100-TotalBminIndex, TotalBminIndex, ElecBmin, ElecBminErr, ProtBmin, ProtBminErr, TotalBmin, TotalBminErr)

String = ', '.join(map(str, List))
print(String)
CSVFile.writelines(String + "\n")
CSVFile.close()
'''
# Open the CSV file for writing
with open(Path + "Results.txt", 'w') as CSVFile:
    # Write the header to the CSV file
    CSVFile.write(
        "Material A, Material B, % A, % B, Electron Dose, Electron Error, Proton Dose, Proton Error, Total Dose, Total Error, \n")
    print("TotalAmin:", TotalAmin, "TotalBmin:", TotalBmin)

    if TotalAmin < TotalBmin:
        elec_dose = str(ufloat(ElecAmin, ElecAminErr)).split('+/-')
        prot_dose = str(ufloat(ProtAmin, ProtAminErr)).split('+/-')
        total_dose = str(ufloat(TotalAmin, TotalAminErr)).split('+/-')
        List = (A, B, TotalAminIndex, 100 - TotalAminIndex, elec_dose[0], elec_dose[1], prot_dose[0], prot_dose[1],
                total_dose[0], total_dose[1])
    elif TotalAmin > TotalBmin:
        elec_dose = str(ufloat(ElecBmin, ElecBminErr)).split('+/-')
        prot_dose = str(ufloat(ProtBmin, ProtBminErr)).split('+/-')
        total_dose = str(ufloat(TotalBmin, TotalBminErr)).split('+/-')
        List = (B, A, 100 - TotalBminIndex, TotalBminIndex, elec_dose[0], elec_dose[1], prot_dose[0], prot_dose[1],
                total_dose[0], total_dose[1])

    String = ', '.join(map(str, List))
    print(String)
    CSVFile.write(String + "\n")
