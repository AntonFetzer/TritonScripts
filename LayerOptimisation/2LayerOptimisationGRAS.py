import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from Dependencies.HistogramPlots import savePdf
from Dependencies.MergeTotalDose import mergeTotalDose
from Dependencies.TotalDose import totalDose
from uncertainties import ufloat

Path = "/scratch/work/fetzera1/GRAS/2D/2LayerOptimisation"
ShieldingDepth = "0.564"  # g/cm2

MatA = "Polyethylene"
MatB = "Lead"
A = "PE"
B = "Pb"
Shield = A + "-" + B
Folder = Shield + "-0.564gcm2-FlatMissionGTO"  # Study folder; equal to Shield for the original runs

# Biased GPS sources (/gps/hist/type biase): GRAS 5 overestimates the quoted
# errors, so take the error from the replica scatter instead.
ReplicaError = True
Months = None  # Divide all doses by this number of months, e.g. 132 for an 11-year fluence

# Optional: python3 -m LayerOptimisation.2LayerOptimisationGRAS <study folder> [--months=N]
Args = [a for a in sys.argv[1:] if not a.startswith("--")]
if Args:
    Folder = Args[0]  # e.g. PE-Pb-0.564gcm2-FS2-A9-flatlog
for a in sys.argv[1:]:
    if a.startswith("--months="):
        Months = float(a.split("=", 1)[1])

Ymax = 1  # Max kRad shown in plots so that every plot has the same scale

Path = os.path.join(Path, Folder) + "/"
ResultsFolder = "Res/"
ElecA = totalDose(Path + ResultsFolder, filename_contains="Electrons2MatA", replica_error=ReplicaError)
ProtA = totalDose(Path + ResultsFolder, filename_contains="Protons2MatA", replica_error=ReplicaError)

ElecB = totalDose(Path + ResultsFolder, filename_contains="Electrons2MatB", replica_error=ReplicaError)
ProtB = totalDose(Path + ResultsFolder, filename_contains="Protons2MatB", replica_error=ReplicaError)

PerUnit = ""
if Months:
    PerUnit = " per month"
    for d in (ElecA, ProtA, ElecB, ProtB):
        d["dose"] = d["dose"] / Months
        d["error"] = d["error"] / Months

x = np.linspace(0, 100, num=101, dtype=int)

print(np.shape(ElecA["dose"]))
print(np.shape(ElecB["dose"]))
print(np.shape(ProtA["dose"]))
print(np.shape(ProtB["dose"]))

ElecB = {key: np.flip(values) for key, values in ElecB.items()}
ProtB = {key: np.flip(values) for key, values in ProtB.items()}


plt.figure(1)
plt.errorbar(x, ElecA["dose"], ElecA["error"], fmt=' ', capsize=2, label="Electrons " + A + " on " + B)
plt.errorbar(x, ElecB["dose"], ElecB["error"], fmt=' ', capsize=2, label="Electrons " + B + " on " + A)

plt.errorbar(x, ProtA["dose"], ProtA["error"], fmt=' ', capsize=2, label="Protons " + A + " on " + B)
plt.errorbar(x, ProtB["dose"], ProtB["error"], fmt=' ', capsize=2, label="Protons " + B + " on " + A)

#plt.ylim(0, Ymax)
plt.title("Dose deposited by trapped particles in 0.5 mm Si \n behind " + ShieldingDepth + " g/cm2 of " + MatA + "-" + MatB + " shielding")
plt.xlabel("Percentage of shielding mass in " + MatA + " [%]")
plt.ylabel("Deposited ionising dose [kRad" + PerUnit + "]")
plt.grid(which='both')
plt.legend()
# plt.yscale("log")
savePdf(Path + Shield + "-Gradient.pdf")
plt.show()

TotalA = mergeTotalDose([ElecA, ProtA])

TotalB = mergeTotalDose([ElecB, ProtB])

plt.figure(2)
plt.errorbar(x, TotalA["dose"], TotalA["error"], fmt=' ', capsize=2, label=MatA + " on top of " + MatB)
plt.errorbar(x, TotalB["dose"], TotalB["error"], fmt=' ', capsize=2, label=MatB + " on top of " + MatA)

#plt.ylim(0, Ymax)
plt.title("Total dose deposited by trapped particles in 0.5 mm Si \n behind " + ShieldingDepth + " g/cm2 of " + MatA + "-" + MatB + " shielding")  # --------
plt.xlabel("Percentage of shielding mass in " + MatA + " [%]")
plt.ylabel("Deposited ionising dose [krad" + PerUnit + "]")
plt.grid(which='both')
plt.legend(loc='upper left')
# plt.yscale("log")
savePdf(Path + Shield + "-GradientSum.pdf")
plt.show()


TotalAmin = np.min(TotalA["dose"])
TotalBmin = np.min(TotalB["dose"])
TotalAminIndex = np.argmin(TotalA["dose"])
TotalBminIndex = np.argmin(TotalB["dose"])

ElecAmin = ElecA["dose"][TotalAminIndex]
ElecBmin = ElecB["dose"][TotalBminIndex]

ElecAminErr = ElecA["error"][TotalAminIndex]
ElecBminErr = ElecB["error"][TotalBminIndex]

ProtAmin = ProtA["dose"][TotalAminIndex]
ProtBmin = ProtB["dose"][TotalBminIndex]

ProtAminErr = ProtA["error"][TotalAminIndex]
ProtBminErr = ProtB["error"][TotalBminIndex]

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
