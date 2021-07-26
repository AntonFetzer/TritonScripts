import numpy as np
import matplotlib.pyplot as plt
from ReadSD2Q import readSDQ2


####### Import and Plot SHIELDOSE Data #########
SDData = readSDQ2("Shieldose/Slab/spenvis_sqo.txt")
# SDDataCollumns = ['Aluminium Thickness', 'Total Dose', 'Electrons', 'Bremsstrahlung', 'Protons']
plt.plot(SDData[:, 0], SDData[:, 4] / 1000, label="SHIELDOSE-2Q trapped Protons")
plt.plot(SDData[:, 0], (SDData[:, 2] + SDData[:, 3]) / 1000, label="SHIELDOSE-2Q trapped Electrons + Bremsstrahlung")

####### Plot 10kRad line #########
CriticalDose = [10] * SDData.shape[0]
plt.plot(SDData[:, 0], CriticalDose, color='k', linewidth=2, label='Critical Dose of 10 krad')

####### Import and Plot MULASSIS Data #########
#MulasData = readMULASSIScsv("GrasCompRAD")
#MulasData1 = readMULASSIScsv("MulasEleTabulatedFull")
#plt.errorbar(MulasData1[0, :], MulasData1[1, :], MulasData1[2, :], capsize=10, fmt='o', label="MULASSIS on SPENVIS")
#MulasData2 = readMULASSIScsv("MulasEleTabulated500keV")
#plt.errorbar(MulasData2[0, :], MulasData2[1, :], MulasData2[2, :], capsize=10, fmt='o', label="MULASSIS on SPENVIS > 500 keV")


####### Import and Plot GRAS Data #########
#GrasData1 = GrasTotalDose("ProtonsTrappedCutoff10MeV", 1e4, 3.381390E+11)
#GrasData1 = GrasTotalDose("ElectronsTrappedCutoff", 1e5, 7.742855E+15)
#plt.plot(GrasData1[0], GrasData1[2], 'd', label="GRAS on SPENVIS")
#GrasData2 = GrasTotalDose("ProtonsTrappedCutoff", 1e6, 8.380532E+14)
#GrasData2 = GrasTotalDose("ElectronsTrappedCutoff500keV", 1e5, 5.886798E+14)
#plt.plot(GrasData2[0], GrasData2[2], 'd', label="GRAS on SPENVIS > 10 MeV")


FolderName = "MulasElectron2e9AlFull"


# ####### Import and Plot G4 Data #########
G4Data = G4TotalDose("FolderName", 2e9, 8.003046E+14)
plt.plot(G4Data[0], G4Data[2], 'o', label="G4 Electrons 2e9 Full Spectrum")
plt.plot(G4Data[0], G4Data[4], 'o', label="G4 Electrons 2e9 Full Spectrum Secondary")

# plt.xlim([0, 17])
# plt.ylim([0.2, 2e2])
plt.yscale("log")
# plt.xscale("log")
plt.grid(which='major')
plt.title("Trapped particle spectra shielded by aluminium")
plt.xlabel("Aluminium Absorber Thickness [mm]")
plt.ylabel("Dose in Si [krad]")
plt.legend()
plt.show()
plt.savefig("../" + FolderName + "/" + FolderName + "TotalDoseCurve.eps", format='eps')



# srun --mem=50G --time=00:15:00 python CompareShieldMulassGras.py