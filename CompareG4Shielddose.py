import matplotlib.pyplot as plt
from Dependencies.ReadSD2Q import readSDQ2


# ------------------------------- Import and Plot SHIELDOSE Data -------------------------------------------------------
SDData = readSDQ2("Dependencies/spenvis_sqo.txt")
# SDDataCollumns = ['Aluminium Thickness', 'Total Dose', 'Electrons', 'Bremsstrahlung', 'Protons']
plt.plot(SDData[:, 0], (SDData[:, 2] + SDData[:, 3]) / 1000, '-.', label="SHIELDOSE-2Q trapped Electrons")
plt.plot(SDData[:, 0], SDData[:, 4] / 1000, '--', label="SHIELDOSE-2Q trapped Protons")
plt.plot(SDData[:, 0], SDData[:, 1] / 1000, '-', label="SHIELDOSE-2Q Total trapped Dose")

####### Plot 10kRad line #########
CriticalDose = [10] * SDData.shape[0]
plt.plot(SDData[:, 0], CriticalDose, color='k', linewidth=2, label='Critical Dose of 10 krad')


FolderName = "MulasElectron2e9AlFull"

# ####### Import and Plot G4 Data #########
G4Data = G4TotalDose("FolderName", 2e9, 8.003046E+14)
plt.plot(G4Data[0], G4Data[2], 'o', label="G4 Electrons 2e9 Full Spectrum")
plt.plot(G4Data[0], G4Data[4], 'o', label="G4 Electrons 2e9 Full Spectrum Secondary")

# plt.xlim([0, 17])
# plt.ylim([0.2, 2e2])
plt.yscale("log")
# plt.xscale("log")
plt.grid(which='both')
plt.title("Trapped particle spectra shielded by aluminium")
plt.xlabel("Aluminium Absorber Thickness [mm]")
plt.ylabel("Dose in Si [krad]")
plt.legend()
plt.show()
#plt.savefig("../" + FolderName + "/" + FolderName + "TotalDoseCurve.eps", format='eps')