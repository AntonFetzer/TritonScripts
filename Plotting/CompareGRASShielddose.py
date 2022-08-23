import numpy as np
import matplotlib.pyplot as plt
from Dependencies.ReadSD2Q import readSDQ2
from natsort import natsorted
import os
import csv
import sys
# setting path
sys.path.append('../GRAS')
from TotalKRadGras import totalkRadGras

Path = "/home/anton/Desktop/triton_work/ShieldingCurves/100AlFull/Res/"

Electrons = totalkRadGras(Path, "Elec")
Protons = totalkRadGras(Path, "Prot")
x = np.linspace(0, 2.5, num=101, dtype=float)
#x = np.linspace(0.05, 2.5, num=50, dtype=float)
#for i in x:
#    print(i)


# ------------------------------- Import and Plot SHIELDOSE Data -------------------------------------------------------
SDData = readSDQ2("../Dependencies/spenvis_sqoA9gcm2.txt")
# SDDataCollumns = ['Aluminium Thickness', 'Total Dose', 'Electrons', 'Bremsstrahlung', 'Protons']
plt.plot(SDData[:, 0], (SDData[:, 2] + SDData[:, 3]) / 1000, '-.', label="SHIELDOSE-2Q Electrons")
plt.plot(SDData[:, 0], SDData[:, 4] / 1000, '--', label="SHIELDOSE-2Q Protons")
plt.plot(SDData[:, 0], SDData[:, 1] / 1000, '-', label="SHIELDOSE-2Q Total Dose")

####### Plot 10kRad line #########
CriticalDose = [1] * SDData.shape[0]
plt.plot(SDData[:, 0], CriticalDose, color='k', linewidth=2, label='Critical Dose of 1 krad per month')

# ----------------------------------------- Plot electrons500kev ---------------------------------------------------------
plt.errorbar(x, Electrons[0], Electrons[1], fmt='C0 ', capsize=3, label="GRAS Electrons Full")

# ----------------------------------------- Plot electronsfull ---------------------------------------------------------
#plt.errorbar(ThickList, Data[:, 1, 2], Data[:, 1, 3], fmt='bD', capsize=7, markersize=5, label="GRAS Electrons > 40 keV")

# ----------------------------------------- Plot protons10MeV -------------------------------------------------------
plt.errorbar(x, Protons[0], Protons[1], fmt='C1 ', capsize=3, label="GRAS Protons Full")

# ----------------------------------------- Plot protonsfull -------------------------------------------------------
#plt.errorbar(ThickList, Data[:, 3, 2], Data[:, 3, 3], fmt='yo', capsize=7, markersize=5, label="GRAS Protons > 100 keV")

Total = Electrons
Total[0] = Electrons[0] + Protons[0]
Total[1] = np.sqrt(Electrons[1] * Electrons[1] + Protons[1] * Protons[1])

# ----------------------------------------- Plot Total -------------------------------------------------------
plt.errorbar(x, Total[0], Total[1], fmt='C2 ', capsize=3, label="GRAS total ionizing dose")


# plt.xlim([0, 17])
# plt.ylim([0.2, 2e2])
plt.yscale("log")
# plt.xscale("log")
plt.grid(which='major')
plt.title("Trapped particle spectra shielded by aluminium")
plt.xlabel("Aluminium Shield Thickness [mm]")
plt.ylabel("Total ionizing dose in silicon [krad]")
plt.legend()
#plt.show()
plt.savefig(Path + "../CompareGRASSHieldose.eps", format='eps', bbox_inches="tight")
