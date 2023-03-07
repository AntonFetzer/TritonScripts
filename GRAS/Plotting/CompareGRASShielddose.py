import numpy as np
import matplotlib.pyplot as plt
from GRAS.Read.ReadSD2Q import readSDQ2
from GRAS.Dependencies.TotalKRadGras import totalkRadGras

Path = "/home/anton/Desktop/triton_work/ShieldingCurves"

Electrons = totalkRadGras(Path + "/AE9-500keV-10mm/Res/", "Elec")
Protons = totalkRadGras(Path + "/AP9-10MeV-10mm/Res/", "Prot")

x = np.linspace(0, 10, num=101, dtype=float)

for i in x:
    print(i)

fig1 = plt.figure(1, [6.4, 4])
# ------------------------------- Import and Plot SHIELDOSE Data -------------------------------------------------------
SDData = readSDQ2("/home/anton/Desktop/triton_work/Spectra/A9/Shieldose/spenvis_sqo.txt")
# SDDataCollumns = ['Aluminium Thickness', 'Total Dose', 'Electrons', 'Bremsstrahlung', 'Protons']
plt.plot(SDData[:, 0], SDData[:, 2] / 1000, 'C0:.', label="SHIELDOSE-2Q Electrons")
plt.plot(SDData[:, 0], SDData[:, 3] / 1000, 'C1-..', label="SHIELDOSE-2Q Bremsstrahlung")
plt.plot(SDData[:, 0], SDData[:, 4] / 1000, 'C3--.', label="SHIELDOSE-2Q Protons")
plt.plot(SDData[:, 0], SDData[:, 1] / 1000, 'C2-.', label="SHIELDOSE-2Q Total Dose")

####### Plot 10kRad line #########
CriticalDose = [1] * SDData.shape[0]
plt.plot(SDData[:, 0], CriticalDose, color='k', linewidth=2, label='1 krad per month')

# ----------------------------------------- Plot electrons500kev ---------------------------------------------------------
plt.errorbar(x, Electrons[0], Electrons[1], fmt='C0 ', capsize=4, label="GRAS Trapped Electrons", elinewidth=1, capthick=1)

# ----------------------------------------- Plot protons10MeV -------------------------------------------------------
plt.errorbar(x, Protons[0], Protons[1], fmt='C3 ', capsize=4, label="GRAS Trapped Protons", elinewidth=1, capthick=1)


Total = Electrons
Total[0] = Electrons[0] + Protons[0]
Total[1] = np.sqrt(Electrons[1] * Electrons[1] + Protons[1] * Protons[1] )

# ----------------------------------------- Plot Total -------------------------------------------------------
plt.errorbar(x, Total[0], Total[1], fmt='C2 ', capsize=4, label="GRAS Total Ionizing Dose", elinewidth=1, capthick=1)


plt.xlim(-0.2, 10.2)
plt.ylim([5e-2, 5e3])
plt.yscale("log")
#plt.xscale("log")
plt.grid(which='major')
plt.title("TID from energetic trapped particles (AE9/AP9) on GTO")
plt.xlabel("Aluminium Shield Thickness [mm]")
plt.ylabel("Total ionizing dose in silicon [krad per month]")
plt.legend()
plt.savefig(Path + "/Plots/A9-GTO-10mm.pdf", format='pdf', bbox_inches="tight")
'''


fig2 = plt.figure(2)
plt.title(Folder)
plt.plot(x, Electrons[3], 'C0.', label="GRAS Electrons > " + str(ElecMeV) + " MeV")
plt.plot(x, Protons[3], 'C1.', label="GRAS Protons > " + str(ProtMeV) + " MeV")
plt.legend()
plt.xlabel("Aluminium Shield Thickness [g/cm2]")
plt.ylabel("Number of Non Zero Entries")
plt.grid(which='both')
plt.yscale("log")
#plt.savefig(Path + Folder + "/NonZeroEntries.eps", format='eps', bbox_inches="tight")

fig3 = plt.figure(3)
plt.title(Folder)
plt.plot(x, 100 * Electrons[1] / Electrons[0], 'C0.', label="GRAS Electrons > " + str(ElecMeV) + " MeV")
plt.plot(x, 100 * Protons[1] / Protons[0], 'C1.', label="GRAS Protons > " + str(ProtMeV) + " MeV")
plt.legend()
plt.xlabel("Aluminium Shield Thickness [g/cm2]")
plt.ylabel("Relative Error in %")
plt.grid(which='both')
#plt.yscale("log")
#plt.savefig(Path + Folder + "/RelativeError.eps", format='eps', bbox_inches="tight")

plt.show()
'''

