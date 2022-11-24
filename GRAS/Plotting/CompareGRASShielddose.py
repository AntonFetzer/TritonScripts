import numpy as np
import matplotlib.pyplot as plt
from GRAS.Read.ReadSD2Q import readSDQ2
from GRAS.Dependencies.TotalKRadGras import totalkRadGras

Path = "/home/anton/Desktop/triton_work/ShieldingCurves/"
Folder = "SuperGTO"
ElecMeV = 0.5
ProtMeV = 10

Electrons = totalkRadGras(Path + Folder + "/Res/", "Elec")
Protons = totalkRadGras(Path + Folder + "/Res/", "Prot")
Sol = totalkRadGras(Path + Folder + "/Res/", "Sol")

Sol = Sol / (30*24*60*60) # Fluence was in per month. Not per second.

x = np.linspace(0, 2.5/0.27, num=101, dtype=float)
#x = np.linspace(0.05, 2.5, num=50, dtype=float)
#for i in x:
#    print(i)

fig1 = plt.figure(1)
# ------------------------------- Import and Plot SHIELDOSE Data -------------------------------------------------------
SDData = readSDQ2("/home/anton/Desktop/triton_work/SuperGTO/spenvis_sqo.txt")
# SDDataCollumns = ['Aluminium Thickness', 'Total Dose', 'Electrons', 'Bremsstrahlung', 'Protons']
plt.plot(SDData[:, 0], (SDData[:, 2] + SDData[:, 3]) / 1000, '-.', label="SHIELDOSE-2Q Electrons", linewidth=1)
plt.plot(SDData[:, 0], SDData[:, 4] / 1000, '--', label="SHIELDOSE-2Q Protons", linewidth=1)
plt.plot(SDData[:, 0], SDData[:, 5] / 1000, '--', label="SHIELDOSE-2Q Solar Protons", linewidth=1)
plt.plot(SDData[:, 0], SDData[:, 1] / 1000, '-', label="SHIELDOSE-2Q Total Dose", linewidth=1)

####### Plot 10kRad line #########
CriticalDose = [1] * SDData.shape[0]
plt.plot(SDData[:, 0], CriticalDose, color='k', linewidth=2, label='1 krad per month')

# ----------------------------------------- Plot electrons500kev ---------------------------------------------------------
plt.errorbar(x, Electrons[0], Electrons[1], fmt='C0 ', capsize=4, label="GRAS Trapped Electrons > " + str(ElecMeV) + " MeV", elinewidth=1, capthick=1)

# ----------------------------------------- Plot protons10MeV -------------------------------------------------------
plt.errorbar(x, Protons[0], Protons[1], fmt='C1 ', capsize=4, label="GRAS Trapped Protons > " + str(ProtMeV) + " MeV", elinewidth=1, capthick=1)

# ----------------------------------------- Plot Solar protons -------------------------------------------------------
plt.errorbar(x, Sol[0], Sol[1], fmt='C2 ', capsize=4, label="GRAS Solar Protons > " + str(ProtMeV) + " MeV", elinewidth=1, capthick=1)

Total = Electrons
Total[0] = Electrons[0] + Protons[0] + Sol[0]
Total[1] = np.sqrt(Electrons[1] * Electrons[1] + Protons[1] * Protons[1] + Sol[1] * Sol[1])

# ----------------------------------------- Plot Total -------------------------------------------------------
plt.errorbar(x, Total[0], Total[1], fmt='C3 ', capsize=4, label="GRAS Total Ionizing Dose", elinewidth=1, capthick=1)


plt.xlim(-0.2, 2.5/0.27)
# plt.ylim([0.2, 2e2])
plt.yscale("log")
# plt.xscale("log")
plt.grid(which='major')
plt.title(Folder)
plt.xlabel("Aluminium Shield Thickness [mm]")
plt.ylabel("Total ionizing dose in silicon [krad]")
plt.legend()
plt.savefig("/home/anton/Desktop/triton_work/ShieldingCurves/SuperGTO/CompareGRASSHieldose-mm.svg", format='svg', bbox_inches="tight")


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