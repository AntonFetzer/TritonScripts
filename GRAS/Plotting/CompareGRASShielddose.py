import numpy as np
import matplotlib.pyplot as plt
from GRAS.Read.ReadSD2Q import readSDQ2
from GRAS.Dependencies.TotalKRadGras import totalkRadGras
from matplotlib.ticker import LogLocator
from matplotlib.ticker import LogLocator, MultipleLocator

Path = "/home/anton/Desktop/triton_work/ShieldingCurves"

Electrons = totalkRadGras(Path + "/05Elec-10Prot/Res/", "Elec")
Protons = totalkRadGras(Path + "/05Elec-10Prot/Res/", "Prot")

x = np.linspace(0, 2.5, num=101, dtype=float)

# for i in x:
#    print(i)

# fig1 = plt.figure(1, [6.4, 4])
# ------------------------------- Import and Plot SHIELDOSE Data -------------------------------------------------------
ConvertKrad = 1 / 1000  #  1/1000 to convert to krad
mmtogcm2 = 0.27  # 0.27 to convert mm to g/cm2 of aliminium

SDData = readSDQ2("/home/anton/Desktop/triton_work/Spectra/A9/Shieldose/spenvis_sqo.txt")
# SDDataCollumns = ['Aluminium Thickness', 'Total Dose', 'Electrons', 'Bremsstrahlung', 'Protons']
plt.plot(SDData[:, 0] * mmtogcm2, ( SDData[:, 2]+SDData[:, 3] ) * ConvertKrad, 'C0--', label="SHIELDOSE-2Q Electron Dose")
#plt.plot(SDData[:, 0] * mmtogcm2, SDData[:, 3] * ConvertKrad, 'C1-..', label="SHIELDOSE-2Q Bremsstrahlung")
plt.plot(SDData[:, 0] * mmtogcm2, SDData[:, 4] * ConvertKrad, 'C1-', label="SHIELDOSE-2Q Proton Dose")
plt.plot(SDData[:, 0] * mmtogcm2, SDData[:, 1] * ConvertKrad, 'C2:', label="SHIELDOSE-2Q Total Dose")

####### Plot 10kRad line #########
CriticalDose = [1] * 101
plt.plot(x, CriticalDose, color='k', linewidth=2, label='1 krad per month')

# ----------------------------------------- Plot electrons500kev ---------------------------------------------------------
plt.errorbar(x, Electrons[0], Electrons[1], fmt='C0 ', capsize=4, label="GRAS Electron Dose", elinewidth=1,
             capthick=1)

# ----------------------------------------- Plot protons10MeV -------------------------------------------------------
plt.errorbar(x, Protons[0], Protons[1], fmt='C1 ', capsize=4, label="GRAS Proton Dose", elinewidth=1, capthick=1)

Total = Electrons
Total[0] = Electrons[0] + Protons[0]
Total[1] = np.sqrt(Electrons[1] * Electrons[1] + Protons[1] * Protons[1])

# ----------------------------------------- Plot Total -------------------------------------------------------
plt.errorbar(x, Total[0], Total[1], fmt='C2 ', capsize=4, label="GRAS Total Ionizing Dose", elinewidth=1, capthick=1)

plt.xlim(-0.05, 2.55)
plt.ylim([5e-2, 5e3])
plt.yscale("log")
# plt.xscale("log")

ax = plt.gca()
ax.xaxis.set_minor_locator(MultipleLocator(0.25))  # Set minor ticks every 0.1 units
ax.yaxis.set_minor_locator(LogLocator(subs='all', numticks=10))

plt.grid(which='both', linestyle='-', linewidth=0.5)
plt.title("Ionising dose from trapped particles (AE9/AP9) on GTO", fontsize=14)
plt.xlabel("Aluminium Shielding Depth [g/cm2]", fontsize=12)
plt.ylabel("Ionizing dose in silicon [krad per month]", fontsize=12)
plt.legend(fontsize=12)
plt.savefig("/home/anton/Desktop/TritonPlots/Paper/ShielddoseComparison.pdf", format='pdf', bbox_inches="tight")
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
'''
# plt.show()
