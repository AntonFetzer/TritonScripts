from GRAS.Dependencies.TotalDose import totalkRadGras
import numpy as np
import matplotlib.pyplot as plt

Path = "/l/triton_work/ShieldingCurves/Carrington/CarringtonElectronDiffPowTabelated-10mm/Res/"

Data = totalkRadGras(Path)

NumTiles = np.shape(Data['dose'])[0]

fig1 = plt.figure(1)

x = np.linspace(0, 9.626, num=NumTiles, endpoint=True)


Hours = [1, 2, 4, 8, 16, 32, 64, 128]
Fluence = [0.9833, 1.934, 3.741, 7.008, 12.35, 19.54, 26.14, 29.13]

for h, hour in enumerate(Hours):
    Factor = Fluence[h]/(30*24)
    plt.errorbar(x, Data['dose']*Factor, Data['error']*Factor, fmt=' ', capsize=5, label="T=" + str(hour) + " hours")


####### Plot 10kRad line #########
CriticalDose = [10 for i in x]
plt.plot(x, CriticalDose, color='k', linewidth=2, label='10 krad')
CriticalDose = [100 for i in x]
plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='100 krad')

#plt.ylim(5e-2, 1e3)
#plt.ylim(0, 120)
# plt.xlim(-0.5, 20)
plt.grid(which="major")
plt.yscale("log")
plt.title("Carrington Electron Spectrum\nDose deposited behind planar aluminium shielding")
plt.xlabel("Aluminium Shielding Thickness [mm]")
plt.ylabel("Ionising Dose [krad]")
plt.legend()
#plt.savefig("/l/TritonPlots/CarringtonShieldingCurveTlin.eps", format='eps', bbox_inches="tight")

plt.show()
