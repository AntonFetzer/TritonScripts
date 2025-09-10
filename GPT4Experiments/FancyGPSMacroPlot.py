import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
from GRAS.Read.ReadGPSMacro import readGPSMacro

file = "/l/triton_work/Spectra/GTO/AE9/AE9500keV.mac"
Electrons = readGPSMacro(file)

file = "/l/triton_work/Spectra/GTO/AP9/AP910MeV.mac"
Protons = readGPSMacro(file)

plt.plot(Electrons['Energy'], Electrons['Flux'], linestyle='--', linewidth=2, label="AE-9 Electron Flux")
plt.plot(Protons['Energy'], Protons['Flux'], linestyle='-', linewidth=2, label="AP-9 Proton Flux")

plt.yscale("log")
plt.xscale("log")

ax = plt.gca()
ax.xaxis.set_minor_locator(LogLocator(subs='all', numticks=10))
ax.yaxis.set_minor_locator(LogLocator(subs='all', numticks=10))

ax.grid(which='both', linestyle='-', linewidth=0.5)

plt.legend(fontsize=12)
plt.title("Differential AP-9 and AE-9 spectra on GTO", fontsize=14)
plt.xlabel("Kinetic energy [MeV]", fontsize=12)
plt.ylabel("Differential Flux [cm-2 s-1 MeV-1]", fontsize=12)

plt.show()
# plt.savefig("/l/TritonPlots/Paper/Spectra.pdf", format='pdf', bbox_inches="tight")

