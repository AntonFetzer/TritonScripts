import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import sys, os
# Ensure the package root (the "Python" folder) is on sys.path so "from Read.* import ..." works
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from Read.ReadSpenvis_tri import readSpenvis_tri

t = sp.Symbol('t')  # time in hours

T_c = 1.23 * 24  # Decay time in hours

f0 = 4.511444e10   # cm-2 s-1

flux = f0 * sp.exp(-t / T_c)

print(flux)

T = sp.Symbol('T')  # time in hours

TotalFluence = sp.integrate(flux, (t, 0, sp.oo)) * 3600  # cm-2. The factor 3600 converts from per hour to per second
# Print total fluence rounded to 4 significant figures
print(f"TotalFluence: {TotalFluence:.4g} cm-2")

Fluence = sp.integrate(flux, (t, 0, T)) * 3600  # cm-2. The factor 3600 converts from per hour to per second
print(Fluence)

flux = sp.lambdify(t, flux, "numpy")
Fluence = sp.lambdify(T, Fluence, "numpy")

Evals = np.linspace(0, 100, num=101) # 100 hours in 1 hour steps

#for t in Evals:
#    print(f"{t:.4}", f"{flux(t):.4}")

for T in Evals:
    print(f"{T:.4}", f"{Fluence(T):.4}")


fig, ax1 = plt.subplots(1, figsize=(5, 3))

plt.plot(Evals/24, flux(Evals), label="Electron Flux", color="C0")
ax1.set_ylabel("Total Electron Flux [cm-2 s-1]", color="C0")
ax1.legend(loc='upper center', bbox_to_anchor=(0.3, 1.0))
ax1.tick_params(axis='y', colors='C0')
plt.grid(which='both')

ax2 = ax1.twinx()
plt.plot(Evals/24, Fluence(Evals), label="Electron Fluence", color="C1")
ax2.set_ylabel("Total Electron Fluence [cm-2]", color="C1")
plt.legend(loc='lower center', bbox_to_anchor=(0.3, 0.0))
ax2.tick_params(axis='y', colors='C1')

#plt.grid(which='both')
#plt.yscale("log")
#plt.xscale("log")
plt.title("Time dependent electron flux and event fluence")
ax1.set_xlabel("Time [days]")

plt.savefig(r"\\data.triton.aalto.fi\work\fetzera1\Spectra\Carrington\Electron-EVT-2026\Electron-EVT-2026TimeDependence.pdf", format='pdf', bbox_inches="tight")
#plt.show()


""" 
ThickList = [1, 2, 4, 8, 16]
DosePerHour = [5.78169891, 1.3798122553287573, 0.3627583749127644, 0.14702413643842302, 0.07671414538190979]
Error = [0.7639131599999995, 0.17434836835347167, 0.024627762742141568, 0.007887851333699764, 0.006097929079940834]

for t, thick in enumerate(ThickList):

    #plt.plot(Evals/24, flux(Evals)*DosePerHour, label=str(thick) + " mm", marker="x")
    plt.errorbar(Evals/24, Fluence(Evals)*DosePerHour[t], Fluence(Evals)*Error[t], fmt=' ', capsize=5, label=str(thick) + " mm")

####### Plot 10kRad line #########
CriticalDose = [10 for i in Evals]
plt.plot(Evals/24, CriticalDose, color='k', linewidth=2, label='10 krad')
CriticalDose = [100 for i in Evals]
plt.plot(Evals/24, CriticalDose, '--', color='k', linewidth=2, label='100 krad')

#plt.ylim(5e-2, 1e3)
# plt.xlim(0, 4.1)
plt.grid()
plt.yscale("log")
#plt.xscale("log")
plt.title("Ionising dose due to Carrington event electron flux\n behind aluminium shielding of various thicknesses")
plt.xlabel("Time since start of Carrington event [days]")
plt.ylabel("Ionising dose due to electrons [krad]")
plt.legend() 

# plt.savefig("/u/02/fetzera1/unix/Desktop/TritonPlots/Carrington/CarringtonTimeDep3D.pdf", format='pdf', bbox_inches="tight")
plt.show()
"""