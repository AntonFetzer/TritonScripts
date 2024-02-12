import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from sympy import *
from GRAS.Read.ReadSpenvis_tri import readSpenvis_tri

t = sp.Symbol('t')  # time in hours

T_c = 1.23 * 24  # Decay time in hours

flux = exp(-t / T_c)

print(flux)

T = sp.Symbol('T')  # time in hours

TotalFluence = sp.integrate(flux, (t, 0, oo))
print(TotalFluence)

Fluence = sp.integrate(flux, (t, 0, T))

print(Fluence)

flux = sp.lambdify(t, flux, "numpy")
Fluence = sp.lambdify(T, Fluence, "numpy")

Evals = np.linspace(1, 100, num=100) # 100 hours in 1 hour steps

#for t in Evals:
#    print(f"{t:.4}", f"{flux(t):.4}")

for T in Evals:
    print(f"{T:.4}", f"{Fluence(T):.4}")

'''
fig, ax1 = plt.subplots()

plt.plot(Evals/24, flux(Evals), label="Relative flux f/f_0", color="C0", marker="x")
ax1.set_ylabel("Relative Flux f/f_0", color="C0")
ax1.legend(loc='upper center')
ax1.tick_params(axis='y', colors='C0')

ax2 = ax1.twinx()
plt.plot(Evals/24, Fluence(Evals), label="Fluence in hours of f_0", color="C1", marker="x")
ax2.set_ylabel("Fluence in hours of f_0", color="C1")
ax2.legend(loc='lower center')
ax2.tick_params(axis='y', colors='C1')

plt.grid(which='both')
#plt.yscale("log")
#plt.xscale("log")
plt.title("Time dependent flux and fluence")
ax1.set_xlabel("Time [days]")
'''



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
