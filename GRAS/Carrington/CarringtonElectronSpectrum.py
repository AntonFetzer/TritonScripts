import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from GRAS.Read.ReadSpenvis_tri import readSpenvis_tri

f0 = 10 ** 10
E0 = 0.13
a = 3.7

E = sp.Symbol('E')

f = 4 * sp.pi * f0 * (E / E0) ** (-a)

fdiff = sp.diff(f, E)
print(fdiff)

f = sp.lambdify(E, f, "numpy")
fdiff = sp.lambdify(E, fdiff, "numpy")

Evals = np.geomspace(0.13, 10, num=10)

print("Differential")
for E in Evals:
    print(f"{E:.4}", f"{-fdiff(E):.4}")
    #print("/gps/hist/point", f"{E:.4}", f"{-fdiff(E):.4}")

print("Integral")
for E in Evals:
    print(f"{E:.4}", f"{f(E):.4}")

file = "/l/triton_work/Spectra/A9/spenvis_tri.txt"

Protons, Electrons = readSpenvis_tri(file)

plt.figure(1)
#plt.plot(Protons[0], Protons[1], label="Proton Integral")
#plt.plot(Protons[0], Protons[2], label="Proton Differential")
plt.plot(Electrons[0], Electrons[1], '-.', label="AE9 Integral")
plt.plot(Electrons[0], Electrons[2], '-.', label="AE9 Differential")

plt.plot(Evals, f(Evals), label="Carrington Integral")
plt.plot(Evals, -fdiff(Evals), label="Carrington Differential")
plt.legend()
plt.grid(which='both')
plt.yscale("log")
plt.xscale("log")
plt.title("Spectra Comparison")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Flux [cm-2 s-1] or [cm-2 s-1 MeV-1]")

#plt.savefig("/l/triton_work/CARRINGTON/SpectrumComparison.eps", format='eps', bbox_inches="tight")

plt.show()
