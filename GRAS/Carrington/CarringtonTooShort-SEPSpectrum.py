import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from GRAS.Read.ReadSpenvis_tri import readSpenvis_tri

J0 = 10 ** 11 / (4 * 3600)  # SEP Millenium event differential fluence [cm-2 MeV-1] ~ per 4 hours?
#J0 = 5* 10 ** 10  # SEP Century event differential fluence [cm-2 MeV-1] ~ per 4 hours?
E0 = 10  # MeV
a = 0.9  # Lol I don't think so

E = sp.Symbol('E')

J = J0 * (E / E0) ** (-a)

Jint = sp.integrate(J, E)
print(Jint)

J = sp.lambdify(E, J, "numpy")
Jint = sp.lambdify(E, Jint, "numpy")

Evals = np.geomspace(10, 50, num=2)

print("Differential:")
for E in Evals:
    print(f"{E:.4}", f"{J(E):.4}")

print("Integral:")
for E in Evals:
    print(f"{E:.4}", f"{Jint(E):.4}")

file = "/l/triton_work/Spectra/A9-GTO/spenvis_tri.txt"

Protons, Electrons = readSpenvis_tri(file)

plt.figure(1)
#plt.plot(Protons[0], Protons[1], label="AP9 Integral flux")
plt.plot(Protons[0], Protons[2], label="Ap9 Differential flux")
#plt.plot(Protons[0], Protons[1]*2000, label="AP9 Integral flux x2000")
plt.plot(Protons[0], Protons[2]*5e4, label="Ap9 Differential flux x2000")
#plt.plot(Electrons[0], Electrons[1], '-.', label="AE9 Integral")
#plt.plot(Electrons[0], Electrons[2], '-.', label="AE9 Differential")

plt.plot(Evals, J(Evals), label="SEP Millenium event differential flux")
#plt.plot(Evals, Jint(Evals), label="SEP Millenium event integral flux")
plt.legend()
plt.grid(which='both')
plt.yscale("log")
plt.xscale("log")
plt.title("Spectra Comparison")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Fluence [cm-2 s-1] or [cm-2 s-1 MeV-1]")

#plt.savefig("/l/triton_work/CARRINGTON/SpectrumComparisonSEP.pdf", format='pdf', bbox_inches="tight")

plt.show()
