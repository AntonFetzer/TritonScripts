import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from GRAS.Read.ReadSpenvis_tri import readSpenvis_tri
from PowerLawFunction import PowerLawFunction

## Input Parameters

F0_low = 10 ** (4.3 - 0.7)  # cm-2 s-1 str-1
F0_mean = 10 ** 4.3  # cm-2 s-1 str-1
F0_high = 10 ** (4.3 + 0.7)  # cm-2 s-1 str-1

E0 = 10

a_low = 0.11 + 0.1
a_mean = 0.11
a_high = 0.11 - 0.1

## Converting Units

F0_low = F0_low * 4 * sp.pi  # cm-2 s-1
F0_mean = F0_mean * 4 * sp.pi  # cm-2 s-1
F0_high = F0_high * 4 * sp.pi  # cm-2 s-1

## Getting the Integral Functions

F_low = PowerLawFunction(F0_low, a_low, E0)
F_mean = PowerLawFunction(F0_mean, a_mean, E0)
F_high = PowerLawFunction(F0_high, a_high, E0)

## Gettting the differential Funticns

f_low = PowerLawFunction(F0_low, a_low, E0, "Diff")
f_mean = PowerLawFunction(F0_mean, a_mean, E0, "Diff")
f_high = PowerLawFunction(F0_high, a_high, E0, "Diff")

Evals = np.geomspace(10, 200, num=100)

print("Integral")
for E in Evals:
    print(f"{E:.4}", f"{F_mean(E):.4}")

#print("Differential")
#for E in Evals:
#    print(f"{E:.4}", f"{-f_high(E):.4}")


'''
## Import AP9 AE9 Data

file = "/home/anton/Desktop/triton_work/Spectra/A9/spenvis_tri.txt"

Protons, Electrons = readSpenvis_tri(file)

## Integral Plot

plt.figure(0)
plt.plot(Protons[0], Protons[1], label="AP9 Integral")
plt.plot(Electrons[0], Electrons[1], '-.', label="AE9 Integral")

plt.plot(Evals, F_high(Evals), label="Carrington Integral High")
plt.plot(Evals, F_mean(Evals), label="Carrington Integral Mean")
plt.plot(Evals, F_low(Evals), label="Carrington Integral Low")
plt.legend()
plt.grid(which='both')
plt.yscale("log")
plt.xscale("log")
plt.title("Integral Spectra Comparison")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Integral Flux [cm-2 s-1]")

plt.savefig("/home/anton/Desktop/triton_work/CARRINGTON/IntegralSpectrumComparison.pdf", format='pdf', bbox_inches="tight")

## Differential Plot

plt.figure(1)
plt.plot(Protons[0], Protons[2], label="AP9 Differential")
plt.plot(Electrons[0], Electrons[2], '-.', label="AE9 Differential")

plt.plot(Evals, -f_low(Evals), label="Carrington Differential Low")
plt.plot(Evals, -f_mean(Evals), label="Carrington Differential Mean")
plt.plot(Evals, -f_high(Evals), label="Carrington Differential High")
plt.legend()
plt.grid(which='both')
plt.yscale("log")
plt.xscale("log")
plt.title("Differential Spectra Comparison")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Differential flux [cm-2 s-1 MeV-1]")

plt.savefig("/home/anton/Desktop/triton_work/CARRINGTON/DifferentialSpectrumComparison.pdf", format='pdf', bbox_inches="tight")
#plt.show()
'''