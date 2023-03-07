import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from GRAS.Read.ReadSpenvis_tri import readSpenvis_tri

E = [10, 30, 60, 100, 200]

ExpectedInt = np.array([5.781E+05, 1.068E+05, 2.890E+04, 1.382E+04, 3.644E+03])
Minus = np.array([2.765E+05, 7.037E+04, 1.759E+04, 5.404E+03, 1.382E+03])
Plus = np.array([1.634E+06, 2.011E+05, 6.786E+04, 5.529E+04, 1.759E+04])
SPE2003 = np.array([1.6e10, 3.5e9, 7e8, 1.5e8, 2e7])

TownsendCarringotn = np.array([1e11, 2e10, 4e9, 9e8, 6e7])


C = 3600*24 # Converts from Flux to daily Fluence

#for i, e in enumerate(E):
#    print(e, SPE2003[i])

Colours = ['C1', 'C0', 'C2', 'C7']


file = "/home/anton/Desktop/triton_work/Spectra/A9/spenvis_tri.txt"

Protons, Electrons = readSpenvis_tri(file)

plt.figure(1)
plt.plot(Protons[0], Protons[1]*C, label="AP9 GTO", color='C3')
#plt.plot(Protons[0], Protons[2], label="Proton Differential")


plt.plot(E, Plus*C, label="+2 Sigma", color='C1')
plt.plot(E, ExpectedInt*C, label="Expected Integral", color='C0')
plt.plot(E, Minus*C, label="-2 SIgma", color='C2')
plt.plot(E, SPE2003, label="SPE 2003", color='C8')
plt.plot(E, TownsendCarringotn, label="Townsend Carrington Estimate", color='C9')


plt.xlim(9, 300)
plt.ylim(1e7, 1e12)
plt.legend()
plt.grid(which='both')
plt.yscale("log")
plt.xscale("log")
plt.title("Spectra Comparison")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Fluence per Day [cm-2]")

plt.savefig("/home/anton/Desktop/TritonPlots/Carrington/SpectrumComparison.pdf", format='pdf', bbox_inches="tight")

#plt.show()