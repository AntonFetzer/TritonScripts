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

AP9Energy = np.array([1.00E+01, 1.50E+01, 2.00E+01, 3.00E+01, 5.00E+01, 6.00E+01, 8.00E+01, 1.00E+02, 1.50E+02, 2.00E+02, 3.00E+02, 4.00E+02, 7.00E+02, 1.20E+03])
AP9IntegralGTO = np.array([1.34E+04, 6.33E+03, 4.03E+03, 2.61E+03, 1.58E+03, 1.33E+03, 9.99E+02, 7.68E+02, 4.18E+02, 2.35E+02, 7.98E+01, 2.82E+01, 3.42E+00, 1.40E-01])
AP9IntegralLEO = np.array([3.92E+01, 3.04E+01, 2.58E+01, 2.09E+01, 1.52E+01, 1.34E+01, 1.06E+01, 8.63E+00, 5.29E+00, 3.30E+00, 1.28E+00, 4.97E-01, 6.44E-02, 4.57E-04])

SAPPHIREEnergy = np.array([10.00, 11.00, 12.00, 14.00, 16.00, 18.00, 20.00, 22.00, 25.00, 28.00, 32.00, 35.00, 40.00, 45.00, 50.00, 55.00, 63.00, 71.00, 80.00, 90.00, 100.00, 110.00, 120.00, 140.00, 160.00, 180.00, 200.00])
SAPPHIRELEOIntegral = np.array([1.460E+06, 1.352E+06, 1.260E+06, 1.111E+06, 9.973E+05, 9.070E+05, 8.353E+05, 7.789E+05, 7.141E+05, 6.657E+05, 6.178E+05, 5.906E+05, 5.563E+05, 5.326E+05, 5.159E+05, 5.039E+05, 4.912E+05, 4.839E+05, 4.791E+05, 4.758E+05, 4.737E+05, 4.723E+05, 4.712E+05, 4.699E+05, 4.663E+05, 4.540E+05, 4.317E+05])


C = 3600*24 # Converts from Flux to daily Fluence

#for i, e in enumerate(E):
#    print(e, SPE2003[i])

Colours = ['C1', 'C0', 'C2', 'C7']



plt.figure(1)

plt.fill_between(E, Plus*C, ExpectedInt*C, color='C1', alpha=0.5)
plt.fill_between(E, ExpectedInt*C, Minus*C, color='C2', alpha=0.5)

plt.plot(E, Plus*C, label="Carrington SEP +2 Sigma", color='C1')

#plt.fill_between(E, ExpectedInt*C, color='C0', alpha=0.5)
plt.plot(E, ExpectedInt*C, label="Carrington SEP EVT", color='C0', linewidth=4)

#plt.fill_between(E, Minus*C, color='C2', alpha=0.5)
plt.plot(E, Minus*C, label="Carrington SEP -2 Sigma", color='C2')

plt.plot(E, TownsendCarringotn, ':', label="Townsend Carrington Estimate", color='C9', linewidth=4)
plt.plot(E, SPE2003, '--', label="SPE 2003", color='C8', linewidth=4)
plt.plot(AP9Energy, AP9IntegralGTO*C, label="Trapped Protons GTO", color='C3')
plt.plot(AP9Energy, AP9IntegralLEO*C, label="Trapped Protons LEO", color='C7')

#plt.plot(SAPPHIREEnergy, SAPPHIRELEOIntegral/(6*30), label="LEO SEP")

plt.xlim(9, 250)
plt.ylim(1e6, 2e11)
plt.legend(loc='lower left')
plt.grid(which='both')
plt.yscale("log")
plt.xscale("log")
plt.title("Spectra Comparison")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Fluence per Day [cm-2]")

plt.savefig("/l/TritonPlots/Carrington/SpectrumComparison.svg", format='svg', bbox_inches="tight")

#plt.show()