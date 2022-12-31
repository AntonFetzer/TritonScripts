import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from GRAS.Read.ReadSpenvis_tri import readSpenvis_tri
from scipy.optimize import curve_fit
from GRAS.Dependencies.TotalGRASHistos import totalGRASHistos

def DiffSpec(E, Int):
    # f = Fi - Fi-1 / dE
    DiffF = []
    for i in range(1, len(Int)):
        DiffF.append((Int[i-1] - Int[i]) / (E[i] - E[i-1]))
    DiffF.append(0)
    return(DiffF)


def PowerLawFunc(x, a, b, c):
    return a * x**(-b) + c

AP9Ener = [1.00E-01, 2.00E-01, 4.00E-01, 6.00E-01, 8.00E-01, 1.00E+00, 2.00E+00, 4.00E+00, 6.00E+00, 8.00E+00, 1.00E+01, 1.50E+01, 2.00E+01, 3.00E+01, 5.00E+01, 6.00E+01, 8.00E+01, 1.00E+02, 1.50E+02, 2.00E+02, 3.00E+02, 4.00E+02, 7.00E+02, 1.20E+03, 2.00E+03]
AP9Int = [5.21E+07, 3.40E+07, 1.61E+07, 8.31E+06, 4.47E+06, 2.77E+06, 7.53E+05, 1.50E+05, 4.79E+04, 2.20E+04, 1.34E+04, 6.33E+03, 4.03E+03, 2.61E+03, 1.58E+03, 1.33E+03, 9.99E+02, 7.68E+02, 4.18E+02, 2.35E+02, 7.98E+01, 2.82E+01, 3.42E+00, 1.40E-01, 0.00E+00]
AP9Diff = [2.11E+08, 1.51E+08, 6.44E+07, 2.91E+07, 1.38E+07, 7.40E+06, 1.45E+06, 1.76E+05, 3.21E+04, 8.61E+03, 3.48E+03, 9.39E+02, 3.54E+02, 1.12E+02, 3.41E+01, 2.24E+01, 1.40E+01, 1.03E+01, 5.33E+00, 2.96E+00, 1.03E+00, 4.08E-01, 5.41E-02, 4.10E-03, 0.00E+00]

#plt.plot(AP9Ener, AP9Int, label="AP9 Integral")
#plt.plot(AP9Ener, AP9Diff, label="AP9 Differential")

Ener = [1, 3, 5, 7, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100]
IntF = [5.00E+08, 1.00E+08, 1.00E+08, 2.50E+07, 2.50E+07, 1.00E+07, 1.00E+07, 3.00E+06, 3.00E+06, 3.00E+06, 1.00E+06, 1.00E+06, 3.00E+05, 3.00E+05, 3.00E+05, 1.00E+05, 1.00E+05, 1.00E+05, 1.00E+05]

NewEner = [1, 4, 8.5, 17.5, 30, 42.5, 55, 85]
NewIntF = [500000000, 100000000, 25000000, 10000000, 3000000, 1000000, 300000, 100000]

#plt.fill_between(Ener, IntF, y2=0, label="Original Spectrum")
plt.plot(Ener, IntF, label="Original Spectrum")
plt.plot(NewEner, NewIntF, label="Simplified Spectrum")
plt.title("Solar Energetic Particle Integral Fluence Spectrum")
plt.xlabel("Kinetic energy [MeV]")
plt.ylabel("Fluene [cm-2]")
plt.yscale("log")
plt.xscale("log")
plt.grid(which='both')

path = "/home/anton/Desktop/triton_work/CARRINGTON/HistogramSEP/Res/"
DoseHist, PrimaryHist = totalGRASHistos(path, "Prot")
lowerID = 0
upperID = 1
meanID = 2
valueID = 3
errorID = 4
entriesID = 5
plt.bar(PrimaryHist[:, lowerID], PrimaryHist[:, entriesID]*5e3, width=PrimaryHist[:, upperID] - PrimaryHist[:, lowerID],
        align='edge', label="Geant4 Histogram")

#######################   Differential Spectrum ########################

AE9DiffTEST = DiffSpec(AP9Ener, AP9Int)
#plt.plot(AP9Ener, AE9DiffTEST, label="AE9 Difftest")

DiffF = DiffSpec(NewEner, NewIntF)
#plt.plot(NewEner, DiffF, label="Differential Fluence")
plt.legend()
plt.savefig("/home/anton/Desktop/triton_work/Spectra/Carrington/SEPIntegralSpectrumPythonFix.eps", format='eps', bbox_inches="tight")
#plt.show()

'''
FitRes = curve_fit(PowerLawFunc, Ener, np.log(IntF), p0=np.asarray([1e8, 1, 1]), bounds=(0, np.inf))

print(FitRes[0])

Fit = PowerLawFunc(Ener, FitRes[0][0], FitRes[0][1], FitRes[0][2])
Fit = np.exp(Fit)

for f in Fit:
    print(f)

plt.plot(Ener, Fit, label="Fit")

plt.savefig("/home/anton/Desktop/triton_work/Spectra/Carrington/SEPIntegralSpectrumPython.eps", format='eps', bbox_inches="tight")
plt.show()
'''