import matplotlib.pyplot as plt
from ReadG4root import readG4root
import math

MeVDat = readG4root("/scratch/work/fetzera1/TEST/SiChip in 1cm cube for Spenvis Test/root/electrons500kevpower.root")  # 0 = Primarykine; 1 = Dose

TotalMeV = sum(MeVDat[1])

print("Total Dose in MeV:", TotalMeV)

# ----------- Norm factors from GPS spectrum ------------
NORM_FACTOR_SPECTRUM = 5.886798E+14  # *********************************************************************************
NORM_FACTOR_ANGULAR = 2.500000E-01
Norm = NORM_FACTOR_SPECTRUM * NORM_FACTOR_ANGULAR
# -----------------------------------------

# --------- Simulated Particle Fluence --------
Npart = 2e9  # *********************************************************************************************************
Radius = 1  # cm *******************************************************************************************************
Area = 4 * math.pi * Radius * Radius  # cm2
SimulatedFluence = Npart / Area
# ------------------------------------

# --------- Mission MeV --------
MissionMeV = TotalMeV / SimulatedFluence * Norm
# -------------------------------

print("MissionMeV:", MissionMeV)  # MeV

# ------ Sensitive Volume ------
SiDensity = 2.33 * 1e-3  # kg/cm3
SiThick = 0.05  # cm
length = 1  # cm *******************************************************************************************************
SiVol = SiThick * length * length
print("SiVol", SiVol)  # 0.05 cm3
SiMass = SiVol * SiDensity
print("SiMass", SiMass)  # kg
# -----------------------------

JperMev = 1.602E-13
MissionJoule = MissionMeV * JperMev
print("MissionJoule", MissionJoule)  # J

Grays = MissionJoule / SiMass
print("Grays", Grays)  # 37763 J/kg

kRads = Grays/10
print("kRads", kRads)  # 3776 krad

# rads reported by GRAS log file : 3.77678e+06
