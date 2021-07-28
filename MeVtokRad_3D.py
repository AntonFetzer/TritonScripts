import math


def MeVtokRad_3D(MeV, NORM_FACTOR_SPECTRUM, Npart, Radius):  # Radius in cm

    # ----------- Norm factors from GPS spectrum ------------
    NORM_FACTOR_ANGULAR = 2.500000E-01
    Norm = NORM_FACTOR_SPECTRUM * NORM_FACTOR_ANGULAR
    # -----------------------------------------

    # --------- Simulated Particle Fluence --------
    Area = 4 * math.pi * Radius * Radius  # cm2
    SimulatedFluence = Npart / Area
    # ------------------------------------

    # --------- Mission MeV --------
    MissionMeV = MeV / SimulatedFluence * Norm
    # -------------------------------

    # print("MissionMeV:", MissionMeV)  # MeV

    # ------ Sensitive Volume ------
    SiDensity = 2.33 * 1e-3  # kg/cm3
    SiThick = 0.05  # cm
    SiLength = 0.3  # cm ***********************************************************************************************
    SiWidth = 0.3  # cm ************************************************************************************************
    SiVol = SiThick * SiLength * SiWidth
    # print("SiVol", SiVol)  # cm3
    SiMass = SiVol * SiDensity
    # print("SiMass", SiMass)  # kg
    # -----------------------------

    JperMev = 1.602E-13
    MissionJoule = MissionMeV * JperMev
    # print("MissionJoule", MissionJoule)  # J

    Grays = MissionJoule / SiMass
    # print("Grays", Grays)  # J/kg

    kRads = Grays / 10
    # print("kRads", kRads)  # krad

    return kRads

