def MeVtokRad_2D(MeV, NORM_FACTOR_SPECTRUM, Npart):
    # ----------- Norm factors from GPS spectrum ------------
    # NORM_FACTOR_SPECTRUM = 8.003046E+14  # ***************************************************************************
    NORM_FACTOR_ANGULAR = 2.500000E-01
    Norm = NORM_FACTOR_SPECTRUM * NORM_FACTOR_ANGULAR
    # -----------------------------------------

    # --------- Simulated Particle Fluence --------
    # Npart = 2e9  # ***************************************************************************************************
    # ------------------------------------

    # --------- Mission MeV --------
    MissionMeV = MeV / Npart * Norm
    # -------------------------------

    print("MissionMeV:", MissionMeV)  # MeV

    # ------ Sensitive Volume ------
    SiDensity = 2.33 * 1e-3  # kg/cm3
    SiThick = 0.05  # cm
    MassPerArea = SiDensity * SiThick  # kg/cm2
    print("SiMass", MassPerArea)  # kg
    # -----------------------------

    JperMev = 1.602E-13
    MissionJoule = MissionMeV * JperMev
    print("MissionJoule", MissionJoule)  # J

    Grays = MissionJoule / MassPerArea
    print("Grays", Grays)  # J/kg

    kRads = Grays / 10
    print("kRads", kRads)  # krad

    return kRads
