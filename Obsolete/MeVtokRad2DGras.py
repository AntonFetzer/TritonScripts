def MeVtokRad_2DGRAS(MeV, NORM_FACTOR_SPECTRUM):

    # ----------- Norm factors from GPS spectrum ------------
    NORM_FACTOR_ANGULAR = 2.500000E-01
    Norm = NORM_FACTOR_SPECTRUM * NORM_FACTOR_ANGULAR
    # -----------------------------------------

    # --------- Mission MeV --------
    MissionMeV = MeV * Norm
    # -------------------------------

    # print("MissionMeV:", MissionMeV)  # MeV

    # ------ Sensitive Volume ------
    SiDensity = 2.33 * 1e-3  # kg/cm3
    SiThick = 0.05  # cm
    MassPerArea = SiDensity * SiThick  # kg/cm2
    # print("SiMass", MassPerArea)  # kg
    # -----------------------------

    JperMev = 1.602E-13
    MissionJoule = MissionMeV * JperMev
    # print("MissionJoule", MissionJoule)  # J

    Grays = MissionJoule / MassPerArea
    # print("Grays", Grays)  # J/kg

    kRads = Grays / 10
    # print("kRads", kRads)  # krad

    return kRads


if __name__ == "__main__":
    x = MeVtokRad_2DGRAS(1, 2.432839E+07)
    #Scale = 30 * 24 * 60 * 60
    print(x)
