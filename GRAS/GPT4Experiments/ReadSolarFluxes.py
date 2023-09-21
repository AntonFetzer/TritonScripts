import numpy as np
import matplotlib.pyplot as plt


def readSpenvis_sef(fileName):
    """
    Function to read SEF data file from SPENVIS.

    Args:
    fileName (str): The path of the SEF file.

    Returns:
    numpy.ndarray: A 3D numpy array containing data for each ion species.
    """
    print("Reading in", fileName)

    FLUENCE_TO_FLUX_CONVERSION = 1.578e+7 # The sef.txt contains the Fluence per 6 months not the flux!
    NUM_SPECIES = 92

    with open(fileName, "r") as f:
        readflag = 0
        Iontable = []
        for line in f:
            if readflag == 0:
                if "Differential Fluence of ','SPECIES" in line:
                    readflag = 1
            elif readflag == 1:
                if "End of File" in line:
                    readflag = 2
                else:
                    Iontable.append(np.fromstring(line, dtype=np.float64, sep=','))


    rows, cols = np.shape(Iontable)

    IonData = np.zeros((rows, cols - 1), dtype=np.float64)
    EnergyPerNucleon = np.zeros(rows, dtype=np.float64)

    for i, line in enumerate(Iontable):
        EnergyPerNucleon[i] = line[0]
        IonData[i] = line[1:cols]  # Reading the Fluence
        IonData[i] /= FLUENCE_TO_FLUX_CONVERSION

    AtomicMass = [1.008, 4.003, 6.941, 9.012, 10.811, 12.011, 14.007, 15.999, 18.998, 20.18, 22.99, 24.305, 26.982,
                  28.086, 30.974, 32.065, 35.453, 39.948, 39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938,
                  55.845, 58.933, 58.693, 63.546, 65.39, 69.723, 72.64, 74.922, 78.96, 79.904, 83.8, 85.468, 87.62,
                  88.906, 91.224, 92.906, 95.94, 98, 101.07, 102.906, 106.42, 107.868, 112.411, 114.818, 118.71, 121.76,
                  127.6, 126.905, 131.293, 132.906, 137.327, 138.906, 140.116, 140.908, 144.24, 145, 150.36, 151.964,
                  157.25, 158.925, 162.5, 164.93, 167.259, 168.934, 173.04, 174.967, 178.49, 180.948, 183.84, 186.207,
                  190.23, 192.217, 195.078, 196.967, 200.59, 204.383, 207.2, 208.98, 209, 210, 222, 223, 226, 227,
                  232.038, 231.036, 238.029]

    Data = np.zeros((NUM_SPECIES, len(Iontable), 3), dtype=float)

    for i in range(NUM_SPECIES):
        Energy = EnergyPerNucleon * AtomicMass[i]
        Data[i, :, 0] = Energy
        Data[i, :, 1] = IonData[:, i]  # Integral Flux
        Data[i, :, 2] = IonData[:, NUM_SPECIES + i]  # Differential Flux

    return Data



if __name__ == "__main__":
    ############ SAPPHIRE Spectra are in Fluence with minimum duration 6 months !!! ###############################
    ############ Check duration and normalisation of the spectrum #################################################
    DataT = readSpenvis_sef("/l/triton_work/Spectra/SAPPHIRE-GTO/spenvis_sef.txt")

    print(DataT)