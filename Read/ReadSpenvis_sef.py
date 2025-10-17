import numpy as np
import matplotlib.pyplot as plt

def readSpenvis_sef_protons(fileName):
    """ 
    Reads the SAPPHIRE solar proton fluence spectrumn from spenvis_sef.txt files.

    Args:
        file (str): Path to the spenvis_sef.txt file.
    
    Returns:
        dict: Dictionary with numpy arrays of fluence spectrum:
            - 'Energy': Energy in MeV
            - 'IFluence': Integral Fluence in cm-2
            - 'DFluence': Differential Fluence in cm-2 MeV-1
    """
    keys = ['Energy', 'IFluence', 'DFluence']
    ProtonTable = {key: [] for key in keys}
    
    ReadFlag = 0
    with open(fileName, 'r') as f:
        for line in f:
            if ReadFlag == 0:
                if "'Proton Exposure Time'" in line:
                    ReadFlag = 1  # Found the proton table
            elif ReadFlag == 1:
                if "'End of Block'" in line:
                    break  # Done reading relevant block
                else:
                    # Parse values
                    values = [value.strip() for value in line.split(',')]
                    ProtonTable['Energy'].append(float(values[0]))
                    ProtonTable['IFluence'].append(float(values[1]))
                    ProtonTable['DFluence'].append(float(values[2]))
    
    # Convert lists inside the dicts to numpy arrays
    for key in keys:
        ProtonTable[key] = np.array(ProtonTable[key])

    return ProtonTable




AtomicMass = [1.008, 4.003, 6.941, 9.012, 10.811, 12.011, 14.007, 15.999, 18.998, 20.18, 22.99, 24.305, 26.982,
              28.086, 30.974, 32.065, 35.453, 39.948, 39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938,
              55.845, 58.933, 58.693, 63.546, 65.39, 69.723, 72.64, 74.922, 78.96, 79.904, 83.8, 85.468, 87.62,
              88.906, 91.224, 92.906, 95.94, 98, 101.07, 102.906, 106.42, 107.868, 112.411, 114.818, 118.71, 121.76,
              127.6, 126.905, 131.293, 132.906, 137.327, 138.906, 140.116, 140.908, 144.24, 145, 150.36, 151.964,
              157.25, 158.925, 162.5, 164.93, 167.259, 168.934, 173.04, 174.967, 178.49, 180.948, 183.84, 186.207,
              190.23, 192.217, 195.078, 196.967, 200.59, 204.383, 207.2, 208.98, 209, 210, 222, 223, 226, 227,
              232.038, 231.036, 238.029]

############ SAPPHIRE Spectra are in Fluence with minimum duration 6 months !!! ###############################

def readSpenvis_sef_Ions(fileName):
    print("Reading in", fileName)
    f = open(fileName, "r")

    # Ignore the Proton Table, because the same data is also in the Ion table.
    # The ion table starts at Hydrogen !!!
    Iontable = []
    readflag = 0

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

    NumSpecies = 92

    IonData = np.zeros((rows, cols-1), dtype=np.float64)
    EnergyPerNucleon = np.zeros(rows, dtype=np.float64)

    for i, line in enumerate(Iontable):
        EnergyPerNucleon[i] = line[0]
        IonData[i] = line[1:cols]     # Reading the Fluence
        #IonData[i] /= (1.578e+7)   # The sef.txt contains the Fluence per 6 months not the flux!

    # print(Energy)

    Data = np.zeros((NumSpecies, len(Iontable), 3), dtype=float)

    for i in range(NumSpecies):
        Energy = EnergyPerNucleon*AtomicMass[i]
        #print("Species:", i)
        #print("AtomicMass:", AtomicMass[i])
        #print("Energy:", Energy)
        #print("IonData:", IonData[:, i+1])
        Data[i, :, 0] = Energy
        Data[i, :, 1] = IonData[:, i]  # Integral Flux
        Data[i, :, 2] = IonData[:, NumSpecies+i]  # Differential Flux

    return Data
# Data[ Z-NUmber , DatapointNum, Energy or Integral Flux or Differential Flux ]


if __name__ == "__main__":
    ############ SAPPHIRE Spectra are in Fluence with minimum duration 6 months !!! ###############################
    ############ Check duration and normalisation of the spectrum #################################################

    Data = readSpenvis_sef_protons("/l/triton_work/Spectra/Carrington/GEO/spenvis_sef.txt")

    print(Data['Energy'])
    print(Data['IFluence'])
    print(Data['DFluence'])

    '''
    DataT = readSpenvis_sef("/l/triton_work/Spectra/ISS/spenvis_sef.txt")

    IntorDiff = 1  # 1 for Int 2 for Diff

    #plt.bar(range(93), DataT[:, 0, 1])
    #print(DataT[:, 0, 1])

    Species = ['H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl',
               'Ar', 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se',
               'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb',
               'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
               'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
               'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U ']

    #print(np.shape(DataT))

    for Specie in range(np.shape(DataT)[0]):
        FluxMaxE = np.max(DataT[Specie, :, IntorDiff] * DataT[Specie, :, 0])
        if FluxMaxE > 1e-10:
            print(f"Species: {Species[Specie]}, Max Flux * Energy: {FluxMaxE}")
            plt.plot(DataT[Specie, :, 0], DataT[Specie, :, IntorDiff], label=Species[Specie])
            #for Energy, Flux in zip(DataT[Specie, :, 0], DataT[Specie, :, IntorDiff]):
                #if Energy > 10:
                    #print(f"{Energy} {Flux}")

    #plt.xlim(10, 1e5)
    #plt.ylim(1e-5, 1e6)
    plt.xscale("log")
    plt.yscale("log")

    if IntorDiff == 1:
        plt.title("Solar Ion Integral Flux")
        plt.ylabel("Integral Flux [cm-2 s-1]")
    elif IntorDiff == 2:
        plt.title("Solar Ion Differential Flux")
        plt.ylabel("Differential Flux [cm-2 s-1 MeV-1]")

    plt.xlabel("Energy [MeV]")
    plt.legend()
    plt.grid(which='both')
    plt.show()

'''