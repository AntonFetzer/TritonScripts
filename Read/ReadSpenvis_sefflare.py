import numpy as np
import matplotlib.pyplot as plt


def readSpenvis_sefflare(file):
    """
    Reads the proton peak flux from Solar Flare spenvis_sefflare.txt files.

    The spenvis_sefflare.txt file contains the proton peak flux data for a solar flare event.
    Energy in                   MeV
    Integral Peak Flux in       m-2 sr-1 s-1
    Differential Peak Flux in   m-2 sr-1 s-1 MeV-1

    Args:
        file (str): The path to the file to be read.

    Returns:
        dict: A dictionary containing the proton peak flux data with keys 'Energy', 'IFlux', and 'DFlux'.
    """
    print("Reading in", file)

    keys = ['Energy', 'IFlux', 'DFlux']

    # Generating the dictionary from the list of keys
    Sefflare = {key: [] for key in keys}

    ReadFlag = 0

    # Open the file and read line by line
    with open(file, 'r') as f:
        for line in f:
            if ReadFlag == 0:
                if "Duration" in line:
                    continue
                elif "Solar particle model" in line:
                    continue
                elif "Proton Exposure Time" in line:
                    ReadFlag = 1
            elif ReadFlag == 1:
                if "End of Block" in line:
                    break
                else:
                    # Split the line into parts based on the comma
                    values = line.split(',')
                    Sefflare['Energy'].append(float(values[0]))
                    Sefflare['IFlux'].append(float(values[1]))
                    Sefflare['DFlux'].append(float(values[2]))

    # Convert lists inside the dicts to numpy arrays
    for key in keys:
        Sefflare[key] = np.array(Sefflare[key])

    ## Print the integtral flux spectrum as a table
    # print("Energy [MeV] | Differential Peak Flux [m^-2 sr-1 s^-1 ]")
    # for i in range(len(Sefflare['Energy'])):
    #     print(Sefflare['Energy'][i], Sefflare['IFlux'][i])

    # Convert the Fluxes from m-2 sr-1 s-1 to cm-2 s-1
    # m-2 to cm-2 --> 1e-4
    # sr-1 to 1 --> 4 pi
    Sefflare['IFlux'] = Sefflare['IFlux'] * 1e-4 * 4 * np.pi
    Sefflare['DFlux'] = Sefflare['DFlux'] * 1e-4 * 4 * np.pi

    return Sefflare


AtomicMass = [1.008, 4.003, 6.941, 9.012, 10.811, 12.011, 14.007, 15.999, 18.998, 20.18, 22.99, 24.305, 26.982,
              28.086, 30.974, 32.065, 35.453, 39.948, 39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938,
              55.845, 58.933, 58.693, 63.546, 65.39, 69.723, 72.64, 74.922, 78.96, 79.904, 83.8, 85.468, 87.62,
              88.906, 91.224, 92.906, 95.94, 98, 101.07, 102.906, 106.42, 107.868, 112.411, 114.818, 118.71, 121.76,
              127.6, 126.905, 131.293, 132.906, 137.327, 138.906, 140.116, 140.908, 144.24, 145, 150.36, 151.964,
              157.25, 158.925, 162.5, 164.93, 167.259, 168.934, 173.04, 174.967, 178.49, 180.948, 183.84, 186.207,
              190.23, 192.217, 195.078, 196.967, 200.59, 204.383, 207.2, 208.98, 209, 210, 222, 223, 226, 227,
              232.038, 231.036, 238.029]


def readSpenvis_sefflare_Ions(fileName, masses=AtomicMass):
    """
    Reads the ion peak flux table from Solar Flare spenvis_sefflare.txt files.

    File units:     Energy in MeV/nuc, Fluxes in m-2 sr-1 s-1 (MeV/nuc)-1
    Output units:   Energy in MeV, IFlux in cm-2 s-1, DFlux in cm-2 s-1 MeV-1

    Args:
        fileName (str): Path to the spenvis_sefflare.txt file.
        masses (list): Mass per species used for the MeV/nuc -> MeV conversion.
                       Defaults to natural-abundance atomic weights; pass integer
                       isotope masses (e.g. 56 for Fe-56) for Geant4 GPS consistency.

    Returns:
        Data[Z-1, DatapointNum, 0..2] = Energy, IFlux, DFlux
    """
    print("Reading in", fileName)

    Iontable = []
    readflag = 0
    with open(fileName, 'r') as f:
        for line in f:
            if readflag == 0:
                if "Differential Flux of ','SPECIES" in line:
                    readflag = 1
            elif readflag == 1:
                if "End of File" in line:
                    break
                else:
                    Iontable.append(np.fromstring(line, dtype=np.float64, sep=','))

    Iontable = np.array(Iontable)
    NumSpecies = 92

    EnergyPerNucleon = Iontable[:, 0]
    IonData = Iontable[:, 1:]

    # Convert the Fluxes from m-2 sr-1 s-1 to cm-2 s-1
    IonData = IonData * 1e-4 * 4 * np.pi

    Data = np.zeros((NumSpecies, len(Iontable), 3), dtype=float)
    for i in range(NumSpecies):
        Data[i, :, 0] = EnergyPerNucleon * masses[i]        # MeV/nuc -> MeV
        Data[i, :, 1] = IonData[:, i]                       # Integral Flux
        Data[i, :, 2] = IonData[:, NumSpecies + i] / masses[i]  # per MeV/nuc -> per MeV

    return Data


if __name__ == "__main__":
    File = "/home/anton/triton_work/Spectra/Carrington/GEO-Extreme/CREME96/spenvis_sefflare.txt"
    '''
    Results = readSpenvis_sefflare(File)

    # Plot the Integral peak Flux
    plt.figure(0)
    plt.plot(Results['Energy'], Results['IFlux'], '.-', label="Integral Peak Flux")
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel("Energy [MeV]")
    plt.ylabel("Integral Peak Flux [cm^-2 sr-1 s^-1]")
    plt.title("Integral Peak Flux")
    plt.legend()
    plt.grid()

    # Plot the Differential peak Flux
    plt.figure(1)
    plt.plot(Results['Energy'], Results['DFlux'], '.-', label="Differential Peak Flux")
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel("Energy [MeV]")
    plt.ylabel("Differential Peak Flux [cm^-2 s^-1 MeV^-1 sr^-1]")
    plt.title("Differential Peak Flux")
    plt.legend()
    plt.grid()

    ## Print the diffeential flux spectrum as a table
    #print("Energy [MeV] | Differential Peak Flux [cm^-2 s^-1 MeV^-1]")
    #for i in range(len(Results['Energy'])):
    #    print(Results['Energy'][i], Results['DFlux'][i])
    '''

    # Print differential ion peak flux tables for selected species,
    # using integer isotope masses (H-1, He-4, C-12, O-16, Fe-56) so the
    # MeV/nuc -> MeV conversion matches the isotope specified in the Geant4 GPS macro.
    SelectedSpecies = {'H-1': (0, 1), 'He-4': (1, 4), 'C-12': (5, 12), 'O-16': (7, 16), 'Fe-56': (25, 56)}

    IsotopeMasses = list(AtomicMass)
    for index, A in SelectedSpecies.values():
        IsotopeMasses[index] = A

    IonDataIso = readSpenvis_sefflare_Ions(File, masses=IsotopeMasses)

    for name, (Z, A) in SelectedSpecies.items():
        print(f"\n{name}: Energy [MeV] | Differential Peak Flux [cm^-2 s^-1 MeV^-1]")
        for E, D in zip(IonDataIso[Z, :, 0], IonDataIso[Z, :, 2]):
            print(f"{E:.17g} {D:.17g}")

    # Compare the ion species and plot the top candidates for integral flux around 100 MeV
    # (natural-abundance atomic weights are fine here, it is only a ranking)
    IonData = readSpenvis_sefflare_Ions(File)

    Species = ['H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl',
               'Ar', 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se',
               'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb',
               'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
               'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
               'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U ']

    EvalEnergy = 200  # MeV
    IntorDiff = 1  # 1 for Int 2 for Diff

    # Integral flux of each species interpolated at EvalEnergy
    FluxAtE = np.zeros(np.shape(IonData)[0])
    for Specie in range(np.shape(IonData)[0]):
        Energy = IonData[Specie, :, 0]
        if Energy[0] <= EvalEnergy <= Energy[-1]:
            FluxAtE[Specie] = np.interp(EvalEnergy, Energy, IonData[Specie, :, 1])

    TopSpecies = np.argsort(FluxAtE)[::-1][:5]

    print(f"\nTop 10 species by integral flux at {EvalEnergy} MeV:")
    for Specie in TopSpecies:
        print(f"Species: {Species[Specie]}, Integral Flux: {FluxAtE[Specie]:.3e} cm^-2 s^-1")
        plt.plot(IonData[Specie, :, 0], IonData[Specie, :, IntorDiff], label=Species[Specie])

    plt.axvline(EvalEnergy, color='grey', linestyle='--', linewidth=1)
    plt.xscale("log")
    plt.yscale("log")

    if IntorDiff == 1:
        plt.title("Solar Ion Integral Peak Flux")
        plt.ylabel("Integral Flux [cm-2 s-1]")
    elif IntorDiff == 2:
        plt.title("Solar Ion Differential Peak Flux")
        plt.ylabel("Differential Flux [cm-2 s-1 MeV-1]")

    plt.xlabel("Energy [MeV]")
    plt.legend()
    plt.grid(which='both')
    plt.show()
