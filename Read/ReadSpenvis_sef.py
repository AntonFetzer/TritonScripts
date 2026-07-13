import numpy as np
import matplotlib.pyplot as plt


def integrated_spectral_energy(energy, differential_spectrum):
    """Return integral(E * dPhi/dE dE) over the available energy range."""
    valid = np.isfinite(energy) & np.isfinite(differential_spectrum)
    energy = energy[valid]
    differential_spectrum = differential_spectrum[valid]
    if energy.size < 2:
        return 0.0

    order = np.argsort(energy)
    energy = energy[order]
    differential_spectrum = differential_spectrum[order]
    return np.trapz(energy * differential_spectrum, energy)


def readSpenvis_sef_duration(fileName):
    """Return the total mission duration in seconds from the MIS_DUR header line.

    The tabulated SAPPHIRE fluences already include the per-energy attenuation
    and exposure-time effects, so dividing by the total mission duration gives
    the mission-averaged flux.
    """
    with open(fileName, 'r') as f:
        for line in f:
            if "'MIS_DUR'" in line:
                days = float(line.split(',')[2])
                return days * 86400.0
    raise ValueError(f"{fileName}: no 'MIS_DUR' header line found")


def readSpenvis_sef_protons(fileName):
    """
    Reads the SAPPHIRE solar proton fluence spectrum from spenvis_sef.txt files
    and converts it to mission-averaged flux.

    Args:
        file (str): Path to the spenvis_sef.txt file.

    Returns:
        dict: Dictionary with numpy arrays of flux spectrum:
            - 'Energy': Energy in MeV
            - 'IFlux': Integral Flux in cm-2 s-1
            - 'DFlux': Differential Flux in cm-2 s-1 MeV-1
            - 'Duration': Mission duration in s (multiply back for fluence)
    """
    keys = ['Energy', 'IFlux', 'DFlux']
    ProtonTable = {key: [] for key in keys}

    Duration = readSpenvis_sef_duration(fileName)

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
                    ProtonTable['IFlux'].append(float(values[1]) / Duration)
                    ProtonTable['DFlux'].append(float(values[2]) / Duration)

    # Convert lists inside the dicts to numpy arrays
    for key in keys:
        ProtonTable[key] = np.array(ProtonTable[key])

    ProtonTable['Duration'] = Duration

    return ProtonTable




AtomicMass = [1.008, 4.003, 6.941, 9.012, 10.811, 12.011, 14.007, 15.999, 18.998, 20.18, 22.99, 24.305, 26.982,
              28.086, 30.974, 32.065, 35.453, 39.948, 39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938,
              55.845, 58.933, 58.693, 63.546, 65.39, 69.723, 72.64, 74.922, 78.96, 79.904, 83.8, 85.468, 87.62,
              88.906, 91.224, 92.906, 95.94, 98, 101.07, 102.906, 106.42, 107.868, 112.411, 114.818, 118.71, 121.76,
              127.6, 126.905, 131.293, 132.906, 137.327, 138.906, 140.116, 140.908, 144.24, 145, 150.36, 151.964,
              157.25, 158.925, 162.5, 164.93, 167.259, 168.934, 173.04, 174.967, 178.49, 180.948, 183.84, 186.207,
              190.23, 192.217, 195.078, 196.967, 200.59, 204.383, 207.2, 208.98, 209, 210, 222, 223, 226, 227,
              232.038, 231.036, 238.029]

def readSpenvis_sef_Ions(fileName, masses=AtomicMass):
    """Read the 92-species SAPPHIRE ion fluence table and return mission-averaged
    flux as Data[Z-1, point, (Energy MeV, IFlux cm-2 s-1, DFlux cm-2 s-1 MeV-1)]."""
    print("Reading in", fileName)

    Duration = readSpenvis_sef_duration(fileName)

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
        # The sef.txt tabulates the mission Fluence (attenuation and exposure
        # already included); divide by the mission duration to get the flux.
        IonData[i] = line[1:cols] / Duration

    # print(Energy)

    Data = np.zeros((NumSpecies, len(Iontable), 3), dtype=float)

    for i in range(NumSpecies):
        Energy = EnergyPerNucleon*masses[i]
        #print("Species:", i)
        #print("AtomicMass:", AtomicMass[i])
        #print("Energy:", Energy)
        #print("IonData:", IonData[:, i+1])
        Data[i, :, 0] = Energy
        Data[i, :, 1] = IonData[:, i]  # Integral Flux
        Data[i, :, 2] = IonData[:, NumSpecies+i] / masses[i]  # per MeV/nuc -> per MeV

    return Data
# Data[ Z-NUmber , DatapointNum, Energy or Integral Flux or Differential Flux ]


if __name__ == "__main__":

    File = "/home/anton/triton_work/Spectra/Carrington/GEO/spenvis_sef.txt"

    # Use integer isotope masses so the MeV/nuc -> MeV conversion matches
    # the isotope specified in the Geant4 GPS macro.
    SelectedSpecies = {'H-1': (0, 1), 'He-4': (1, 4), 'C-12': (5, 12), 'O-16': (7, 16),
                       'Fe-56': (25, 56), 'Si-28': (13, 28), 'Ca-40': (19, 40), 'Ni-58': (27, 58)}

    IsotopeMasses = list(AtomicMass)
    for index, A in SelectedSpecies.values():
        IsotopeMasses[index] = A

    IonDataIso = readSpenvis_sef_Ions(File, masses=IsotopeMasses)

    for name, (Z, A) in SelectedSpecies.items():
        print(f"\n{name}: Energy [MeV] | Differential Flux [cm^-2 s^-1 MeV^-1]")
        for E, D in zip(IonDataIso[Z, :, 0], IonDataIso[Z, :, 2]):
            print(f"{E:.17g} {D:.17g}")

    # Compare species using natural-abundance atomic weights for ranking.
    IonData = readSpenvis_sef_Ions(File)

    Species = ['H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl',
               'Ar', 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se',
               'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb',
               'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
               'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
               'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U ']

    IntorDiff = 1  # 1 for Int 2 for Diff

    # Rank spectra by total particle energy flux, not by their value at one
    # arbitrarily selected energy: integral(E * differential flux dE).
    SpectralEnergy = np.zeros(np.shape(IonData)[0])
    for Specie in range(np.shape(IonData)[0]):
        SpectralEnergy[Specie] = integrated_spectral_energy(
            IonData[Specie, :, 0], IonData[Specie, :, 2]
        )

    TopSpecies = np.argsort(SpectralEnergy)[::-1][:20]

    print(f"\nTop {len(TopSpecies)} species by total energy flux:")
    for Specie in TopSpecies:
        print(f"Species: {Species[Specie]}, Energy Flux: {SpectralEnergy[Specie]:.3e} MeV cm^-2 s^-1")
        plt.plot(IonData[Specie, :, 0], IonData[Specie, :, IntorDiff], label=Species[Specie])

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
