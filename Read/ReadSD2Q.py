import numpy as np
import matplotlib.pyplot as plt

def readSDQ2(fileName):
    """
    Read the TID vs shielding thickness curves from SHIELDOSE-2Q and store them in a dictionary of numpy arrays

    Args:
        file: Path to the .csv file containing the shielding curves

    Returns:
        Dictionary containing the TID vs shielding thickness curves as 1D numpy arrays
    """
    ListofCollumns = ['Thickness', 'Total','Electrons','Bremsstrahlung','Trapped Protons','Solar Protons']

    # Generating the dictionary from the list of Collumns
    SDData = {key: [] for key in ListofCollumns}

    print("Reading in", fileName)
    f = open(fileName, "r")

    ReadFlag = 0

    for line in f:
        if ReadFlag == 0:
            if "Dose in Si" in line:
                ReadFlag = 1
        elif ReadFlag == 1:
            if "End of File" in line:
                ReadFlag = 2
            else:
                values = line.split(',')
                for i in range(len(values)):
                    SDData[ListofCollumns[i]].append(float(values[i]))

    # Convert lists inside the dicts to numpy arrays
    for key in ListofCollumns:
        SDData[key] = np.array(SDData[key])

    # The default dose unit in SHIELDOSE-2Q is rad.
    # --> divide by 1000 to get to krad.
    for key in ['Total','Electrons','Bremsstrahlung','Trapped Protons','Solar Protons']:
        SDData[key] = SDData[key] / 1000

    return SDData


if __name__ == "__main__":
    Data = readSDQ2("/l/triton_work/Spectra/OLD/GTO/Shieldose/spenvis_sqo.txt")

    print("Aluminium Thicknesses", Data['Thickness'])
    print("Total Dose", Data['Total'])
    print("Electrons", Data['Electrons'])
    print("Bremsstrahlung", Data['Bremsstrahlung'])
    print("Trapped Protons", Data['Trapped Protons'])
    print("Solar Protons", Data['Solar Protons'])

    plt.plot(Data['Thickness'], Data['Electrons'], '.-', label="Electrons")
    plt.plot(Data['Thickness'], Data['Bremsstrahlung'], '.-', label="Bremsstrahlung")
    plt.plot(Data['Thickness'], Data['Trapped Protons'], '.-', label="Trapped Protons")
    if len(Data['Solar Protons']) > 0:
        plt.plot(Data['Thickness'], Data['Solar Protons'], '.-', label="Solar Protons")
    plt.plot(Data['Thickness'], Data['Total'], '.-', label="Total Dose")

    plt.yscale("log")
    # plt.xscale("log")
    plt.grid(which='major')
    plt.title("Shielddose")
    plt.xlabel("Aluminium Shield Thickness [mm]")
    plt.ylabel("Total ionizing dose in silicon [krad]")
    plt.legend()

    plt.show()