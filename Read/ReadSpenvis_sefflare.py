import numpy as np
import matplotlib.pyplot as plt


def readSpenvis_sefflare(file):
    """
    Reads the proton peak flux from Solar Flare spenvis_sefflare.txt files.

    The spenvis_sefflare.txt file contains the proton peak flux data for a solar flare event.
    Energy in MeV
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


if __name__ == "__main__":
    File = "/u/02/fetzera1/unix/Desktop/triton_work/Spectra/Carrington/GEO-Extreme/SAPPHIRE100yearPeakFlux/spenvis_sefflare.txt"

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

    plt.show()
