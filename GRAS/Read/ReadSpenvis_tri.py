import numpy as np
import matplotlib.pyplot as plt
import csv


def readSpenvis_tri(file):

    keys = ["Energy", "Integral", "Differential"]
    Protons = {key: [] for key in keys}
    Electrons = {key: [] for key in keys}

    ReadFlag = 0

    # Open the file and read line by line
    print("Reading in File: " + file)
    with open(file, 'r') as f:
        for line in f:
            if ReadFlag == 0 and "'Differential Flux'" in line:
                ReadFlag = 1
            elif ReadFlag == 1:
                if "'End of Block'" in line:
                    ReadFlag = 2
                else:
                    # Split the line into parts based on the comma
                    values = line.split(',')
                    Protons["Energy"].append(float(values[0]))
                    Protons["Integral"].append(float(values[1]))
                    Protons["Differential"].append(float(values[2]))
            elif ReadFlag == 2 and "'Differential Flux'" in line:
                ReadFlag = 3
            elif ReadFlag == 3:
                if "'End of File'" in line:
                    ReadFlag = 4
                else:
                    # Split the line into parts based on the comma
                    values = line.split(',')
                    Electrons["Energy"].append(float(values[0]))
                    Electrons["Integral"].append(float(values[1]))
                    Electrons["Differential"].append(float(values[2]))

    # Remove trailing zeros
    # Check if in the last row collumn 1 or 2 contain zeros
    # If yes, remove the last row
    if Protons["Energy"][-1] == 0 or Protons["Integral"][-1] == 0:
        for key in keys:
            Protons[key].pop(-1)
    
    if Electrons["Energy"][-1] == 0 or Electrons["Integral"][-1] == 0:
        for key in keys:
            Electrons[key].pop(-1)

    # Convert lists inside the dicts to numpy arrays
    for key in keys:
        Protons[key] = np.array(Protons[key])
        Electrons[key] = np.array(Electrons[key])

    # Print out the proton data
    for i in range(len(Protons["Energy"])):
        print(f"{Protons['Energy'][i]:.2f} {Protons['Integral'][i]:.2e} {Protons['Differential'][i]:.2e}")

    # Print out the electron data
    for i in range(len(Electrons["Energy"])):
        print(f"{Electrons['Energy'][i]:.2f} {Electrons['Integral'][i]:.2e} {Electrons['Differential'][i]:.2e}")

    return Protons, Electrons



if __name__ == "__main__":

    file = "/l/triton_work/Spectra/ISS/spenvis_tri.txt"

    Protons, Electrons = readSpenvis_tri(file)

    plt.figure(1)
    plt.plot(Protons['Energy'], Protons['Integral'], label="Integral Protons")
    plt.plot(Electrons['Energy'], Electrons['Integral'], label="Integral Electrons")
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Energy [MeV]")
    plt.ylabel("Integral Flux [cm-2 s-1]")
    plt.title("Integral Flux of Protons and Electrons")
    plt.grid(which='both')
    plt.legend()

    plt.figure(2)
    plt.plot(Protons['Energy'], Protons['Differential'], label="Differential Protons")
    plt.plot(Electrons['Energy'], Electrons['Differential'], label="Differential Electrons")
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Energy [MeV]")
    plt.ylabel("Differential Flux [cm-2 s-1 MeV-1]")
    plt.title("Differential Flux of Protons and Electrons")
    plt.grid(which='both')
    plt.legend()

    plt.show()
