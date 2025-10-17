import numpy as np
import matplotlib.pyplot as plt
import csv


def readGPSMacro(file):

    keys = ['Energy', 'Flux']
    # Generating the dictionary from the list of keys
    GPS = {key: [] for key in keys}

    ReadFlag = 0

    print("Reading in File: " + file)
    with open(file, 'r') as f:
        for line in f:
            if ReadFlag == 0:
                if '/gps/hist/type arb' in line:
                    ReadFlag = 1
            elif ReadFlag == 1:
                if "/gps/hist/inter" in line:
                    ReadFlag = 2
                else:
                    # Split the line into parts based on the comma
                    values = line.split('  ')
                    GPS['Energy'].append(float(values[1]))
                    GPS['Flux'].append(float(values[2]))

    # Convert lists inside the dicts to numpy arrays
    for key in keys:
        GPS[key] = np.array(GPS[key])

    return GPS


if __name__ == "__main__":
    
    Data = readGPSMacro("/l/triton_work/Spectra/Carrington/Electron/CarringtonElectronINTEGRALPowTabelated.mac")
    
    plt.plot(Data['Energy'], Data['Flux'])

    plt.yscale("log")
    plt.xscale("log")
    plt.grid(which="both")
    plt.legend()
    plt.title("GPS")
    plt.xlabel("Kinetic energy [MeV]")
    plt.ylabel("Ingetgal or Differential Flux")

    plt.show()
