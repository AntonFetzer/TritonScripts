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
    
    NewE = readGPSMacro("/l/triton_work/Spectra/FS2/GTO-11year-AE9.mac")
    OldE = readGPSMacro("/l/triton_work/Spectra/OLD/GTO/AE9/AE9Mission.mac")
    AE8 = readGPSMacro("/l/triton_work/Spectra/FS2/GTO-11year-AE8.mac")
    
    NewP = readGPSMacro("/l/triton_work/Spectra/FS2/GTO-11year-AP9.mac")    
    OldP = readGPSMacro("/l/triton_work/Spectra/OLD/GTO/AP9/AP9Mission.mac")
    AP8 = readGPSMacro("/l/triton_work/Spectra/FS2/GTO-11year-AP8.mac")

    # plt.plot(NewE['Energy'], NewE['Flux']/(11*12), label="FS2 AE9 11 year mean per month")
    plt.plot(OldE['Energy'], OldE['Flux'], label="97 percentile AE9 per month")

    # plt.plot(NewP['Energy'], NewP['Flux']/(11*12), label="FS2 AP9 11 year mean per month")
    plt.plot(OldP['Energy'], OldP['Flux'], label="97 percentile AP9 per month")

    plt.plot(AE8['Energy'], AE8['Flux']/(11*12), label="AE8 per month")
    plt.plot(AP8['Energy'], AP8['Flux']/(11*12), label="AP8 per month")

    plt.yscale("log")
    plt.xscale("log")
    plt.grid(which="both")
    plt.legend()
    plt.title("GPS")
    plt.xlabel("Kinetic energy [MeV]")
    plt.ylabel("Ingetgal or Differential Flux")


    # Plot ratio of old to new
    # plt.plot(OldE['Energy'], OldE['Flux']/(NewE['Flux']/(11*12)), label="AE9 Old / New")
    # plt.plot(OldP['Energy'], OldP['Flux']/(NewP['Flux']/(11*12)), label="AP9 Old / New")

    # plt.yscale("log")
    # plt.xscale("log")
    # plt.grid(which="both")
    # plt.legend()
    # plt.title("GPS Ratio Old to New")
    # plt.xlabel("Kinetic energy [MeV]")
    # plt.ylabel("Ratio")

    plt.show()
