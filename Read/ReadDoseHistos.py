import numpy as np
import matplotlib.pyplot as plt
import csv


def readDoseHistos(file):
    """
    Read the dose histograms from a .csv file and store them in a dictionary
    Args:
        file: Path to the .csv file containing the dose histograms
    Returns:
        Tuple of two dictionaries:
            - DoseHistDict: Dictionary containing the dose histogram
            - PrimaryHistDict: Dictionary containing the primary energy histogram
    """
    
    DoseHistDict = {
        'lower': [],    #   rad/s
        'upper': [],    #   rad/s
        'mean': [],     #   rad/s
        'value': [],    #   counts/cm2
        'error': [],    #   counts/cm2
        'entries': []   #   Num
    }

    PrimaryHistDict = {
        'lower': [],    #   MeV
        'upper': [],    #   MeV
        'mean': [],     #   MeV
        'value': [],    #   rad/s
        'error': [],    #   rad/s
        'entries': []   #   Num
    }

    ReadFlag = 0
    with open(file, 'r') as f:
        reader = csv.reader(f)
        for line in reader:
            if ReadFlag == 0:
                if "'Bin entries'" in line:
                    ReadFlag = 1
            elif ReadFlag == 1:
                if "'End of Block'" in line:
                    ReadFlag = 2
                else:
                    for i, key in enumerate(DoseHistDict):
                        DoseHistDict[key].append(float(line[i]))
            elif ReadFlag == 2:
                if "'Bin entries'" in line:
                    ReadFlag = 3
            elif ReadFlag == 3:
                if "'End of Block'" in line:
                    ReadFlag = 4
                else:
                    for i, key in enumerate(PrimaryHistDict):
                        PrimaryHistDict[key].append(float(line[i]))

    # Convert lists to numpy arrays
    for key in DoseHistDict:
        DoseHistDict[key] = np.array(DoseHistDict[key])
        PrimaryHistDict[key] = np.array(PrimaryHistDict[key])

    # Dose data is in rad/s --> multiply with number of seconds in a month to get to dose per months.
    # Dose is given per generated particle --> need to divide by the number of files
    # Dose is given in rad --> divide by 1000 to get to krad.
    ScaleFactor = 30 * 24 * 60 * 60 / 1000
    DoseHistDict['lower'] *= ScaleFactor
    DoseHistDict['upper'] *= ScaleFactor
    DoseHistDict['mean'] *= ScaleFactor

    PrimaryHistDict['value'] *= ScaleFactor
    PrimaryHistDict['error'] *= ScaleFactor

    return DoseHistDict, PrimaryHistDict



if __name__ == "__main__":
    path = "/l/triton_work/Histograms/Cobalt-60/Res/"
    file = "Electrons_627965_97150.csv"

    DoseHist, PrimaryHist = readDoseHistos(path + file)

    #print("DoseHist Shape", np.shape(PrimaryHist))
    #PrimaryHist = merge_bins(PrimaryHist, 10)
    #print("DoseHist Shape", np.shape(PrimaryHist))

    NumberEntries = sum(DoseHist['entries'])
    TotalDose = sum(DoseHist['mean'] * DoseHist['entries'])

    plt.figure(1)
    plt.bar(DoseHist['lower'], DoseHist['entries'], width=DoseHist['upper'] - DoseHist['lower'], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("Dose Depositions Histogram\n" + f"{NumberEntries:.2}" + " entries " + f"{TotalDose:.2}" + " krad total dose ?!?")
    plt.xlabel("Dose [krad per Month]")
    plt.ylabel("Number of entries per dose bin")

    NumberEntries = sum(PrimaryHist['entries'])
    TotalDose = sum(PrimaryHist['value'])

    plt.figure(2)
    plt.bar(PrimaryHist['lower'], PrimaryHist['value'], width=PrimaryHist['upper'] - PrimaryHist['lower'], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("Dose deposited VS primary kinetic energy\n" + f"{NumberEntries:.2}" + " enties " + f"{TotalDose:.2}" + " krad total dose")
    plt.xlabel("Kinetic energy [MeV]")
    plt.ylabel("Dose [krad per Month]")

    NumberEntries = sum(PrimaryHist['entries'])

    plt.figure(3)
    plt.bar(PrimaryHist['lower'], PrimaryHist['entries'], width=PrimaryHist['upper'] - PrimaryHist['lower'],
            align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("Particle count VS primary kinetic energy\n" + f"{NumberEntries:.2}" + " enties ")
    plt.xlabel("Kinetic energy [MeV]")
    plt.ylabel("Number of entries")

    #plt.savefig(path + "../Hist.pdf", format='pdf', bbox_inches="tight")
    plt.show()

