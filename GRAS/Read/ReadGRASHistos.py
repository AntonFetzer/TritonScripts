import numpy as np
import matplotlib.pyplot as plt
import csv


def readGRASHistos(file):

    DoseHist = []
    PrimaryHist = []

    ReadFlag = 0
    #print("Reading in File: " + file)
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
                    DoseHist.append([float(x) for x in line])
            elif ReadFlag == 2:
                if "'Bin entries'" in line:
                    ReadFlag = 3
            elif ReadFlag == 3:
                if "'End of Block'" in line:
                    ReadFlag = 4
                else:
                    PrimaryHist.append([float(x) for x in line])

    DoseHist = np.asarray(DoseHist)
    PrimaryHist = np.asarray(PrimaryHist)

    # Dose data is in rad/s --> multiply with number of seconds in a month to get to dose per months.
    # Dose is given per generated particle --> need to divide by the number of files
    # Dose is given in rad --> divide by 1000 to get to krad.
    ScaleFactor = 30 * 24 * 60 * 60 / 1000

    DoseHist[:, 0] = DoseHist[:, 0] * ScaleFactor
    DoseHist[:, 1] = DoseHist[:, 1] * ScaleFactor
    DoseHist[:, 2] = DoseHist[:, 2] * ScaleFactor

    PrimaryHist[:, 3] = PrimaryHist[:, 3] * ScaleFactor
    PrimaryHist[:, 4] = PrimaryHist[:, 4] * ScaleFactor

    return DoseHist, PrimaryHist



if __name__ == "__main__":
    file = "/home/anton/Desktop/triton_work/CARRINGTON/HistogramAE9Test/Res0mm/Electrons1Tile0mm_179_27.csv"

    DoseHist, PrimaryHist = readGRASHistos(file)

    print("DoseHist Shape", np.shape(DoseHist))

                    #   Dose            Primary
    lowerID = 0     #   rad/s           MeV
    upperID = 1     #   rad/s           MeV
    meanID = 2      #   rad/s           MeV
    valueID = 3     #   counts/cm2      rad/s
    errorID = 4     #   counts/cm2      rad/s
    entriesID = 5   #   Num             Num

    # plt.pyplot.bar(x, height, width=0.8, bottom=None, *, align='center', DoseHist=None, **kwargs)

    NumberEntries = sum(DoseHist[:, entriesID])
    TotalDose = sum(DoseHist[:, meanID] * DoseHist[:, entriesID])
    '''
    plt.figure(1)
    plt.bar(DoseHist[:, lowerID], DoseHist[:, entriesID], width=DoseHist[:, upperID] - DoseHist[:, lowerID], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("Dose Depositions Histogram\n" + f"{NumberEntries:.2}" + " entries " + f"{TotalDose:.2}" + " krad total dose ?!?")
    plt.xlabel("Dose [krad per Month]")
    plt.ylabel("Number of entries per dose bin")
    '''

    NumberEntries = sum(PrimaryHist[:, entriesID])
    TotalDose = sum(PrimaryHist[:, valueID])

    plt.figure(2)
    plt.bar(PrimaryHist[:, lowerID], PrimaryHist[:, valueID], width=PrimaryHist[:, upperID] - PrimaryHist[:, lowerID], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("Dose deposited VS primary kinetic energy\n" + f"{NumberEntries:.2}" + " enties " + f"{TotalDose:.2}" + " krad total dose")
    plt.xlabel("Kinetic energy [MeV]")
    plt.ylabel("Dose [krad per Month]")

    NumberEntries = sum(PrimaryHist[:, entriesID])

    plt.figure(3)
    plt.bar(PrimaryHist[:, lowerID], PrimaryHist[:, entriesID], width=PrimaryHist[:, upperID] - PrimaryHist[:, lowerID], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("Particle count VS primary kinetic energy\n" + f"{NumberEntries:.2}" + " enties ")
    plt.xlabel("Kinetic energy [MeV]")
    plt.ylabel("Number of entries")

    plt.show()
