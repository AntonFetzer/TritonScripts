import numpy as np
import matplotlib.pyplot as plt

def readGRASLETHistos(file):
    keys = ['lower', 'upper', 'mean', 'value', 'error', 'entries']

    # Initialize dictionaries for LET and Eff histograms
    LET = {key: [] for key in keys}
    Eff = {key: [] for key in keys}

    ReadFlag = 0

    # Open the file and read line by line
    with open(file, 'r') as f:
        for line in f:
            # Process the file based on the current state of ReadFlag
            if ReadFlag == 0:
                if "'Bin entries'" in line:
                    ReadFlag = 1
            elif ReadFlag == 1:
                if "'End of Block'" in line:
                    ReadFlag = 2
                else:
                    values = [float(x) for x in line.split(',')]
                    for i, key in enumerate(keys):
                        LET[key].append(values[i])
            elif ReadFlag == 2:
                if "'Bin entries'" in line:
                    ReadFlag = 3
            elif ReadFlag == 3:
                if "'End of Block'" in line:
                    ReadFlag = 4
                else:
                    values = [float(x) for x in line.split(',')]  # Adjust split as necessary
                    for i, key in enumerate(keys):
                        Eff[key].append(values[i])

    # Convert lists inside the dicts to numpy arrays
    for key in keys:
        LET[key] = np.array(LET[key])
        Eff[key] = np.array(Eff[key])

    ## Convert bin edges from [MeV/cm] to [MeV cm2 mg-1]
    C = 2330  # to convert from MeV/cm to Mev cm2 mg-1 in silicon with silicion density being 2.33 g/cm3

    for key in ['lower', 'upper', 'mean']:
        LET[key] = LET[key] / C
        Eff[key] = Eff[key] / C

    # The LET spectra are in 'counts' per MeV/cm bin ???
    # If the input spcectrum is in [cm-2 s-1] then the counts are in [s-1] per interface area between the shield and detector
    # This area is 1e6 cm2 in the 1Tile.gdml file written in December 2022
    # To get the values in [cm-2 s-1] they have to be divided by 1e6
    AreaNormFactor = 1e6

    for key in ['value', 'error']:
        LET[key] = LET[key] / AreaNormFactor
        Eff[key] = Eff[key] / AreaNormFactor

    # Adjust for mission duration by dividing by the mission duration in seconds
    # MissionDuration = 5 * 365.25 * 24 * 3600 # 5 years in seconds

    # for key in ['value', 'error']:
    #     LET[key] = LET[key] / MissionDuration
    #     Eff[key] = Eff[key] / MissionDuration

    return LET, Eff


if __name__ == "__main__":
    file = "/l/triton_work/LET/Foresail1-Hercules/FS1-SolarProtons/0mm/Res/LET_28556_4417.csv"

    LET, Eff = readGRASLETHistos(file)

    # For demonstration, replace np.shape(LETHist) with more appropriate operation
    print("LETHist Number of bins", len(LET['entries']))

    # Calculate sums and total LET by entries for LETHist
    NumberEntriesLETHist = np.sum(LET['entries'])
    TotalLETbyEntries = np.sum(LET['mean'] * LET['entries'])

    # Plotting LET Histogram by entries
    plt.figure(0)
    plt.bar(LET['lower'], LET['entries'], width=LET['upper'] - LET['lower'], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("LET Histogram " + f"{NumberEntriesLETHist:.2e}" + " entries\nTotal LET by Entries " + f"{TotalLETbyEntries:.2e}" + " [MeV cm2 mg-1]")
    plt.xlabel("LET [MeV cm2 mg-1]")
    plt.ylabel("Number of entries per LET bin")

    # Calculate total LET by values for LETHist
    TotalLETbyValues = np.sum(LET['mean'] * LET['value'])

    # Plotting LET Histogram by values and error bars
    plt.figure(1)
    plt.bar(LET['lower'], LET['value'], width=LET['upper'] - LET['lower'], align='edge', alpha=0.3)
    plt.errorbar(LET['mean'], LET['value'], LET['error'], fmt=' ', capsize=5, elinewidth=1, capthick=1, label="LET Histogram")
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("LET Histogram " + f"{NumberEntriesLETHist:.2e}" + " entries\nTotal LET by Values " + f"{TotalLETbyValues:.2e}" + " [MeV cm2 mg-1]")
    plt.xlabel("LET [MeV cm2 mg-1]")
    plt.ylabel("Rate per LET bin [cm-2 s-1]")

    # Calculate sums and total Eff by entries for EffHist
    NumberEntriesEffHist = np.sum(Eff['entries'])
    TotalEffbyEntries = np.sum(Eff['mean'] * Eff['entries'])

    # Plotting Eff Histogram by entries
    plt.figure(2)
    plt.bar(Eff['lower'], Eff['entries'], width=Eff['upper'] - Eff['lower'], align='edge')
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title(f"EffLET Histogram {NumberEntriesEffHist:.2e} entries\nTotal EffLET by Entries {TotalEffbyEntries:.2e} MeV cm2 mg-1")
    plt.xlabel("EffLET [MeV cm2 mg-1]")
    plt.ylabel("Number of entries per EffLET bin")

    # Calculate total EffLET by values for EffHist
    TotalEffLETbyValues = np.sum(Eff['mean'] * Eff['value'])

    # Plotting Eff Histogram by values
    plt.figure(3)
    plt.bar(Eff['lower'], Eff['value'], width=Eff['upper'] - Eff['lower'], align='edge', alpha=0.3)
    plt.errorbar(Eff['mean'], Eff['value'], Eff['error'], fmt=' ', capsize=5, elinewidth=1, capthick=1, label="Eff Histogram")
    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title(f"EffLET Histogram {NumberEntriesEffHist:.2e} entries\nTotal EffLET by Values {TotalEffLETbyValues:.2e} MeV cm2 mg-1")
    plt.xlabel("EffLET [MeV cm2 mg-1]")
    plt.ylabel("Rate per LET bin [cm-2 s-1]")

    plt.show()
