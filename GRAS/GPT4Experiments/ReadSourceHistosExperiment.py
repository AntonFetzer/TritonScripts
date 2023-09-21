import os
import pandas as pd
import matplotlib.pyplot as plt


def readSourceHistos(file):
    # Find start and end of the block with histogram data
    start_line = None
    end_line = None
    with open(file, 'r') as f:
        for i, line in enumerate(f):
            if "'Bin entries'" in line:
                start_line = i + 1  # The actual data starts from the next line
            if "'End of Block'" in line and start_line is not None:
                end_line = i
                break

    if start_line is None or end_line is None:
        raise ValueError("Could not locate the data block in the file.")

    # Read only the relevant block of data
    hist_data = pd.read_csv(file, skiprows=start_line, nrows=end_line - start_line, header=None, engine='python')
    hist_data.columns = ['lower', 'upper', 'mean', 'value', 'error', 'entries']

    return hist_data


if __name__ == "__main__":
    Path = "/l/triton_work/SourceHistograms/Test/ISO-GTO-Fe/Res/"
    File = "Source_OnlyEnergy.csv"
    FullPath = os.path.join(Path, File)

    data = readSourceHistos(FullPath)

    # Plot the histogram
    plt.bar(data['lower'], data['value'], width=data['upper'] - data['lower'], yerr=data['error'], align='edge')

    plt.xlabel('Energy [MeV]')
    plt.ylabel('Counts')
    plt.title('Kinetic Energy Spectrum')
    plt.xscale('log')  # because the 'X_AXIS_SCALE' is 'log'
    plt.yscale('log')
    plt.show()
