import os
import pandas as pd
import matplotlib.pyplot as plt


def readSourceHistos(file):
    # List of known histograms in the expected order
    histo_names = ['Energy', 'Momentum', 'Phi', 'Theta', 'Weights']

    # Initialize a dictionary to store the histograms
    histograms = {}

    # Initialize tracking variables
    start_line = None
    end_line = None
    current_histo_idx = 0  # Index to track the current histogram name in histo_names list

    with open(file, 'r') as f:
        for i, line in enumerate(f):
            if "'Bin entries'" in line:
                # Start of a histogram data block
                start_line = i + 1
            elif "'End of Block'" in line and start_line is not None:
                # End of a histogram data block
                end_line = i
                hist_data = pd.read_csv(file, skiprows=start_line, nrows=end_line - start_line, header=None, engine='python')
                hist_data.columns = ['lower', 'upper', 'mean', 'value', 'error', 'entries']

                # Convert columns to numpy arrays and store in dictionary
                hist_dict = {col: hist_data[col].to_numpy() for col in hist_data.columns}

                # Check if we have reached beyond the known histogram names
                if current_histo_idx >= len(histo_names):
                    raise ValueError("More histogram blocks in file than expected.")

                # Assign the current histogram name from the list of known names
                current_histo_name = histo_names[current_histo_idx]

                # Store the data in the dictionary
                histograms[current_histo_name] = hist_dict

                # Reset tracking variables for next histogram block
                start_line = None
                end_line = None

                # Move to the next histogram name in the list
                current_histo_idx += 1

    if len(histograms) == 0:
        raise ValueError("Could not locate any data blocks in the file.")

    # Check if we found all expected histograms
    if current_histo_idx < len(histo_names):
        raise ValueError("Fewer histogram blocks in file than expected.")

    return histograms


if __name__ == "__main__":
    Path = "/l/triton_work/SourceHistograms/Carrington/CarringtonElectronINTEGRALPowTabelated/Res/"
    # Find the first .csv file in the specified directory
    FullPath = next((os.path.join(Path, f) for f in os.listdir(Path) if f.endswith('.csv')), None)
    # Print the full path of the file
    print("Reading in file: ", FullPath)
    # Construct the path for the 'Plot' directory
    PlotDir = os.path.join(os.path.dirname(Path), '../Plot')

    histograms = readSourceHistos(FullPath)

    # Set the grid to be behind the plot
    plt.rc('axes', axisbelow=True)

    # Plot the 'Energy' histogram
    data = histograms['Energy']

    plt.figure(0)
    plt.bar(data['lower'], data['value'], width=data['upper'] - data['lower'], yerr=data['error'], align='edge')
    plt.xlabel('Energy [MeV]')
    plt.ylabel('Counts per bin')
    plt.title('Kinetic Energy Spectrum')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(which='both')

    print("Saving Energy Histogram: ", os.path.join(PlotDir, 'Energy Source Spectrum.pdf'))
    plt.savefig(os.path.join(PlotDir, 'Energy Source Spectrum.pdf'), format='pdf', bbox_inches='tight')

    # Plot the 'Momentum' histogram
    data = histograms['Momentum']
    plt.figure(1)
    plt.bar(data['lower'], data['value'], width=data['upper'] - data['lower'], yerr=data['error'], align='edge')
    plt.xlabel('Momentum [MeV]')
    plt.ylabel('Counts per bin')
    plt.title('Momentum Spectrum')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(which='both')

    print("Saving Momentum Histogram: ", os.path.join(PlotDir, 'Momentum Source Spectrum.pdf'))
    plt.savefig(os.path.join(PlotDir, 'Momentum Source Spectrum.pdf'), format='pdf', bbox_inches='tight')

    # Plot the 'Phi' histogram
    data = histograms['Phi']
    plt.figure(2)
    plt.bar(data['lower'], data['value'], width=data['upper'] - data['lower'], yerr=data['error'], align='edge')
    plt.xlabel('Phi [degrees]')
    plt.ylabel('Counts per bin')
    plt.title('Phi Distribution')
    plt.grid(which='both')

    print("Saving Phi Histogram: ", os.path.join(PlotDir, 'Phi Source Distribution.pdf'))
    plt.savefig(os.path.join(PlotDir, 'Phi Source Distribution.pdf'), format='pdf', bbox_inches='tight')

    # Plot the 'Theta' histogram
    data = histograms['Theta']
    plt.figure(3)
    plt.bar(data['lower'], data['value'], width=data['upper'] - data['lower'], yerr=data['error'], align='edge')
    plt.xlabel('Theta [degrees]')
    plt.ylabel('Counts per bin')
    plt.title('Theta Distribution')
    plt.grid(which='both')

    print("Saving Theta Histogram: ", os.path.join(PlotDir, 'Theta Source Distribution.pdf'))
    plt.savefig(os.path.join(PlotDir, 'Theta Source Distribution.pdf'), format='pdf', bbox_inches='tight')

    # Plot the 'Weights' histogram
    data = histograms['Weights']
    plt.figure(4)
    plt.bar(data['lower'], data['value'], width=data['upper'] - data['lower'], yerr=data['error'], align='edge')
    plt.xlabel('Weights')
    plt.ylabel('Counts')
    plt.title('Weights Distribution')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(which='both')

    print("Saving Weights Histogram: ", os.path.join(PlotDir, 'Weights Source Distribution.pdf'))
    plt.savefig(os.path.join(PlotDir, 'Weights Source Distribution.pdf'), format='pdf', bbox_inches='tight')

    #plt.show()

