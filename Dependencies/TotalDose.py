import os
import numpy as np
from Read.ReadDose import readDose
import matplotlib.pyplot as plt
import sys
from uncertainties import ufloat, ufloat_fromstr


def totalDose(path, filename_contains=None, replica_error=False, min_replicas=5):
    """
    Reads the TID values of all sensitive volumes of all GRAS CSV output files in a folder and calculates the total dose for each of the tiles.

    Args:
        path (str): The path to the folder containing the CSV files.
        filename_contains (str or None): Optional substring used to select CSV
            filenames in legacy result folders that contain several particle
            datasets. Prefer passing a directory containing only one dataset.
        replica_error (bool): If False (default), 'error' is pooled from the
            errors GRAS quotes in each file. If True, 'error' is the error of the
            entry-weighted mean estimated from the between-file scatter of the
            replica files instead. Use this for importance-weighted sources
            (GPS '/gps/hist/type biase'): GRAS 5.02.01 GRASStatDouble::fill
            treats the event weight as a frequency weight, so its quoted errors
            are wrong whenever the weights are not 1. The files must be
            replicas of one configuration that differ only by their seed.
        min_replicas (int): Minimum number of files required when
            replica_error is True; fewer raises ValueError. The scatter-based
            error is itself uncertain by about 1/sqrt(2*(files-1)).

    Returns:
        dict: A dictionary containing the following keys:
            - 'dose': A numpy array containing the total dose values in kRad for each tile.
            - 'error': A numpy array containing the absolute total dose error in kRad for each tile.
            - 'entries': A numpy array containing the number of entries for each tile.
            - 'non-zeros': A numpy array containing the number of non-zero entries for each tile.
    """

    print("\nReading in all csv files in folder:", path)

    # Get list of all csv files in Path
    Files = [
        f for f in os.listdir(path)
        if f.endswith('.csv')
        and (filename_contains is None or filename_contains in f)
    ]

    if not Files:
        # Raise an error if no files are found
        selector = (
            f" containing {filename_contains!r}" if filename_contains is not None else ""
        )
        sys.exit("ERROR !!! No csv files" + selector + " found in folder: " + path)

    if replica_error and len(Files) < max(min_replicas, 2):
        raise ValueError(
            f"replica_error needs at least {max(min_replicas, 2)} replica files, "
            f"found {len(Files)} in {path}")

    keys = ['dose', 'error', 'entries', 'non-zeros']
    TID = {key: [] for key in keys}
    # Ionising dose in kRad
    # Number of entries is the same for each tile, because each detected particle causes an entry at all tiles. Most of these entries are therefore empty.
    # Number of non-zero entries is different for each tile.

    # RawData is a list of dictionaries.
    # Each dictionary contains the dose, error, entries and non-zero entries for each file.
    RawData = [readDose(os.path.join(path, file)) for file in Files]

    # Initialize numpy arrays for dose, error, entries, and non-zeros using the first file to get the correct shape.
    for key in keys:
        TID[key] = np.zeros_like(RawData[0][key])

    for raw in RawData:
        #print(raw['dose'] * raw['entries'])
        TID['dose'] += raw['dose'] * raw['entries']  # Weighted average of all the dose results. Multiply dose by number of entries, then later divide by the total number of entries.
        TID['entries'] += raw['entries']  # Summed up number of entries from all files
        TID['non-zeros'] += raw['non-zeros']  # Summed up number of Non Zero Entries from all files

    # Calculate the weighted average of the dose.
    TID['dose'] = TID['dose'] / TID['entries']  # Divide by summed up number of entries to get the weighted average of the dose results

    # Pooled statistical error: err = sqrt(sum(N_i^2 * sigma_i^2)) / sum(N_i).
    # The between-file scatter term of the full ANOVA identity would enter with
    # weight N_i (never N_i^2, which double-counts the statistical error and
    # inflates the result by up to sqrt(2)); at that weight it is a relative
    # O(1/N_i) correction (~1e-10 for typical runs), so it is dropped here and
    # the scatter is used as a consistency diagnostic below instead.
    InternalErrSq = np.zeros_like(TID['error'])
    ExternalErrSq = np.zeros_like(TID['error'])
    for raw in RawData:
        InternalErrSq += raw['entries'] ** 2 * raw['error'] ** 2
        ExternalErrSq += raw['entries'] ** 2 * (raw['dose'] - TID['dose']) ** 2

    TID['error'] = np.sqrt(InternalErrSq) / TID['entries']

    # Consistency diagnostic (Birge ratio): for files that differ only by their
    # random seed, the observed between-file scatter is fully explained by the
    # quoted statistical errors, so external/internal ~ 1.  A ratio well above 1
    # means the files genuinely differ (configuration, build, or a corrupt csv)
    # -- a bug signal, not extra statistics, so halt and let the user look.
    if len(RawData) > 1:
        Birge = np.sqrt(np.divide(ExternalErrSq, InternalErrSq,
                                  out=np.zeros_like(InternalErrSq), where=InternalErrSq > 0))
        MaxBirge = max(Birge)
        print("Maximum Birge ratio (scatter / quoted errors):", f"{MaxBirge:.2f}", "at Tile number", np.argmax(Birge))
        if MaxBirge > 2:
            print("ERROR !!! Inconsistent files: the between-file scatter exceeds the quoted statistical errors")
            print("  The files differ by more than their random seed, or a csv file is corrupt.")
            # Prompt only on a terminal: an open pipe without input (ssh, a pipeline)
            # never delivers EOF and would block input() indefinitely.
            if sys.stdin.isatty():
                input("  Press Enter to continue anyway...")
            else:
                print("  (non-interactive run, continuing)")

    # Error of the entry-weighted mean from the replica scatter. File i is the
    # mean of N_i independent events, Var(y_i) = s^2 / N_i, and
    # E[sum(N_i * (y_i - mean)^2)] = (k - 1) * s^2 for k files, so
    #   Var(mean) = s^2 / N = sum(N_i * (y_i - mean)^2) / ((k - 1) * N).
    # This needs no per-file error, so it holds for weighted events whose
    # quoted GRAS error does not. Weight N_i, not N_i^2 as in the Birge sums.
    if len(RawData) > 1:
        ScatterSq = np.zeros_like(TID['error'])
        for raw in RawData:
            ScatterSq += raw['entries'] * (raw['dose'] - TID['dose']) ** 2
        ReplicaError = np.sqrt(ScatterSq / ((len(RawData) - 1) * TID['entries']))

        # Quoted errors far above the scatter mean the per-file errors are
        # overestimated, as GRAS does for importance weights below 1. For
        # consistent files each tile's squared ratio is ~chi2(k-1)/(k-1), so the
        # mean over T tiles is ~chi2(T(k-1))/(T(k-1)). Warn when the RMS ratio is
        # below 0.5 and below the 0.1 % quantile of that distribution
        # (Wilson-Hilferty), which keeps chance low scatter of few files quiet.
        Quoted = TID['error'] > 0
        if np.any(Quoted):
            RMSRatio = np.sqrt(np.mean((ReplicaError[Quoted] / TID['error'][Quoted]) ** 2))
            Dof = np.count_nonzero(Quoted) * (len(RawData) - 1)
            Spread = np.sqrt(2 / (9 * Dof))
            LowQuantile = max(1 - Spread ** 2 - 3.09 * Spread, 0) ** 3
            LowRatio = RMSRatio < min(0.5, np.sqrt(LowQuantile))
        else:
            RMSRatio, LowRatio = np.nan, False
        if replica_error or LowRatio:
            print("RMS ratio of replica-scatter to quoted error:", f"{RMSRatio:.2f}",
                  f"({len(RawData)} files)")
        if LowRatio:
            print("Warning !!! The quoted statistical errors are overestimated: the between-file scatter is far below them.")
            print("  Expected for importance-weighted sources (GRAS treats event weights as frequency weights);")
            print("  use totalDose(..., replica_error=True) for the error.")

        if replica_error:
            print("Using the replica-scatter error; its own relative uncertainty is about",
                  f"{1 / np.sqrt(2 * (len(RawData) - 1)) * 100:.0f} %")
            TID['error'] = ReplicaError

    # Entries is the same for all tiles, because each tile gets an entry if a particle is detected in any tile.
    TotalEntries = TID['entries'][0]
    # Check if the number of entries is the same for all tiles.
    if not np.all(TID['entries'] == TID['entries'][0]):
        sys.exit("ERROR !!! Number of entries is not the same for all tiles")

    NumTiles = len(TID['dose'])
    MinNZE = min(TID['non-zeros'])   # Minimum number of non-zero entries
    LowestTile = np.argmin(TID['non-zeros'])
    ParticlesPerTile = TID['entries'][0] / NumTiles

    print("Total number of particles:", f"{TotalEntries:.3g}")
    print("Average number of Particles per tile:", f"{ParticlesPerTile:.3g}")
    print("Minimum number of non Zero Entries:", f"{MinNZE:.3g}", "at Tile number", str(LowestTile))


    #if MinNZE < 100:
    #    print("Warning !!! Tile", LowestTile, "has only", MinNZE, "Non-Zero entries !!!")

    if TID['dose'][LowestTile] == 0:
        print("Warning !!! Tile", LowestTile, "has zero dose !!!")
        RelativeError = np.zeros_like(TID['dose'])
        for i in range(NumTiles):
            if TID['dose'][i] != 0:
                RelativeError[i] = TID['error'][i] / TID['dose'][i]
            else:
                RelativeError[i] = 0
    else:
        RelativeError = TID['error'] / TID['dose']
    
    MaxRelativeError = max(RelativeError)
    MaxRelativeErrorTile = np.argmax(RelativeError)
    print("Maximum relative error =", f"{MaxRelativeError * 100:.2f} % at Tile number {MaxRelativeErrorTile}")

    PartNumRequired = (MaxRelativeError / 0.01) ** 2 * TotalEntries
    ParticleMultipleRequired = PartNumRequired / TotalEntries

    print(f"The number of particles required to achieve 1% error is: {PartNumRequired:.3g} or {ParticleMultipleRequired:.3g} times the number of particles in the run")

    if MaxRelativeError > 1:
        print(f"Warning !!! Tile {MaxRelativeErrorTile} has {MaxRelativeError * 100:.2f} % relative error!!!")

    print("Lowest Dose result:", f"{min(TID['dose']):.3g} +- {TID['error'][np.argmin(TID['dose'])]:.3g} kRad")

    # Print datatypes and shapes of the arrays
    #for key in keys:
    #    print(key, type(TID[key]), TID[key].shape)
        # print data type of the array entries
    #    print(key, TID[key].dtype)
        
    return TID


if __name__ == "__main__":
    Path = "/l/triton_work/2MatThickAl/OLD-GTO/AE9500keV"
    # Optional: python3 -m Dependencies.TotalDose <path> [--replica-error]
    # --replica-error takes the error from the replica scatter (biased GPS sources).
    Args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if Args:
        Path = Args[0]
    ReplicaError = "--replica-error" in sys.argv[1:]

    # Find all subdirectories in the given path that contain a "Res" subfolder
    # and calculate the total dose for each of them
    for root, dirs, files in os.walk(Path):
        if 'Res' in dirs:
            Path = os.path.join(root, 'Res/')
            print("Calculating total dose for path:", Path)

            # Extract folder name from root
            folder_name = os.path.basename(root)
            print("Folder name:", folder_name)
            
            # Calculate the total dose for the given path
            try:
                Results = totalDose(Path, replica_error=ReplicaError)
            except ValueError as error:
                print("Warning !!! Skipping", folder_name + ":", error)
                continue

            NumTiles = len(Results['dose'])


            # Print the dose results in csv format scientifically rounded and safe them to a csv file
            output_file = os.path.join(Path, "../TotalDose_" + folder_name + ".csv")
            with open(output_file, 'w') as f:
                Header = "Tile, Dose [kRad], Error [kRad], Non-Zero Entries"
                print("\n" + Header)
                f.write(Header + "\n")
            
                for i in range(NumTiles):
                    DoseWithError = ufloat(Results['dose'][i], Results['error'][i])
                    # Force uncertainties-style rounding (2 sig digits in uncertainty)
                    DoseWithErrorString = f"{DoseWithError:.2u}"            # rounded

                    # Read the rounded value back instead of splitting the string
                    # on '+/-'. At dose magnitudes like 1.4e-10 kRad the package
                    # emits a factored form, "(1.3975+/-0.0026)e-10", and a plain
                    # split yields "(1.3975", which float() rejects. Plotting/
                    # RadExTID.py already does this round trip; this keeps the two
                    # writers of TotalDose_<folder>.csv in agreement.
                    DoseRounded = ufloat_fromstr(DoseWithErrorString)
                    RoundedDose = DoseRounded.n
                    RoundedError = DoseRounded.s

                    print(f"{i}, {RoundedDose}, {RoundedError}, {Results['non-zeros'][i]}")
                    f.write(f"{i},{RoundedDose},{RoundedError},{Results['non-zeros'][i]}\n")


            # Plot the dose with error bars
            plt.figure(0)
            plt.errorbar(np.arange(NumTiles), Results['dose'], yerr=Results['error'], fmt=' ', capsize=5, elinewidth=1, capthick=1, label='Dose')
            # Add horizontal line at 1 kRad
            plt.axhline(y=1, color='r', linestyle='--', label='1 kRad')
            plt.title('Dose per tile ' + folder_name)
            plt.xlabel('Tile number')
            plt.ylabel('Dose [kRad]')
            plt.yscale('log')
            plt.grid(which='both')
            plt.legend()

            plt.savefig(Path + "../Dose_" + folder_name + ".pdf", format='pdf', bbox_inches="tight")

            
            # Plot the relative error
            plt.figure(1)
            plt.plot(100 * Results['error'] / Results['dose'], '.', label='Relative Error')
            # Add horizontal line at 1%
            plt.axhline(y=1, color='r', linestyle='--', label='1% error')
            plt.title('Relative Error in %')
            plt.xlabel('Tile number')
            plt.ylabel('Relative Error [%]')
            plt.ylim(bottom=0)
            plt.grid(which='both')
            plt.legend()

            plt.savefig(Path + "../Error_" + folder_name + ".pdf", format='pdf', bbox_inches="tight")

            # Plot the number of non-zero entries
            plt.figure(2)
            plt.plot(Results['non-zeros'], '.', label='Non-zero entries')
            # Add horizontal line at 1
            plt.axhline(y=1, color='r', linestyle='--', label='1 entry')
            plt.title('Number of non-zero entries')
            plt.xlabel('Tile number')
            plt.ylabel('Number of non-zero entries')
            plt.yscale('log')
            plt.grid(axis='x', which='both')
            plt.grid(axis='y', which='major')
            plt.legend()

            plt.savefig(Path + "../NonZeros_" + folder_name + ".pdf", format='pdf', bbox_inches="tight")
            
            # plt.show()
            plt.close('all')
