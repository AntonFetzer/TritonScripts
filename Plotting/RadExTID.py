import os
import argparse
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from Dependencies.AggregateRun import aggregateRun
import matplotlib.pyplot as plt
from uncertainties import ufloat, ufloat_fromstr

parser = argparse.ArgumentParser(description="Aggregate and plot RadEx TID results.")
parser.add_argument("path", nargs="?", default="/l/triton_work/RadEx/RadEx-ThickPCB-6mm")
parser.add_argument("--fluence", type=float, help="Incident fluence in particles/cm2; overrides folder-name convention")
args = parser.parse_args()
if args.fluence is not None and (not np.isfinite(args.fluence) or args.fluence <= 0):
    parser.error("--fluence must be finite and positive")
Path = args.path

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
        Results = aggregateRun(Path)

        NumTiles = len(Results['dose'])

        fluence = args.fluence if args.fluence is not None else (2e12 if 'MeVElectron' in folder_name else 1.0)
        Results['dose'] *= fluence
        Results['error'] *= fluence
        print(f"Applied incident fluence: {fluence:g} particles/cm2")


        # Print the dose results in csv format scientifically rounded and safe them to a csv file
        output_file = os.path.join(Path, "../TotalDose_" + folder_name + ".csv")
        with open(output_file, 'w') as f:
            Header = "Tile, Dose [kRad], Error [kRad], Non-Zero Entries"
            print("\n" + Header)
            f.write(Header + "\n")
        
            for i in range(NumTiles):
                DoseWithError = ufloat(Results['dose'][i], Results['error'][i])

                # This performs uncertainties-style rounding to match the display (2 sig digits on uncertainty)
                s = f"{DoseWithError:.2u}"

                # Re-create a ufloat from the rounded string, then extract numeric parts safely
                DoseRounded = ufloat_fromstr(s)

                RoundedDose = DoseRounded.n
                RoundedError = DoseRounded.s

                RoundedDoseString = f"{DoseRounded.n:.12g}"
                RoundedErrorString = f"{DoseRounded.s:.12g}"

                print(f"{i}, {RoundedDoseString}, {RoundedErrorString}, {Results['non-zeros'][i]}")
                f.write(f"{i},{RoundedDoseString},{RoundedErrorString},{Results['non-zeros'][i]}\n")


        # Plot the dose with error bars
        plt.figure(0)
        plt.errorbar(np.arange(NumTiles), Results['dose'], yerr=Results['error'], fmt=' ', capsize=5, elinewidth=1, capthick=1, label='Dose')
        # Add horizontal line at 1 kRad
        plt.axhline(y=1, color='r', linestyle='--', label='1 kRad')
        plt.title('Dose per tile ' + folder_name, fontsize=10)
        plt.xlabel('Tile number')
        plt.ylabel('Dose [kRad]')
        plt.yscale('log')
        plt.grid(which='both')
        plt.legend()

        plt.savefig(Path + "../Dose_" + folder_name + ".pdf", format='pdf', bbox_inches="tight")

        
        # Plot the relative error
        plt.figure(1)
        plt.plot(Results['relative_error_percent'], '.', label='Relative Error')
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
