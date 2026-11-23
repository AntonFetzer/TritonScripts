import sys, os, re
# Ensure the package root (the "Python" folder) is on sys.path so "from Read.* import ..." works
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import numpy as np
from Dependencies.TotalDose import totalDose
from Read.ReadSD2Q import readSDQ2
from Dependencies.MergeTotalDose import mergeTotalDose
from Dependencies.CombineMeasurements import ratioWithError
import matplotlib.pyplot as plt
from uncertainties import ufloat, ufloat_fromstr
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.ticker import FuncFormatter, MaxNLocator


Path = "/scratch/work/fetzera1/GRAS/2D/2MatThickAl/MEO-11y-mean"

ElectronFolder = "MEO-AE9-11y-mean"
ProtonFolder = "MEO-AP9-11y-mean"  # None when the study has no trapped-proton run
SolarFolder = "MEO-SAPPHIRE-11y-mean"  # None when the study has no solar-proton run
TotalFolder = "MEO-total-11y-mean"

Datset = "T" #E - Electron, P - Proton, S - Solar, T - Total

sd2q_file = "/scratch/work/fetzera1/Spectra/Carrington/MEO/spenvis_sqo.txt"
HeatmapTitle = "Dose relative to same mass of Al on MEO (11-year mean)"
# Factor from the run's dose to kRad per month: 30*24*3600 for per-second spectra (OLD-GTO), 1/(11*12) for 11-year fluences, 1/(5*12) for 5-year fluences
MonthFactor = 1/(11*12)
# Biased electron source (flat-log GPS bias): GRAS 5 overestimates the quoted errors, so take the electron error from the replica scatter
ReplicaError = True

print("Calculating total dose for Folder name:", Path)
    
# Calculate the total dose for the given path and dataset
if Datset == "E":
    Results = totalDose(os.path.join(Path, ElectronFolder, "Res"), replica_error=ReplicaError)
    folder = ElectronFolder
elif Datset == "P":
    Results = totalDose(os.path.join(Path, ProtonFolder, "Res"))
    folder = ProtonFolder
elif Datset == "S":
    Results = totalDose(os.path.join(Path, SolarFolder, "Res"))
    folder = SolarFolder
elif Datset == "T":
    ElectronResults = totalDose(os.path.join(Path, ElectronFolder, "Res"), replica_error=ReplicaError)
    Components = [ElectronResults]
    if ProtonFolder:
        Components.append(totalDose(os.path.join(Path, ProtonFolder, "Res")))
    if SolarFolder:
        Components.append(totalDose(os.path.join(Path, SolarFolder, "Res")))
    Results = mergeTotalDose(Components)
    folder = TotalFolder

NumTiles = len(Results['dose'])


# Convert from flux to 1 month fluence
# Results['dose'] *= 30*24*3600  # seconds in 1 month
# Results['error'] *= 30*24*3600

# Convert from 11 year mean to 1 month fluence
#Results['dose'] /= (11*12)
#Results['error'] /= (11*12)


# Print the dose results in csv format scientifically rounded and safe them to a csv file
output_file = os.path.join(Path, folder, "Plot/TotalDose_" + folder + ".csv")
with open(output_file, 'w') as f:
    Header = "Tile, Dose [kRad], Error [kRad], Non-Zero Entries"
    #print("\n" + Header)
    f.write(Header + "\n")

    for i in range(NumTiles):
        DoseWithError = ufloat(Results['dose'][i], Results['error'][i])

        # This performs uncertainties-style rounding to match the display (2 sig digits on uncertainty)
        s = f"{DoseWithError:.2u}"

        # Re-create a ufloat from the rounded string, then extract numeric parts safely
        DoseRounded = ufloat_fromstr(s)

        RoundedDose = DoseRounded.n
        RoundedError = DoseRounded.s

        RoundedDoseString = f"{DoseRounded.n:.2g}"
        RoundedErrorString = f"{DoseRounded.s:.2g}"

        #print(f"{i}, {RoundedDoseString}, {RoundedErrorString}, {Results['non-zeros'][i]}")
        f.write(f"{i},{RoundedDoseString},{RoundedErrorString},{Results['non-zeros'][i]}\n")


# Plot the dose with error bars
plt.figure(0)
plt.errorbar(np.arange(NumTiles), Results['dose'], yerr=Results['error'], fmt=' ', capsize=5, elinewidth=1, capthick=1, label='Dose')
# Add horizontal line at 1 kRad
plt.axhline(y=1, color='r', linestyle='--', label='1 kRad')
plt.title('Dose per tile ' + folder)
plt.xlabel('Tile number')
plt.ylabel('Dose [kRad]')
plt.yscale('log')
plt.grid(which='both')
plt.legend()

plt.savefig(Path + "/" + folder + "/Plot/Dose_" + folder + ".pdf", format='pdf', bbox_inches="tight")


# Plot the relative error
plt.figure(1)
plt.plot(100 * Results['error'] / Results['dose'], '.', label='Relative Error')
# Add horizontal line at 1%
plt.axhline(y=1, color='r', linestyle='--', label='1% error')
plt.title('Relative Error in %')
plt.xlabel('Tile number')
plt.ylabel('Relative Error [%]')
plt.grid(which='both')
plt.legend()

plt.savefig(Path + "/" + folder + "/Plot/Error_" + folder + ".pdf", format='pdf', bbox_inches="tight")

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

plt.savefig(Path + "/" + folder + "/Plot/NonZeros_" + folder + ".pdf", format='pdf', bbox_inches="tight")



# Plot heat maps for each tile
# The tiles are arranges in 12 by 32 grid
TilesX = 12
TilesY = 32

Map = np.zeros((TilesY, TilesX))
ErrorMap = np.zeros((TilesY, TilesX))

for x in range(TilesX):
    for y in range(TilesY):
        Map[y, x] = Results['dose'][x + y * TilesX]
        ErrorMap[y, x] = Results['error'][x + y * TilesX]

# Save map as CSV file with header and thickness column
map_csv_file = os.path.join(Path, folder, "Plot/DoseHeatmap_" + folder + ".csv")
# Create header with material percentages
header = "Shielding Depth [g/cm2]," + ",".join([f"{i}% PE" if i < 110 else "Al" for i in range(0, 120, 10)]) + ",Al Thickness [mm]"
# Create thickness column (0.5 mm increasing by 10% per tile in Y)
thicknesses = np.array([0.5 * (1.1 ** y) for y in range(TilesY)])  # in mm
# Convert thickness to shielding depth in g/cm² assuming aluminum density (2.7 g/cm³)
al_density = 2.699  # g/cm³
shielding_depth = thicknesses * al_density / 10  # convert mm to cm, then multiply by density
# Combine thickness column, shielding depth column, and Map
MapWithThickness = np.column_stack([shielding_depth, Map, thicknesses])
# Save to CSV
with open(map_csv_file, 'w') as f:
    f.write(header + '\n')
    np.savetxt(f, MapWithThickness, delimiter=',', fmt='%.6e')


# Plot shielding curves for the 12 different material compositions. Dose vs shielding depth
plt.figure(figsize=(10, 6))
for x in [11]:  # 0%, 50%, 100% PE and Al]:
    doses = Map[:, x]
    plt.plot(shielding_depth, doses, marker='.', label=f'{x*10}% PE' if x < 11 else 'Al')

# Add shielddose data for comparison
sd2q_data = readSDQ2(sd2q_file)
#plt.plot(sd2q_data['Thickness'] * 0.27, (sd2q_data['Electrons'] + sd2q_data['Bremsstrahlung'])/ (11*12), '.-', label="Shieldose Electrons")
if Datset == "E":
    plt.plot(sd2q_data['Thickness'] * 0.27, sd2q_data['Electrons'], '.-', label="Shieldose Electrons Only")
    plt.plot(sd2q_data['Thickness'] * 0.27, sd2q_data['Bremsstrahlung'], '.-', label="Shieldose Bremsstrahlung")
if Datset == "P" or Datset == "S":
    plt.plot(sd2q_data['Thickness']/2.699, sd2q_data['Trapped Protons'], '.-', label="Shieldose Protons")
    if len(sd2q_data['Solar Protons']):  # some SHIELDOSE-2Q files (OLD-GTO) have no solar protons
        plt.plot(sd2q_data['Thickness'] * 0.27, sd2q_data['Solar Protons'], '.-', label="Shieldose Solar Protons")
if Datset == "T":
    plt.plot(sd2q_data['Thickness'] * 0.27, sd2q_data['Total'], '.-', label="Shieldose Total Dose")

plt.title('Dose vs Shielding Depth ' + folder)
plt.xlabel('Shielding depth [g/cm²]')
plt.ylabel('Dose [kRad]')
plt.yscale('log')
#plt.xscale('log')
plt.xlim(0, 2.7)
#plt.ylim(1e-2, 2e2)
plt.grid(which='both')
plt.legend(title='Material Composition', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig(Path + "/" + folder + "/Plot/DoseShieldingCurves_" + folder + ".pdf", format='pdf', bbox_inches="tight")
plt.close()


# Create a custom diverging colormap: green-white-red
# Green for values < 1 (lower than reference), red for values > 1 (higher than reference)
colors_low = ['#228B22', '#FFFFFF']  # Green to white
colors_high = ['#FFFFFF', '#FF0000']  # White to red
n_bins = 256
cmap_low = LinearSegmentedColormap.from_list('green_white', colors_low)
cmap_high = LinearSegmentedColormap.from_list('white_red', colors_high)

# Combine the two colormaps
colors_combined = list(cmap_low(np.linspace(0, 1, n_bins//2))) + list(cmap_high(np.linspace(0, 1, n_bins//2)))
cmap_diverging = LinearSegmentedColormap.from_list('green_white_red', colors_combined, N=n_bins)

# One heat map per reference: Al (last column), 100% PE (column 10), 100% Pb (column 0, 0% PE),
# and the lowest of these three single materials at each depth, which shows the true multilayer advantage
for RefColumns, RefName, FileSuffix in [([-1], "Al", ""), ([10], "PE", "PE"), ([0], "Pb", "Pb"),
                                        ([0, 10, -1], "best single material", "Best")]:
    # Normalize each row to the reference column (the lowest one when there are several)
    MapNormalized = np.zeros_like(Map)
    for row in range(TilesY):
        reference_value = np.min(Map[row, RefColumns])
        if reference_value != 0:
            MapNormalized[row, :] = Map[row, :] / reference_value
        else:
            MapNormalized[row, :] = Map[row, :]

    # Set the normalization to center at 1.0 (100% reference); extend the scale just past 1 when all values lie on one side
    norm = TwoSlopeNorm(vmin=min(np.min(MapNormalized), 0.99), vcenter=1.0, vmax=max(np.max(MapNormalized), 1.01))

    plt.figure(figsize=(8, 6))
    im = plt.imshow(MapNormalized, cmap=cmap_diverging, aspect='auto', norm=norm)
    # Show the dose change against the reference: a dose of 60 % of the reference reads -40 %
    cbar = plt.colorbar(im, label='Dose change relative to same mass of ' + RefName + ' (%)')
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{(x - 1) * 100:+.0f}%' if abs(x - 1) >= 0.005 else '0%'))
    # Ticks for the green (reduction) and red (increase) halves separately, so both get a similar density
    Low, High = (norm.vmin - 1) * 100, (norm.vmax - 1) * 100
    Change = np.union1d(MaxNLocator(nbins=6, integer=True).tick_values(Low, 0), MaxNLocator(nbins=6, integer=True).tick_values(0, High))
    cbar.set_ticks([1 + c / 100 for c in Change if Low <= c <= High])

    # Mark the lowest entry with an X
    min_idx = np.unravel_index(np.argmin(MapNormalized), MapNormalized.shape)
    plt.plot(min_idx[1], min_idx[0], 'x', color='black', markersize=12, markeredgewidth=2)

    # Write the best multilayer point relative to the best single material, with its dose per month
    if FileSuffix == "Best":
        row, col = min_idx
        RefColumn = [c % TilesX for c in RefColumns][np.argmin(Map[row, RefColumns])]
        if col == RefColumn:
            Ratio, RatioError = 1.0, 0.0
        else:
            Ratio, RatioError = ratioWithError(Map[row, col], ErrorMap[row, col], Map[row, RefColumn], ErrorMap[row, RefColumn])
        MaterialLabels = [f'{i}% PE' for i in range(0, 110, 10)] + ['Al']
        RefLabels = {0: "Pb", 10: "PE", TilesX - 1: "Al"}
        with open(os.path.join(Path, folder, "Plot/BestMultilayer_" + folder + ".csv"), 'w') as f:
            f.write("Folder,Relative Dose [%],Relative Dose Error [%],Shielding Depth [g/cm2],Material,"
                    "TID [kRad per month],TID Error [kRad per month],Best Single Material\n")
            f.write(f"{folder},{100 * Ratio:.4g},{100 * RatioError:.2g},{shielding_depth[row]:.4g},{MaterialLabels[col]},"
                    f"{Map[row, col] * MonthFactor:.4g},{ErrorMap[row, col] * MonthFactor:.2g},{RefLabels[RefColumn]}\n")

    # Label the minimum next to the X; place it on the left when the X is in the right half
    MinLabel = f"{(MapNormalized[min_idx] - 1) * 100:+.0f}%" if abs(MapNormalized[min_idx] - 1) >= 0.005 else "0%"
    OnRight = min_idx[1] < TilesX / 2
    plt.annotate(MinLabel, (min_idx[1], min_idx[0]), xytext=(12 if OnRight else -12, 0), textcoords='offset points',
                 ha='left' if OnRight else 'right', va='center',
                 bbox=dict(boxstyle='round', facecolor='white', edgecolor='black', alpha=0.9))

    # Set x-axis labels to show PE percentages and Al
    x_labels = [f'{i}% PE' for i in range(0, 110, 10)] + ['Al']
    plt.xticks(range(TilesX), x_labels, rotation=45, ha='right')

    # Set y-axis labels to show shielding depth in g/cm²
    y_ticks = np.linspace(0, TilesY-1, 10, dtype=int)
    y_labels = [f'{shielding_depth[i]:.2f}' for i in y_ticks]
    plt.yticks(y_ticks, y_labels)

    Title = HeatmapTitle.replace("same mass of Al", "same mass of " + RefName) if HeatmapTitle else 'Dose Heatmap ' + RefName + ' ' + folder
    if FileSuffix == "Best":
        Title = re.sub(r"\s*\([^)]*\)", "", Title)  # the longest title: drop the bracketed spectrum details
    plt.title(Title)
    plt.xlabel('Shielding Material (% PE)')
    plt.ylabel('Shielding depth [g/cm²]')
    plt.savefig(Path + "/" + folder + "/Plot/DoseHeatmap" + FileSuffix + "_" + folder + ".pdf", format='pdf', bbox_inches="tight")
    plt.close()

# Heat map of the absolute dose per month of every configuration, centred at 1 kRad per month
MapMonth = Map * MonthFactor
# Logarithmic colour scale: colour log10 of the dose, so 1 kRad per month sits at 0
LogMapMonth = np.log10(MapMonth)
# When all doses lie on one side of 1 kRad per month, extend the scale just past 1 so the centre stays inside it
norm = TwoSlopeNorm(vmin=min(np.min(LogMapMonth), -0.01), vcenter=0.0, vmax=max(np.max(LogMapMonth), 0.01))
plt.figure(figsize=(8, 6))
im = plt.imshow(LogMapMonth, cmap=cmap_diverging, aspect='auto', norm=norm)
cbar = plt.colorbar(im, label='Dose [kRad per month]')
cbar.set_ticks(np.log10([m * 10.0**e for e in range(-4, 4) for m in (1, 2, 5)
                         if norm.vmin <= np.log10(m * 10.0**e) <= norm.vmax]))
cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{10**x:.4g}'))
cbar.ax.yaxis.set_minor_locator(plt.NullLocator())
plt.xticks(range(TilesX), x_labels, rotation=45, ha='right')
plt.yticks(y_ticks, y_labels)
plt.title(HeatmapTitle.replace("Dose relative to same mass of Al", "Absolute dose") if HeatmapTitle else 'Absolute Dose Heatmap ' + folder)
plt.xlabel('Shielding Material (% PE)')
plt.ylabel('Shielding depth [g/cm²]')
plt.savefig(Path + "/" + folder + "/Plot/DoseHeatmapAbsolute_" + folder + ".pdf", format='pdf', bbox_inches="tight")
plt.close()
