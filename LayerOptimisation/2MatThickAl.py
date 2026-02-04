import os
import numpy as np
from Dependencies.TotalDose import totalDose
from Read.ReadSD2Q import readSDQ2
from Dependencies.MergeTotalDose import mergeTotalDose
import matplotlib.pyplot as plt
from uncertainties import ufloat, ufloat_fromstr
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.ticker import FuncFormatter


Path = "/l/triton_work/2MatThickAl/Apophis-1y-GTO/"

ElectronFolder = "Apophis-AE9-95p-1y"
ProtonFolder = "Apophis-AP9-95p-1y"
SolarFolder = "Apophis-SAPPHIRE-95p-1y"
TotalFolder = "Apophis-A9-total-95p-1y-mean"

Datset = "P" # Electron, P - Proton, S - Solar, T - Total

print("Calculating total dose for Folder name:", Path)
    
# Calculate the total dose for the given path and dataset
if Datset == "E":
    Results = totalDose(os.path.join(Path, ElectronFolder, "Res"))
    folder = ElectronFolder
elif Datset == "P":
    Results = totalDose(os.path.join(Path, ProtonFolder, "Res"))
    folder = ProtonFolder
elif Datset == "S":
    Results = totalDose(os.path.join(Path, SolarFolder, "Res"))
    folder = SolarFolder
elif Datset == "T":
    ElectronResults = totalDose(os.path.join(Path, ElectronFolder, "Res"))
    ProtonResults = totalDose(os.path.join(Path, ProtonFolder, "Res"))
    SolarResults = totalDose(os.path.join(Path, SolarFolder, "Res"))
    Results = mergeTotalDose([ElectronResults, ProtonResults, SolarResults])
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

for x in range(TilesX):
    for y in range(TilesY):
        Map[y, x] = Results['dose'][x + y * TilesX]

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
sd2q_file = "/l/triton_work/Spectra/Apophis/spenvis_sqo.txt"
sd2q_data = readSDQ2(sd2q_file)
#plt.plot(sd2q_data['Thickness'] * 0.27, (sd2q_data['Electrons'] + sd2q_data['Bremsstrahlung'])/ (11*12), '.-', label="Shieldose Electrons")
if Datset == "E":
    plt.plot(sd2q_data['Thickness'] * 0.27, sd2q_data['Electrons'], '.-', label="Shieldose Electrons Only")
    plt.plot(sd2q_data['Thickness'] * 0.27, sd2q_data['Bremsstrahlung'], '.-', label="Shieldose Bremsstrahlung")
if Datset == "P" or Datset == "S":
    plt.plot(sd2q_data['Thickness']/2.699, sd2q_data['Trapped Protons'], '.-', label="Shieldose Protons")
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


# Normalize each row to its last entry (reference material Al)
MapNormalized = np.zeros_like(Map)
for row in range(TilesY):
    reference_value = Map[row, -1]  # Last entry is Al reference
    if reference_value != 0:
        MapNormalized[row, :] = Map[row, :] / reference_value
    else:
        MapNormalized[row, :] = Map[row, :]

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

# Set the normalization to center at 1.0 (100% reference)
norm = TwoSlopeNorm(vmin=np.min(MapNormalized), vcenter=1.0, vmax=np.max(MapNormalized))

plt.figure(figsize=(8, 6))
im = plt.imshow(MapNormalized, cmap=cmap_diverging, aspect='auto', norm=norm)
cbar = plt.colorbar(im, label='Dose relative to same mass of Al (%)')
cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{int(x*100)}%'))

# Mark the lowest entry with an X
min_idx = np.unravel_index(np.argmin(MapNormalized), MapNormalized.shape)
plt.plot(min_idx[1], min_idx[0], 'x', color='black', markersize=12, markeredgewidth=2)

# Set x-axis labels to show PE percentages and Al
x_labels = [f'{i}% PE' for i in range(0, 110, 10)] + ['Al']
plt.xticks(range(TilesX), x_labels, rotation=45, ha='right')

# Set y-axis labels to show shielding depth in g/cm²
y_ticks = np.linspace(0, TilesY-1, 10, dtype=int)
y_labels = [f'{shielding_depth[i]:.2f}' for i in y_ticks]
plt.yticks(y_ticks, y_labels)

plt.title('Dose Heatmap ' + folder)
plt.xlabel('Shielding Material (% PE)')
plt.ylabel('Shielding depth [g/cm²]')
plt.savefig(Path + "/" + folder + "/Plot/DoseHeatmap_" + folder + ".pdf", format='pdf', bbox_inches="tight")
plt.close()
