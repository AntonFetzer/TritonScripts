from GRAS.Dependencies.TotalDose import totalkRadGras
import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat, unumpy

ProtElecOrTotal = "Total"  # "Protons", "Electrons", "Total

Proton_Path = "/l/triton_work/2LayerOptimisation10TilesPE/Protons/"
Electron_Path = "/l/triton_work/2LayerOptimisation10TilesPE/Electrons/"
res_suffix = "/Res/"

Labels = [
    "LH2",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "lN2",
    "lO2",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "lAr",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "lBr",
    "lKr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "lXe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf"
]

if ProtElecOrTotal == "Protons":
    Path = Proton_Path
    OutputPath = Proton_Path
elif ProtElecOrTotal == "Electrons":
    Path = Electron_Path
    OutputPath = Electron_Path
elif ProtElecOrTotal == "Total":
    Path = Electron_Path  # Set path first to Electron path. Proton path will be accessed later
    OutputPath = Electron_Path + "../"
else:
    print("ERROR: ProtElecOrTotal not set correctly")
    exit()

Paths = []
for i, label in enumerate(Labels, start=1):
    Paths.append(Path + str(i) + "-" + label + res_suffix)

Data = []
for i, path in enumerate(Paths):
    Data.append(totalkRadGras(path, ""))

if ProtElecOrTotal == "Total":
    ElectronData = Data

    ProtonPaths = []
    for i, label in enumerate(Labels, start=1):
        ProtonPaths.append(Proton_Path + str(i) + "-" + label + res_suffix)

    ProtonData = []
    for i, path in enumerate(ProtonPaths):
        ProtonData.append(totalkRadGras(path, ""))

    TotalData = []

    print(np.size(ElectronData))

    for i in range(len(ElectronData)):
        TotalDose = ElectronData[i][0] + ProtonData[i][0]
        TotalError = np.sqrt(ElectronData[i][1]**2 + ProtonData[i][1]**2)
        TotalEntries = ElectronData[i][2] + ProtonData[i][2]
        TotalNonZeroEntries = ElectronData[i][3] + ProtonData[i][3]

        TotalData.append(np.asarray([TotalDose, TotalError, TotalEntries, TotalNonZeroEntries]))

    Data = TotalData

Colours = ['C0', 'C1', 'C2', 'C8', 'C3', 'C9', 'C7', 'k', 'C4', 'C5', 'C6',
           'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'white',
           '#800000', '#808000', '#800080', '#008080',  # some HTML hex colors
           (0.1, 0.2, 0.3), (0.4, 0.5, 0.6), (0.7, 0.8, 0.9),  # some RGB colors
           (0.1, 0.2, 0.3, 0.4), (0.5, 0.6, 0.7, 0.8),  # some RGBA colors
           'C0', 'C1', 'C2', 'C8', 'C3', 'C9', 'C7', 'k', 'C4', 'C5', 'C6',
           'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'white']

Colours.extend(Colours)
Colours.extend(Colours)

NumTiles = np.shape(Data[0])[1]

x = np.linspace(0, (NumTiles - 1) * 10, num=NumTiles, dtype=int, endpoint=True)

for i, data in enumerate(Data):
    plt.errorbar(x, data[0], data[1], fmt='', markersize=5, capsize=5, label=Labels[i], color=Colours[i], linestyle='-')

####### Plot 1kRad line #########
CriticalDose = [1 for i in x]
plt.plot(x, CriticalDose, color='k', linewidth=2, label='1 krad per month')
CriticalDose = [10 for i in x]
plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='10 krad per month')

plt.ylim(0, 1)
plt.title("Ionising Dose per month behind multilayer shielding")
plt.xlabel("Percentage of shielding mass in Polyethylene [%]")
plt.ylabel("Ionising Dose per month [krad]")
plt.grid(which='both')
# plt.legend()

plt.savefig(OutputPath + "/2LayerOpt.pdf", format='pdf', bbox_inches="tight")
# plt.show()


# Data export ###############################################################

# Initialize lists to store additional calculated values
min_dose_idx = []
max_dose_idx = []

min_dose_data = []
max_dose_data = []

dose_at_0_data = []
dose_at_100_data = []

# For CSV File data processing
for dataset in Data:
    doses = dataset[0]
    errors = dataset[1]

    # Find min and max doses
    min_dose_idx.append(np.argmin(doses))
    max_dose_idx.append(np.argmax(doses))

    # Append the min and max doses with their errors to the lists
    min_dose_data.append((doses[min_dose_idx[-1]], errors[min_dose_idx[-1]]))
    max_dose_data.append((doses[max_dose_idx[-1]], errors[max_dose_idx[-1]]))

    dose_at_0_data.append((doses[0], errors[0]))
    dose_at_100_data.append((doses[-1], errors[-1]))

with open(OutputPath + 'PE-Materials.csv', 'w') as f:
    # Write the header to the CSV file
    header = ("Z,Material,Min Dose,Error,At PE %,Max Dose,Error,At PE %,Dose at 0% PE,Error,Dose at 100% PE,Error\n")
    f.write(header)

    # Loop over the Labels and Data to populate the CSV file
    for i, _ in enumerate(Data):
        # Convert the doses and errors to strings with ufloat for formatting
        min_dose_str = str(ufloat(*min_dose_data[i])).split('+/-')
        max_dose_str = str(ufloat(*max_dose_data[i])).split('+/-')
        dose_0_str = str(ufloat(*dose_at_0_data[i])).split('+/-')
        dose_100_str = str(ufloat(*dose_at_100_data[i])).split('+/-')

        line = (f"{i+1},"f"{Labels[i]},"
                f"{min_dose_str[0]},{min_dose_str[1]},{min_dose_idx[i] * 10},"
                f"{max_dose_str[0]},{max_dose_str[1]},{max_dose_idx[i] * 10},"
                f"{dose_0_str[0]},{dose_0_str[1]},"
                f"{dose_100_str[0]},{dose_100_str[1]}\n")
        f.write(line)
