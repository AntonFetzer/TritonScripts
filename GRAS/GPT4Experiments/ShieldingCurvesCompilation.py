from GRAS.Dependencies.TotalDose import totalkRadGras
import numpy as np
import matplotlib.pyplot as plt

Data = []
Labels = []

## A9 Trapped particles

A9_Paths = [
    "/l/triton_work/ShieldingCurves/MultilayerPaper/AE9-GTO/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/AP9-GTO/Res/"
]

Trapped_Data = []  # We'll use a separate list to keep the intermediate data.
for i, path in enumerate(A9_Paths):
    Trapped_Data.append(totalkRadGras(path, ""))

# Calculate total dose and error for Trapped particles
Total_Trapped_Dose = sum([data[0] for data in Trapped_Data])
Total_Trapped_Error = np.sqrt(sum([data[1]**2 for data in Trapped_Data]))
# Create a new dataset for Trapped with the same structure but with Entries and NonZeroEntries set to None
Total_Trapped_Data = [Total_Trapped_Dose, Total_Trapped_Error]
# Append this dataset to Data and the corresponding label to Labels
Data.append(Total_Trapped_Data)
Labels.append("Total Trapped Particles")


## Solar Energetic Particles

SEP_Paths = [
    "/l/triton_work/ShieldingCurves/MultilayerPaper/SolarGTO-Protons/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/SolarGTO-He/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/SolarGTO-C/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/SolarGTO-O/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/SolarGTO-Ne/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/SolarGTO-Mg/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/SolarGTO-Si/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/SolarGTO-Fe/Res/"
]

SEP_Data = []
for i, path in enumerate(SEP_Paths):
    SEP_Data.append(totalkRadGras(path, ""))

# Calculate total dose and error for SEP particles
Total_SEP_Dose = sum([data[0] for data in SEP_Data])
Total_SEP_Error = np.sqrt(sum([data[1]**2 for data in SEP_Data]))
# Create a new dataset for SEP with the same structure but with Entries and NonZeroEntries set to None
Total_SEP_Data = [Total_SEP_Dose, Total_SEP_Error]
# Append this dataset to Data and corresponding label to Labels
Data.append(Total_SEP_Data)
Labels.append("Solar Energetic Particles")

## Cosmic Particles

Cosmic_Paths = [
    "/l/triton_work/ShieldingCurves/MultilayerPaper/ISO-GTO-Protons-mission/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/ISO-GTO-He-mission/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/ISO-GTO-C-mission/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/ISO-GTO-N-mission/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/ISO-GTO-O-mission/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/ISO-GTO-Mg-mission/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/ISO-GTO-Si-mission/Res/",
    "/l/triton_work/ShieldingCurves/MultilayerPaper/ISO-GTO-Fe-mission/Res/"
]

Cosmic_Data = []
for i, path in enumerate(Cosmic_Paths):
    Cosmic_Data.append(totalkRadGras(path, ""))

seconds_in_a_month = 60 * 60 * 24 * 30.44  # number of seconds in a month

# Calculate total dose and error for Cosmic particles
Total_Cosmic_Dose = sum([data[0] for data in Cosmic_Data])/ seconds_in_a_month
Total_Cosmic_Error = np.sqrt(sum([data[1]**2 for data in Cosmic_Data])) / seconds_in_a_month
# Create a new dataset for Cosmic with the same structure
Total_Cosmic_Data = [Total_Cosmic_Dose, Total_Cosmic_Error]
# Append this dataset to Data and corresponding label to Labels
Data.append(Total_Cosmic_Data)
Labels.append("Cosmic Particles")


Colours = ['C0', 'C1', 'C2', 'C8', 'C3', 'C9', 'C7', 'k', 'C4', 'C5', 'C6',
           'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'white',
           '#800000', '#808000', '#800080', '#008080',  # some HTML hex colors
           (0.1, 0.2, 0.3), (0.4, 0.5, 0.6), (0.7, 0.8, 0.9),  # some RGB colors
           (0.1, 0.2, 0.3, 0.4), (0.5, 0.6, 0.7, 0.8)]  # some RGBA colors

plt.figure(figsize=(5, 4.2))
NumTiles = np.shape(Data[0])[1]
x = np.linspace(0, 2.5, num=NumTiles, endpoint=True)

for i, data in enumerate(Data):
    plt.errorbar(x, data[0], data[1], fmt='.', markersize=3, capsize=5, label=Labels[i], color=Colours[i], linestyle='')

####### Plot 10kRad line #########
CriticalDose = [1 for i in x]
plt.plot(x, CriticalDose, color='k', linewidth=2, label='1 krad per month')


plt.xlim(0.25, 2.5)
plt.ylim(1e-4, 2e2)
plt.grid(which='both', linestyle='-', linewidth=0.5)
plt.yscale("log")
plt.title("Ionising Dose per month behind shielding on GTO")
plt.xlabel("Aluminium Shielding Depth [g/cm2]")
plt.ylabel("Ionising Dose per month [krad]")
plt.legend()

plt.show()
#plt.savefig("/l/triton_work/ShieldingCurves/MultilayerPaper/ShieldingCurvesTypes.pdf", format='pdf', bbox_inches="tight")