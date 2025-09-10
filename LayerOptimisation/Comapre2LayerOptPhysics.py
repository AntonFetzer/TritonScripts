from GRAS.Dependencies.TotalDose import totalDose
import numpy as np
import matplotlib.pyplot as plt

Path = "/l/triton_work/2LayerOptimisationPhysicsTest10Tiles/Pe-Pb-Electrons/"
res_suffix = "/Res/"


Labels = [
    "Template",
    ############# Models that all produce the same results as the default #################
    "FTFP_BERT",
    "QBBC",
    # "QGSP_BERT",
    # "QGSP_BIC",
    # "binary_ion",
    # "qmd_ion",
    # "incl_ion",
    # "decay",
    # "raddecay",
    # "gamma_nuc",
    # "stopping",
    # "binary_had",
    # "bertini",
    # "bertini_preco",
    # "elastic",
    # "elasticCHIPS",
    # "binary",
    # "firsov",

    "EM3andQBBC",
    "QBBCandEM3",

    ############# Models that dont work at all #################
    # "QGSP_BIC_HP",
    # "QGSP_QMD_HP",
    # "Shielding",
    # "bertini_hp",
    # "binary_hp",
    # "elasticHP",



    ############# Models that produce significantly different results as the default #################
    # "em_standard",  # Not actually the standard
    # "em_standard_opt1",  # low precision
    # "em_standard_opt2",  # experimental
    "em_standard_opt3",
    "em_standard_opt4",  # Seems to be the default EM ???

    # "em_penelope",  # Produces slightly higher dose
    # "em_livermore",

    # "em_lowenergy",
    # "em_standard_space",  # Low performance and significant deviations
    # "em_standard_remizovich", # No Proton results

    # "em_standardSS",  # terrible performance No Proton results
    # "em_standardWVI",  # slow
    # "em_standardNR",  # terrible performance No Proton results

    # Does not seem to work at all
    # "secondary_generator",
    # "rmc_em_standard",  #No Proton results. This is reverse Monte Carlo. DO NOT USE !!!
    # "QGSP_BERT_HP", 
    # 
]


Paths = [Path + label + res_suffix for label in Labels]


Colours = ['C0', 'C1', 'C2', 'C8', 'C3', 'C9', 'C7', 'k', 'C4', 'C5', 'C6',
           'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'white',
           '#800000', '#808000', '#800080', '#008080',  # some HTML hex colors
           (0.1, 0.2, 0.3), (0.4, 0.5, 0.6), (0.7, 0.8, 0.9),  # some RGB colors
           (0.1, 0.2, 0.3, 0.4), (0.5, 0.6, 0.7, 0.8),  # some RGBA colors
           'C0', 'C1', 'C2', 'C8', 'C3', 'C9', 'C7', 'k', 'C4', 'C5', 'C6',
           'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'white']

List_of_Datasets = []
i = 0
while i < len(Paths):
    try:
        List_of_Datasets.append(totalDose(Paths[i]))
        i += 1
    except:
        print(f"Failed to read in data from {Paths[i]}")
        # Remove the label and colour for the failed dataset
        Paths.pop(i)
        Labels.pop(i)

keys = ['dose', 'error', 'entries', 'non-zeros']

NumTiles = len(List_of_Datasets[0]['dose'])

plt.figure(1)

x = np.linspace(0, NumTiles - 1, num=NumTiles, dtype=int, endpoint=True)

for i, TID in enumerate(List_of_Datasets):
    plt.errorbar(x, TID['dose'], TID['error'], fmt='', markersize=5, capsize=5, label=Labels[i], color=Colours[i], linestyle='-')

####### Plot 1kRad line #########
# CriticalDose = [1 for i in x]
# plt.plot(x, CriticalDose, color='k', linewidth=2, label='1 krad per month')
# CriticalDose = [10 for i in x]
# plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='10 krad per month')

# plt.ylim(0, 1)
# plt.yscale("log")
# plt.xscale("log")
plt.title("Ionising Dose per month behind multilayer shielding")
plt.xlabel("Percentage of shielding mass in top layer [%]")
plt.ylabel("Ionising Dose per month [krad]")
plt.grid(which='both')
plt.legend()

plt.savefig(Path + "/2LayerOpt.pdf", format='pdf', bbox_inches="tight")
# plt.show()

plt.figure(2)

# Calculate the difference from the first dataset and plot it
for i, TID in enumerate(List_of_Datasets[1:]):  # Skip the first dataset as it's the reference
    TIDDiff = TID['dose'] - List_of_Datasets[0]['dose'] # Subtracting the TID of the first dataset
    plt.plot(x, TIDDiff, '-', label=Labels[i + 1], color=Colours[i + 1])  # i+1 as the first dataset is skipped

plt.title("Difference from the first dataset")
plt.xlabel("Percentage of shielding mass in top layer [%]")
plt.ylabel("Difference in Ionising Dose per month [krad]")
plt.legend()
plt.grid(which='both')
plt.savefig(Path + "/Difference.pdf", format='pdf', bbox_inches="tight")

plt.figure(3)

# Calculate the relative deviation with the first dataset and plot it
for i, TID in enumerate(List_of_Datasets[1:]):  # Skip the first dataset as it's the reference
    TIDRatio = ((TID['dose'] / List_of_Datasets[0]['dose']) - 1) * 100  # Dividing by the reference dataset
    plt.plot(x, TIDRatio, '-', label=Labels[i + 1], color=Colours[i + 1])  # i+1 as the first dataset is skipped

# plt.yscale("log")
plt.title("Relative Deviation from the first dataset")
plt.xlabel("Percentage of shielding mass in top layer [%]")
plt.ylabel("Relative Deviation [%]")
plt.legend()
plt.grid()
plt.savefig(Path + "/RelativeDeviation.pdf", format='pdf', bbox_inches="tight")

# Initialize lists to store additional calculated values
chi_squared_values = []
min_non_zero_entries = []
max_relative_error = []

# For CSV File creation and writing
for i, TID in enumerate(List_of_Datasets[1:]):  # Skip the first List_of_Datasetsset as it's the reference
    O = TID['dose']  # Observed values (Dose for the current model)
    E = List_of_Datasets[0]['dose']  # Expected values (Dose for the reference model)

    # Calculate the chi-squared value
    chi_squared = np.sum(((O - E) ** 2) / E)
    chi_squared_values.append(chi_squared)

    # Calculate the minimum number of non-zero entries
    min_non_zero_entries.append(np.min(TID['non-zeros']))

    # Calculate the maximum relative error
    max_relative_error.append(np.max(TID['error'] / TID['dose']))

# Open a CSV file to write the List_of_Datasets
with open(Path +'../Plotting/model_entries.csv', 'w') as f:
    # Write the header to the CSV file
    f.write("Model Name,Number of Entries,Chi-Squared,Min Non-Zero Entries,Max Relative Error\n")

    # Write the reference model to the CSV file
    f.write(f"{Labels[0]},{int(sum(List_of_Datasets[0]['entries']))},0,{np.min(List_of_Datasets[0]['non-zeros'])},{np.max(List_of_Datasets[0]['error'] / List_of_Datasets[0]['dose'])}\n")

    # Loop over the Labels and List_of_Datasets to populate the CSV file
    for i, TID in enumerate(List_of_Datasets[1:]):  # Skip the first List_of_Datasetsset as it's the reference
        f.write(f"{Labels[i+1]},{int(sum(TID['entries']))},{chi_squared_values[i]},{min_non_zero_entries[i]},{max_relative_error[i]}\n") 
        # i+1 as the first List_of_Datasetsset is skipped
