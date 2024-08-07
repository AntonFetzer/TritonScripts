from GRAS.Dependencies.TotalDose import totalDose
import numpy as np
import matplotlib.pyplot as plt

Path = "/l/triton_work/2LayerOptimisationPhysicsTest10Tiles/Pe-Pb-Protons/"
res_suffix = "/Res/"


Labels = [
    "Template",
    ############# Models that all produce the same results as the default #################
    # "FTFP_BERT",
    # "QBBC",
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
    # "em_standard_opt3",
    "em_standard_opt4",  # Seems to be the default EM ???

    "em_penelope",  # Produces slightly higher dose
    "em_livermore",

    "em_lowenergy",  # Does not work for protons and huge deviations for electrons
    # "em_standard_space",  # Low performance
    # "em_standard_remizovich",

    # "em_standardSS",  # terrible performance
    "em_standardWVI",  # slow
    # "em_standardNR",  # terrible performance

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

Data = []
for i, path in enumerate(Paths):
    Data.append(totalDose(path))

keys = ['dose', 'error', 'entries', 'non-zeros']

# Compute the weighted average for the Dose.
total_entries = np.sum([data['entries'] for data in Data], axis=0)
weighted_doses = np.sum([data['dose'] * data['entries'] for data in Data], axis=0) / total_entries

# New Reference Dataset
new_reference = {
    keys[0]: weighted_doses,
    keys[1]: np.zeros_like(weighted_doses),
    keys[2]: np.zeros_like(weighted_doses),
    keys[3]: np.zeros_like(weighted_doses)
}

# Insert the new reference dataset at the beginning of the Data list
Data.insert(0, new_reference)

# Update Labels and Colors to include the new reference
Labels = ["Weighted Average Reference"] + Labels
Colours = ['grey'] + Colours

NumTiles = np.shape(Data[0]['dose'])[0]

plt.figure(1)

x = np.linspace(0, NumTiles - 1, num=NumTiles, dtype=int, endpoint=True)

for i, data in enumerate(Data):
    plt.errorbar(x, data['dose'], data['error'], fmt='', markersize=5, capsize=5, label=Labels[i], color=Colours[i], linestyle='-')

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

# Calculate the difference with the weighted average and plot it
for i, data in enumerate(Data[1:]):  # Skip the first dataset as it's the new reference
    DataDiff = data['dose'] - data['dose'][0] # Subtracting the reference dataset
    plt.plot(x, DataDiff, '-', label=Labels[i + 1], color=Colours[i + 1])  # i+1 as the first dataset is skipped

plt.title("Difference from the average")
plt.xlabel("Percentage of shielding mass in top layer [%]")
plt.ylabel("Difference in Ionising Dose per month [krad]")
plt.legend()
plt.grid(which='both')
plt.savefig(Path + "/Difference.pdf", format='pdf', bbox_inches="tight")

plt.figure(3)

# Calculate the relative deviation with the reference dataset (weighted average) and plot it
for i, data in enumerate(Data[1:]):  # Skip the first dataset as it's the new reference
    DataRatio = ((data['dose'] / data['dose'][0]) - 1) * 100  # Dividing by the reference dataset
    plt.plot(x, DataRatio, '-', label=Labels[i + 1], color=Colours[i + 1])  # i+1 as the first dataset is skipped

# plt.yscale("log")
plt.title("Relative Deviation from the average")
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
for i, data in enumerate(Data[1:]):  # Skip the first dataset as it's the new reference
    O = data['dose']  # Observed values (Dose for the current model)
    E = data['dose'][0]  # Expected values (Dose for the reference model)

    # Calculate the chi-squared value
    chi_squared = np.sum(((O - E) ** 2) / E)
    chi_squared_values.append(chi_squared)

    # Calculate the minimum number of non-zero entries
    min_non_zero_entries.append(np.min(data['non-zeros']))  # data['non-zeros'] is NonZeroEntries for the current model

    # Calculate the maximum relative error
    max_relative_error.append(np.max(data['error'] / data['dose']))

# Open a CSV file to write the data
with open(Path +'../Plotting/model_entries.csv', 'w') as f:
    # Write the header to the CSV file
    f.write("Model Name,Number of Entries,Chi-Squared,Min Non-Zero Entries,Max Relative Error\n")

    # Loop over the Labels and Data to populate the CSV file
    for i, data in enumerate(Data[1:]):  # Skip the first dataset as it's the reference
        f.write(
            f"{Labels[i + 1]},{int(data['entries'][0])},{chi_squared_values[i]},{min_non_zero_entries[i]},{max_relative_error[i]}\n")  # i+1 as the first dataset is skipped
