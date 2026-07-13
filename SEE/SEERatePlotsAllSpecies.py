import csv
import matplotlib.pyplot as plt
import numpy as np

Directories = {
    "CREME96": "/home/anton/triton_work/GRAS/LET_Histograms/Carrington-CREME96-All-Species/",
    "SAPPHIRE": "/home/anton/triton_work/GRAS/LET_Histograms/Carrington-SAPPHIRE-All-Species/",
}

CrossectionNames = ["Cypress CY62167GE30-45ZXI", "NanoXplore SEU"]

# Define default errorbar style
default_errorbar_style = {
    'capsize': 4,        # Size of the error bar caps
    'elinewidth': 2,     # Width of the error bar lines
    'capthick': 2,       # Thickness of the error bar caps
}

FigureNumber = 0
for ModelName, Directory in Directories.items():
  for CrossectionName in CrossectionNames:
    FigureNumber += 1

    # Initialize the final data structure
    Dict = {}

    # Read the CSV file
    with open(Directory + "SEERates_" + CrossectionName + ".csv", mode='r') as file:
        csv_reader = csv.DictReader(file)
        CrossectionName = None

        for row in csv_reader:
            # Extract relevant data
            dataset_name = row['Data']
            shielding = row['Shielding']
            see_rate = float(row['SEE_Rate'])
            see_error = float(row['SEE_Error'])
            entries_contributing_to_see = float(row['Entries_Contributing_To_SEE'])

            # Initialize the inner dictionary if it doesn't exist
            if dataset_name not in Dict:
                Dict[dataset_name] = {'Shielding': [], 'SEE_Rate': [], 'SEE_Error': [], 'Entries_Contributing_To_SEE': []}

            if CrossectionName:
                if CrossectionName != row['Crossection']:
                    raise ValueError("Crossection name is not the same for all entries")
            else:
                CrossectionName = row['Crossection']

            # Store the data dictionary
            Dict[dataset_name]['Shielding'].append(shielding)
            Dict[dataset_name]['SEE_Rate'].append(see_rate)
            Dict[dataset_name]['SEE_Error'].append(see_error)
            Dict[dataset_name]['Entries_Contributing_To_SEE'].append(entries_contributing_to_see)

    # Convert the lists to numpy arrays
    for key in Dict:
        for inner_key in Dict[key]:
            Dict[key][inner_key] = np.array(Dict[key][inner_key])

    plt.figure(FigureNumber, figsize=(5, 7))

    # Sort the datasets by their SEE rate at 4mm shielding, so the legend
    # is ordered from highest to lowest contribution
    Shielding4mmIndex = np.where(Dict[list(Dict.keys())[0]]['Shielding'] == '4mm')[0][0]
    sorted_names = sorted(Dict.keys(), key=lambda name: Dict[name]['SEE_Rate'][Shielding4mmIndex], reverse=True)
    Top3Names = sorted_names[:3]
    RestNames = sorted_names[3:]

    # Combine all species into a total rate
    Total = {}
    Total['Shielding'] = Dict[sorted_names[0]]['Shielding']
    Total['SEE_Rate'] = np.sum([Dict[name]['SEE_Rate'] for name in sorted_names], axis=0)
    Total['SEE_Error'] = np.sqrt(np.sum([np.square(Dict[name]['SEE_Error']) for name in sorted_names], axis=0))

    # Combine all species except the top 3 at 4mm shielding
    Total_Ex3 = {}
    Total_Ex3['Shielding'] = Dict[sorted_names[0]]['Shielding']
    Total_Ex3['SEE_Rate'] = np.sum([Dict[name]['SEE_Rate'] for name in RestNames], axis=0)
    Total_Ex3['SEE_Error'] = np.sqrt(np.sum([np.square(Dict[name]['SEE_Error']) for name in RestNames], axis=0))

    plt.errorbar(Total['Shielding'], Total['SEE_Rate'], yerr=Total['SEE_Error'],
                 label="Total", color='black', linewidth=3, **default_errorbar_style)

    for i, name in enumerate(Top3Names):
        label = name.replace("GEO-" + ModelName + "-Solar-", "")
        plt.errorbar(Dict[name]['Shielding'], Dict[name]['SEE_Rate'],
                     yerr=Dict[name]['SEE_Error'], label=label,
                     color='C' + str(i % 10), **default_errorbar_style)

    plt.errorbar(Total_Ex3['Shielding'], Total_Ex3['SEE_Rate'], yerr=Total_Ex3['SEE_Error'],
                 label="Total excluding top 3", color='C5', linestyle='--', **default_errorbar_style)

    plt.yscale("log")
    plt.grid()
    plt.title(CrossectionName + "\nGEO " + ModelName + " Solar Ion SEU Rate Estimates - All Species")

    plt.xlabel("Aluminium Shielding Thickness")
    plt.ylabel("Single Event Upset Rate [s-1 bit-1]")

    plt.legend(loc='lower left')

    plt.savefig(Directory + CrossectionName + " Rates All Species.pdf", format='pdf', bbox_inches="tight")
    # plt.show()
