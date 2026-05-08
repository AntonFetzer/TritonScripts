import csv
import matplotlib.pyplot as plt
import numpy as np

Directory = "/l/triton_work/Fluence_Histograms/CarringtonShielded/"

CrossectionName = "NanoXplore Proton SEU"  
# CrossectionName = "Cypress CY62167GE30-45ZXI Proton SEU"

# Define default errorbar style
default_errorbar_style = {
    'capsize': 4,        # Size of the error bar caps
    'elinewidth': 2,     # Width of the error bar lines
    'capthick': 2,       # Thickness of the error bar caps
}

Expected = 'blue' # Blue
PlusColor = 'C1'  # Green
MinusColor = 'C2' # Orange
LEOColor = 'C8'   # yellow
GEOColor = 'C7'   # grey
VAPColor = 'C3'   # Red
CREME96Color = 'C6' # Pink
SAP_PeakColor = 'C9' # Turquoise

# Initialize the final data structure
Dict = {}

# Read the CSV file
with open(Directory + "/SEERates_" + CrossectionName + ".csv", mode='r') as file:
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


#,linestyle='', capsize=1.5, capthick=2, elinewidth=1

plt.figure(0, figsize=(5, 7))
########################### Carrington SEP ###################################
Shielding = Dict['Carrington-SEP-Expected-Int']['Shielding']
ExpectedInt = Dict['Carrington-SEP-Expected-Int']['SEE_Rate']
Minus = Dict['Carrington-SEP-Minus2Sigma-Int']['SEE_Rate']
Plus = Dict['Carrington-SEP-Plus2Sigma-Int']['SEE_Rate']
Error = Dict['Carrington-SEP-Expected-Int']['SEE_Error']

plt.errorbar(Shielding, ExpectedInt, yerr=Error, label="GEO EVT 1-in-150 year peak SEP"
             , color=Expected, **default_errorbar_style)

plt.fill_between(Shielding, ExpectedInt, Plus, color=PlusColor, alpha=0.5, label="GEO EVT 1-in-150 year peak SEP +2\u03C3") 
plt.fill_between(Shielding, ExpectedInt, Minus, color=MinusColor, alpha=0.5, label="GEO EVT 1-in-150 year peak SEP -2\u03C3")


# Other EVT plots. Do not plot for now.
# plt.errorbar(Dict['Carrington-SEP-Plus2Sigma-Int']['Shielding'], Dict['Carrington-SEP-Plus2Sigma-Int']['SEE_Rate'],
#             yerr=Dict['Carrington-SEP-Plus2Sigma-Int']['SEE_Error'], label="Carrington SEP EVT +2\u03C3"
#             , color=PlusColor, **default_errorbar_style)
# plt.errorbar(Dict['Carrington-SEP-Minus2Sigma-Int']['Shielding'], Dict['Carrington-SEP-Minus2Sigma-Int']['SEE_Rate'],
#             yerr=Dict['Carrington-SEP-Minus2Sigma-Int']['SEE_Error'], label="Carrington SEP EVT -2\u03C3"
#             , color=MinusColor, **default_errorbar_style)

# Plot the GEO Solar Proton 5min Peak Flux
plt.errorbar(Dict['GEO-SolarProton-5minPeakFlux']['Shielding'], Dict['GEO-SolarProton-5minPeakFlux']['SEE_Rate'],
              yerr=Dict['GEO-SolarProton-5minPeakFlux']['SEE_Error'], marker='o', label="GEO CREME96 peak 5 min SEP"
              , color=CREME96Color, linestyle=':', **default_errorbar_style)


########################### LEO ###################################
# plt.errorbar(Dict['LEO-trapped-proton']['Shielding'], Dict['LEO-trapped-proton']['SEE_Rate'], 
#              yerr=Dict['LEO-trapped-proton']['SEE_Error'], label="LEO trapped proton", capsize=5, capthick=2)
# plt.errorbar(Dict['LEO-solar-proton']['Shielding'], Dict['LEO-solar-proton']['SEE_Rate'], 
#              yerr=Dict['LEO-solar-proton']['SEE_Error'], label="LEO Solar Proton", capsize=5, capthick=2)
# plt.errorbar(Dict['LEO-cosmic-proton']['Shielding'], Dict['LEO-cosmic-proton']['SEE_Rate'], 
#              yerr=Dict['LEO-cosmic-proton']['SEE_Error'], label="LEO Cosmic Proton", capsize=5, capthick=2)

# Combine the LEO trapped proton, Solar Proton and Cosmic data
Dict['LEO'] = {}
Dict['LEO']['Shielding'] = Dict['LEO-solar-proton']['Shielding']

Dict['LEO']['SEE_Rate'] = (  Dict['LEO-trapped-proton']['SEE_Rate'] 
                           + Dict['LEO-solar-proton']['SEE_Rate'] 
                           + Dict['LEO-cosmic-proton']['SEE_Rate'] )

Dict['LEO']['SEE_Error'] = ( np.sqrt(np.square(Dict['LEO-trapped-proton']['SEE_Error'])
                                     + np.square(Dict['LEO-solar-proton']['SEE_Error']) 
                                     + np.square(Dict['LEO-cosmic-proton']['SEE_Error']) ) )

plt.errorbar(Dict['LEO']['Shielding'], Dict['LEO']['SEE_Rate'], yerr=Dict['LEO']['SEE_Error'],
                label="LEO total 11-year mean", color=LEOColor, linestyle='-.', **default_errorbar_style)

########################### GEO ###################################
# plt.errorbar(Dict['GEO-solar-proton']['Shielding'], Dict['GEO-solar-proton']['SEE_Rate'],
                # yerr=Dict['GEO-solar-proton']['SEE_Error'], label="GEO Solar Proton", capsize=5, capthick=2)
# plt.errorbar(Dict['GEO-cosmic-proton']['Shielding'], Dict['GEO-cosmic-proton']['SEE_Rate'],
#                 yerr=Dict['GEO-cosmic-proton']['SEE_Error'], label="GEO Cosmic Proton", capsize=5, capthick=2)

# # Combine the GEO Solar Proton, Cosmic and trapped proton data
Dict['GEO'] = {}
Dict['GEO']['Shielding'] = Dict['GEO-solar-proton']['Shielding']

Dict['GEO']['SEE_Rate'] = ( Dict['GEO-solar-proton']['SEE_Rate']  
                           + Dict['GEO-cosmic-proton']['SEE_Rate'] )

Dict['GEO']['SEE_Error'] = ( np.sqrt(np.square(Dict['GEO-solar-proton']['SEE_Error'])  
                                     + np.square(Dict['GEO-cosmic-proton']['SEE_Error'])) )

plt.errorbar(Dict['GEO']['Shielding'], Dict['GEO']['SEE_Rate'], yerr=Dict['GEO']['SEE_Error'],
                label="GEO total 11-year mean", color=GEOColor, **default_errorbar_style)


########################### VAP ###################################
# plt.errorbar(Dict['VAP-trapped-proton']['Shielding'], Dict['VAP-trapped-proton']['SEE_Rate'],
                # yerr=Dict['VAP-trapped-proton']['SEE_Error'], label="VAP trapped proton", capsize=5, capthick=2)
# plt.errorbar(Dict['VAP-solar-proton']['Shielding'], Dict['VAP-solar-proton']['SEE_Rate'],
#                 yerr=Dict['VAP-solar-proton']['SEE_Error'], label="VAP Solar Proton", capsize=5, capthick=2)
# plt.errorbar(Dict['VAP-cosmic-proton']['Shielding'], Dict['VAP-cosmic-proton']['SEE_Rate'],
#                 yerr=Dict['VAP-cosmic-proton']['SEE_Error'], label="VAP Cosmic Proton", capsize=5, capthick=2)

# Combine the VAP Solar Proton, Cosmic and trapped proton data
Dict['VAP'] = {}
Dict['VAP']['Shielding'] = Dict['VAP-solar-proton']['Shielding']

Dict['VAP']['SEE_Rate'] = (+ Dict['VAP-trapped-proton']['SEE_Rate'] 
                           + Dict['VAP-solar-proton']['SEE_Rate'] 
                           + Dict['VAP-cosmic-proton']['SEE_Rate'] )

Dict['VAP']['SEE_Error'] = ( np.sqrt(np.square(Dict['VAP-trapped-proton']['SEE_Error'])
                                     + np.square(Dict['VAP-solar-proton']['SEE_Error']) 
                                     + np.square(Dict['VAP-cosmic-proton']['SEE_Error']) ) )

plt.errorbar(Dict['VAP']['Shielding'], Dict['VAP']['SEE_Rate'], yerr=Dict['VAP']['SEE_Error']
                , label="VAP total 11-year mean", color=VAPColor, linestyle='--', **default_errorbar_style)


plt.axhline(1/(8e+3 * 60 * 60 * 24), linestyle='-', label="1 upset per kByte per day", color='black')
plt.axhline(1/(8e+6 * 60 * 60 * 24), linestyle='--', label="1 upset per MByte per day", color='black')
# plt.axhline(1/(8e+9 * 60 * 60 * 24), linestyle=':', label="1 upset per GByte per day", color='black')

# plt.ylim(1e-17, 2e-11)
plt.ylim(2e-13, 3e-6)
plt.yscale("log")
plt.grid()
plt.title(CrossectionName + "\nProton Fluence Cross Section SEU Rate Estimates")

plt.xlabel("Aluminium Shielding Thickness")
plt.ylabel("Single Event Upset Rate [s-1 bit-1]")

# Create a new order for the legend
handles, labels = plt.gca().get_legend_handles_labels()
order = [0, 4, 1, 5, 7, 6, 8, 2, 3]  # Adjust this list to reorder as needed
plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='lower right')

plt.savefig("/l/triton_work/Fluence_Histograms/CarringtonShielded/" + CrossectionName + " Rates.pdf", format='pdf', bbox_inches="tight")
# plt.show()