import csv
import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unc

Expected = 'blue' # Blue
PlusColor = 'C1'  # Orange
MinusColor = 'C2' # Green
LEOColor = 'C8'   # Yellow
MEOColor = 'C9'   # Turquoise
VAPColor = 'C3'   # Red
GEOColor = 'C7'   # Grey

# Initialize the final data structure
Dict = {}

# Read the CSV file
with open('/l/triton_work/LET_Histograms/Carrington/SEERatesnanoXplore.csv', mode='r') as file:
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

        # Remove 0mm shielding
        if shielding == '0mm':
            continue

        # Store the data dictionary
        Dict[dataset_name]['Shielding'].append(shielding)
        Dict[dataset_name]['SEE_Rate'].append(see_rate)
        Dict[dataset_name]['SEE_Error'].append(see_error)
        Dict[dataset_name]['Entries_Contributing_To_SEE'].append(entries_contributing_to_see)

# Convert the lists to numpy arrays
for key in Dict:
    for inner_key in Dict[key]:
        Dict[key][inner_key] = np.array(Dict[key][inner_key])

# Print the final data structure for verification
# for key in Dict:
#     print(key)
#     print(Dict[key].keys())
#     print(Dict[key]['Shielding'])


# print(Dict['Carrington-SEP-Expected-Int'].keys())
# ['Shielding', 'SEE_Rate', 'SEE_Error', 'EntriesContributingToSEE']

# print(Dict['Carrington-SEP-Expected-Int']['Shielding'])
# ['1mm', '2mm', '4mm', '8mm', '16mm']

plt.figure(0, figsize=(5, 7))
########################### Carrington SEP ###################################

plt.errorbar(Dict['Carrington-SEP-Expected-Int']['Shielding'], Dict['Carrington-SEP-Expected-Int']['SEE_Rate'], 
             yerr=Dict['Carrington-SEP-Expected-Int']['SEE_Error'], capsize=5, elinewidth=1, capthick=2, label="Carrington SEP EVT SEU Rate Estimate", color='blue')

# plt.fill_between(Dict['Carrington-SEP-Expected-Int']['Shielding'], 
#                  Dict['Carrington-SEP-Expected-Int']['SEE_Rate'], 
#                  Dict['Carrington-SEP-Plus2Sigma-Int']['SEE_Rate'], color='C1', alpha=0.5, label="Carrington SEP EVT +2 Sigma")
# plt.fill_between(Dict['Carrington-SEP-Expected-Int']['Shielding'],
#                  Dict['Carrington-SEP-Expected-Int']['SEE_Rate'],
#                  Dict['Carrington-SEP-Minus2Sigma-Int']['SEE_Rate'], color='C2', alpha=0.5, label="Carrington SEP EVT -2 Sigma")


# Other EVT plots. Do not plot for now.
plt.errorbar(Dict['Carrington-SEP-Plus2Sigma-Int']['Shielding'], Dict['Carrington-SEP-Plus2Sigma-Int']['SEE_Rate'],
             yerr=Dict['Carrington-SEP-Plus2Sigma-Int']['SEE_Error'], capsize=5, elinewidth=1, capthick=2, label="Carrington SEP EVT +2 Sigma", color='C1')
plt.errorbar(Dict['Carrington-SEP-Minus2Sigma-Int']['Shielding'], Dict['Carrington-SEP-Minus2Sigma-Int']['SEE_Rate'],
             yerr=Dict['Carrington-SEP-Minus2Sigma-Int']['SEE_Error'], capsize=5, elinewidth=1, capthick=2, label="Carrington SEP EVT -2 Sigma", color='C2')

# Carrington Electrons
# plt.errorbar(Dict['CarringtonElectronINTEGRALPowTabelated']['Shielding'], Dict['CarringtonElectronINTEGRALPowTabelated']['SEE_Rate'],
                # yerr=Dict['CarringtonElectronINTEGRALPowTabelated']['SEE_Error'], capsize=5, elinewidth=1, capthick=2, label="Carrington Electron SEU Rate", color='C3')

# Plot the GEO Solar Proton 5min Peak Flux
plt.errorbar(Dict['GEO-SolarProton-5minPeakFlux']['Shielding'], Dict['GEO-SolarProton-5minPeakFlux']['SEE_Rate'],
              yerr=Dict['GEO-SolarProton-5minPeakFlux']['SEE_Error']
                , capsize=5, elinewidth=1, capthick=2, label="GEO Solar Proton 5min Peak SEU Rate", color='C7', linestyle=':')


########################### LEO ###################################
# plt.errorbar(Dict['LEO-electron']['Shielding'], Dict['LEO-electron']['SEE_Rate'],
#               yerr=Dict['LEO-electron']['SEE_Error'], label="LEO Electron", capsize=5, capthick=2)
# plt.errorbar(Dict['LEO-trapped-proton']['Shielding'], Dict['LEO-trapped-proton']['SEE_Rate'], 
#              yerr=Dict['LEO-trapped-proton']['SEE_Error'], label="LEO trapped proton", capsize=5, capthick=2)
# plt.errorbar(Dict['LEO-solar-proton']['Shielding'], Dict['LEO-solar-proton']['SEE_Rate'], 
#              yerr=Dict['LEO-solar-proton']['SEE_Error'], label="LEO Solar Proton", capsize=5, capthick=2)
# plt.errorbar(Dict['LEO-cosmic-proton']['Shielding'], Dict['LEO-cosmic-proton']['SEE_Rate'], 
#              yerr=Dict['LEO-cosmic-proton']['SEE_Error'], label="LEO Cosmic Proton", capsize=5, capthick=2)
# plt.errorbar(Dict['LEO-cosmic-iron']['Shielding'], Dict['LEO-cosmic-iron']['SEE_Rate'], 
#              yerr=Dict['LEO-cosmic-iron']['SEE_Error'], label="LEO Cosmic Iron", capsize=5, capthick=2)

# Combine the LEO trapped proton, Solar Proton and Cosmic data
Dict['LEO'] = {}
Dict['LEO']['Shielding'] = Dict['LEO-solar-proton']['Shielding']
Dict['LEO']['SEE_Rate'] = ( Dict['LEO-electron']['SEE_Rate'] 
                           + Dict['LEO-trapped-proton']['SEE_Rate'] 
                           + Dict['LEO-solar-proton']['SEE_Rate'] 
                           + Dict['LEO-cosmic-proton']['SEE_Rate'] 
                           + Dict['LEO-cosmic-iron']['SEE_Rate'] )
Dict['LEO']['SEE_Error'] = ( np.sqrt(np.square(Dict['LEO-trapped-proton']['SEE_Error'])
                                     + np.square(Dict['LEO-solar-proton']['SEE_Error']) 
                                     + np.square(Dict['LEO-cosmic-proton']['SEE_Error']) 
                                     + np.square(Dict['LEO-cosmic-iron']['SEE_Error'])) )

plt.errorbar(Dict['LEO']['Shielding'], Dict['LEO']['SEE_Rate'], yerr=Dict['LEO']['SEE_Error']
                , capsize=5, elinewidth=1, capthick=2, label="LEO Total Average SEU Rate", color=LEOColor, linestyle='-.')

########################### GEO ###################################
# plt.errorbar(Dict['GEO-electron']['Shielding'], Dict['GEO-electron']['SEE_Rate'],
#                 yerr=Dict['GEO-electron']['SEE_Error'], label="GEO Electron", capsize=5, capthick=2)
# plt.errorbar(Dict['GEO-trapped-proton']['Shielding'], Dict['GEO-trapped-proton']['SEE_Rate'],
#                 yerr=Dict['GEO-trapped-proton']['SEE_Error'], label="GEO trapped proton", capsize=5, capthick=2)
# plt.errorbar(Dict['GEO-solar-proton']['Shielding'], Dict['GEO-solar-proton']['SEE_Rate'],
                # yerr=Dict['GEO-solar-proton']['SEE_Error'], label="GEO Solar Proton", capsize=5, capthick=2)
# plt.errorbar(Dict['GEO-cosmic-proton']['Shielding'], Dict['GEO-cosmic-proton']['SEE_Rate'],
#                 yerr=Dict['GEO-cosmic-proton']['SEE_Error'], label="GEO Cosmic Proton", capsize=5, capthick=2)
# plt.errorbar(Dict['GEO-cosmic-iron']['Shielding'], Dict['GEO-cosmic-iron']['SEE_Rate'],
#                 yerr=Dict['GEO-cosmic-iron']['SEE_Error'], label="GEO Cosmic Iron", capsize=5, capthick=2)


# # Combine the GEO Solar Proton, Cosmic and trapped proton data
Dict['GEO'] = {}
Dict['GEO']['Shielding'] = Dict['GEO-solar-proton']['Shielding']
Dict['GEO']['SEE_Rate'] = ( Dict['GEO-solar-proton']['SEE_Rate'] 
                           + Dict['GEO-cosmic-iron']['SEE_Rate'] 
                           + Dict['GEO-cosmic-proton']['SEE_Rate'] )
Dict['GEO']['SEE_Error'] = ( np.sqrt(np.square(Dict['GEO-solar-proton']['SEE_Error']) 
                                     + np.square(Dict['GEO-cosmic-iron']['SEE_Error']) 
                                     + np.square(Dict['GEO-cosmic-proton']['SEE_Error'])) )

plt.errorbar(Dict['GEO']['Shielding'], Dict['GEO']['SEE_Rate'], yerr=Dict['GEO']['SEE_Error']
                , capsize=5, elinewidth=1, capthick=2, label="GEO Total Average SEU Rate", color=GEOColor)


########################### VAP ###################################
# plt.errorbar(Dict['VAP-electron']['Shielding'], Dict['VAP-electron']['SEE_Rate'],
#                 yerr=Dict['VAP-electron']['SEE_Error'], label="VAP Electron", capsize=5, capthick=2)
# plt.errorbar(Dict['VAP-trapped-proton']['Shielding'], Dict['VAP-trapped-proton']['SEE_Rate'],
                # yerr=Dict['VAP-trapped-proton']['SEE_Error'], label="VAP trapped proton", capsize=5, capthick=2)
# plt.errorbar(Dict['VAP-solar-proton']['Shielding'], Dict['VAP-solar-proton']['SEE_Rate'],
#                 yerr=Dict['VAP-solar-proton']['SEE_Error'], label="VAP Solar Proton", capsize=5, capthick=2)
# plt.errorbar(Dict['VAP-cosmic-proton']['Shielding'], Dict['VAP-cosmic-proton']['SEE_Rate'],
#                 yerr=Dict['VAP-cosmic-proton']['SEE_Error'], label="VAP Cosmic Proton", capsize=5, capthick=2)
# plt.errorbar(Dict['VAP-cosmic-iron']['Shielding'], Dict['VAP-cosmic-iron']['SEE_Rate'],
#                 yerr=Dict['VAP-cosmic-iron']['SEE_Error'], label="VAP Cosmic Iron", capsize=5, capthick=2)

# Combine the VAP Solar Proton, Cosmic and trapped proton data
Dict['VAP'] = {}
Dict['VAP']['Shielding'] = Dict['VAP-solar-proton']['Shielding']
Dict['VAP']['SEE_Rate'] = ( Dict['VAP-electron']['SEE_Rate'] 
                           + Dict['VAP-trapped-proton']['SEE_Rate'] 
                           + Dict['VAP-solar-proton']['SEE_Rate'] 
                           + Dict['VAP-cosmic-proton']['SEE_Rate'] 
                           + Dict['VAP-cosmic-iron']['SEE_Rate'] )
Dict['VAP']['SEE_Error'] = ( np.sqrt(np.square(Dict['VAP-trapped-proton']['SEE_Error'])
                                     + np.square(Dict['VAP-solar-proton']['SEE_Error']) 
                                     + np.square(Dict['VAP-cosmic-proton']['SEE_Error']) 
                                     + np.square(Dict['VAP-cosmic-iron']['SEE_Error'])) )

plt.errorbar(Dict['VAP']['Shielding'], Dict['VAP']['SEE_Rate'], yerr=Dict['VAP']['SEE_Error']
                , capsize=5, elinewidth=1, capthick=2, label="Van-Allen-Belt Probes Average SEU Rate", color=VAPColor, linestyle='-.')


# plt.axhline(1/8e+3, linestyle='--', label="1 Upset per kB per Second", color='black')
# plt.axhline(1/8e+6, linestyle=':', label="1 Upset per MByte per Second", color='black')
# plt.axhline(1/8e+9, linestyle='-.', label="1 Upset per GByte per Second", color='black')
plt.axhline(1/8e+12, linestyle=':', label="1 Upset per TByte per Second", color='black')
# plt.axhline(1/8e+15, linestyle=':', label="1 Upset per PByte per Second", color='black')
plt.axhline(1/(8e+6 * 60 * 60 * 24), linestyle='--', label="1 Upset per MByte per day", color='black')

#plt.ylim(1.5e-14, 4e-9)
# plt.ylim(2e-18, 3e-12)
plt.yscale("log")
plt.grid()
plt.title("NanoXplore NG-Medium SEU Rate Estimates")

plt.xlabel("Aluminium Shielding Thickness")
plt.ylabel("Single Event Upset Rate [s-1 bit-1]")

# Manually add legend entries in the desired order
# handles, labels = plt.gca().get_legend_handles_labels()

# Create a new order for the legend
# order = [
#     labels.index("Carrington SEP EVT +2 Sigma"), 
#     labels.index("Carrington SEP EVT SEU Rate Estimate"), 
#     labels.index("Carrington SEP EVT -2 Sigma"),
    # labels.index("1 Upset per GByte per Second"),
    # labels.index("1 Upset per MByte per day"),
    # labels.index("1 Upset per TByte per Second"),
    # labels.index("1 Upset per PByte per Second"),
    # labels.index("Van-Allen-Belt Probes Average SEU Rate"),
    # labels.index("GEO Solar Proton 5min Peak SEU Rate"),
    # labels.index("GEO Total Average SEU Rate"),
    # labels.index("LEO LEO Total Average SEU Rate"),
    # labels.index("Carrington Electron SEU Rate")
    # ]

# plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='lower left', bbox_to_anchor=(0, 0.25))
plt.legend()

plt.savefig("/l/triton_work/LET_Histograms/Carrington/" + CrossectionName + " Rates.pdf", format='pdf', bbox_inches="tight")
# plt.show()