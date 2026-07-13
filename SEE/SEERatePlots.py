import csv
import matplotlib.pyplot as plt
import numpy as np

Directory = "/home/anton/triton_work/GRAS/LET_Histograms/Carrington/"

CrossectionNames = ["Cypress CY62167GE30-45ZXI", "NanoXplore SEU"]  

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

for CrossectionName in CrossectionNames:


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

    plt.figure(figsize=(5, 7))
    ########################### Carrington SEP ###################################
    # Shielding = Dict['Carrington-SEP-Expected-Int']['Shielding']
    # ExpectedInt = Dict['Carrington-SEP-Expected-Int']['SEE_Rate']
    # Minus = Dict['Carrington-SEP-Minus2Sigma-Int']['SEE_Rate']
    # Plus = Dict['Carrington-SEP-Plus2Sigma-Int']['SEE_Rate']
    # Error = Dict['Carrington-SEP-Expected-Int']['SEE_Error']

    # plt.errorbar(Shielding, ExpectedInt, yerr=Error, label="GEO EVT 1-in-150 year peak SEP"
    #              , color=Expected, **default_errorbar_style)

    # plt.fill_between(Shielding, ExpectedInt, Plus, color=PlusColor, alpha=0.5, label="GEO EVT 1-in-150 year peak SEP +2\u03C3") 
    # plt.fill_between(Shielding, ExpectedInt, Minus, color=MinusColor, alpha=0.5, label="GEO EVT 1-in-150 year peak SEP -2\u03C3")

    # Other EVT plots. Do not plot for now.
    # plt.errorbar(Dict['Carrington-SEP-Plus2Sigma-Int']['Shielding'], Dict['Carrington-SEP-Plus2Sigma-Int']['SEE_Rate'],
    #             yerr=Dict['Carrington-SEP-Plus2Sigma-Int']['SEE_Error'], label="Carrington SEP EVT +2\u03C3"
    #             , color=PlusColor, **default_errorbar_style)
    # plt.errorbar(Dict['Carrington-SEP-Minus2Sigma-Int']['Shielding'], Dict['Carrington-SEP-Minus2Sigma-Int']['SEE_Rate'],
    #             yerr=Dict['Carrington-SEP-Minus2Sigma-Int']['SEE_Error'], label="Carrington SEP EVT -2\u03C3"
    #             , color=MinusColor, **default_errorbar_style)

    # Carrington Electrons
    # plt.errorbar(Dict['CarringtonElectronINTEGRALPowTabelated']['Shielding'], Dict['CarringtonElectronINTEGRALPowTabelated']['SEE_Rate'],
    #                 yerr=Dict['CarringtonElectronINTEGRALPowTabelated']['SEE_Error'], capsize=5, elinewidth=1, capthick=2, label="Carrington Electron SEU Rate", color='C3')



    # # Plot the GEO Solar Proton 5min Peak flux
    # plt.errorbar(Dict['GEO-SolarProton-5minPeakFlux']['Shielding'], Dict['GEO-SolarProton-5minPeakFlux']['SEE_Rate'],
    #               yerr=Dict['GEO-SolarProton-5minPeakFlux']['SEE_Error'], marker='o', label="GEO CREME96 peak 5 min SEP Proton"
    #               , color='C1', linestyle=':', **default_errorbar_style)

    # plt.errorbar(Dict['GEO-CREME96-Solar-Helium']['Shielding'], Dict['GEO-CREME96-Solar-Helium']['SEE_Rate'],
    #               yerr=Dict['GEO-CREME96-Solar-Helium']['SEE_Error'], marker='o', label="GEO CREME96 peak 5 min SEP Helium"
    #               , color='C2', linestyle='--', **default_errorbar_style)

    # plt.errorbar(Dict['GEO-CREME96-Solar-Oxygen']['Shielding'], Dict['GEO-CREME96-Solar-Oxygen']['SEE_Rate'],
    #               yerr=Dict['GEO-CREME96-Solar-Oxygen']['SEE_Error'], marker='o', label="GEO CREME96 peak 5 min SEP Oxygen"
    #                 , color='C3', linestyle='-.', **default_errorbar_style)

    # plt.errorbar(Dict['GEO-CREME96-Solar-Iron']['Shielding'], Dict['GEO-CREME96-Solar-Iron']['SEE_Rate'],
    #               yerr=Dict['GEO-CREME96-Solar-Iron']['SEE_Error'], marker='o', label="GEO CREME96 peak 5 min SEP Iron"
    #                 , color='C4', linestyle='-', **default_errorbar_style)

    # plt.errorbar(Dict['GEO-CREME96-Solar-Carbon']['Shielding'], Dict['GEO-CREME96-Solar-Carbon']['SEE_Rate'],
    #               yerr=Dict['GEO-CREME96-Solar-Carbon']['SEE_Error'], marker='o', label="GEO CREME96 peak 5 min SEP Carbon"
    #                 , color='C5', linestyle=':', **default_errorbar_style)

    # plt.errorbar(Dict['GEO-CREME96-Solar-Silicon']['Shielding'], Dict['GEO-CREME96-Solar-Silicon']['SEE_Rate'],
    #               yerr=Dict['GEO-CREME96-Solar-Silicon']['SEE_Error'], marker='o', label="GEO CREME96 peak 5 min SEP Silicon"
    #                 , color='C6', linestyle='--', **default_errorbar_style)

    # plt.errorbar(Dict['GEO-CREME96-Solar-Calcium']['Shielding'], Dict['GEO-CREME96-Solar-Calcium']['SEE_Rate'],
    #               yerr=Dict['GEO-CREME96-Solar-Calcium']['SEE_Error'], marker='o', label="GEO CREME96 peak 5 min SEP Calcium"
    #                 , color='C7', linestyle='-.', **default_errorbar_style)

    # Dict['GEO-CREME96-Solar-Total'] = {}
    # Dict['GEO-CREME96-Solar-Total']['Shielding'] = Dict['GEO-SolarProton-5minPeakFlux']['Shielding']
    # Dict['GEO-CREME96-Solar-Total']['SEE_Rate'] = (Dict['GEO-SolarProton-5minPeakFlux']['SEE_Rate']
    #                                                + Dict['GEO-CREME96-Solar-Helium']['SEE_Rate'] 
    #                                                + Dict['GEO-CREME96-Solar-Oxygen']['SEE_Rate'] 
    #                                                + Dict['GEO-CREME96-Solar-Iron']['SEE_Rate']
    #                                                + Dict['GEO-CREME96-Solar-Carbon']['SEE_Rate']
    #                                                + Dict['GEO-CREME96-Solar-Silicon']['SEE_Rate']
    #                                                + Dict['GEO-CREME96-Solar-Calcium']['SEE_Rate'])
    # Dict['GEO-CREME96-Solar-Total']['SEE_Error'] = (np.sqrt(np.square(Dict['GEO-SolarProton-5minPeakFlux']['SEE_Error']) 
    #                                                         + np.square(Dict['GEO-CREME96-Solar-Helium']['SEE_Error']) 
    #                                                         + np.square(Dict['GEO-CREME96-Solar-Oxygen']['SEE_Error']) 
    #                                                         + np.square(Dict['GEO-CREME96-Solar-Iron']['SEE_Error'])
    #                                                         + np.square(Dict['GEO-CREME96-Solar-Carbon']['SEE_Error'])
    #                                                         + np.square(Dict['GEO-CREME96-Solar-Silicon']['SEE_Error'])
    #                                                         + np.square(Dict['GEO-CREME96-Solar-Calcium']['SEE_Error'])))

    # plt.errorbar(Dict['GEO-CREME96-Solar-Total']['Shielding'], Dict['GEO-CREME96-Solar-Total']['SEE_Rate'],
    #               yerr=Dict['GEO-CREME96-Solar-Total']['SEE_Error'], marker='o', label="GEO CREME96 peak 5 min SEP Total"
    #                 , color='C0', linestyle='-', **default_errorbar_style)


    ########################### LEO ###################################
    # plt.errorbar(Dict['LEO-electron']['Shielding'], Dict['LEO-electron']['SEE_Rate'],
    #               yerr=Dict['LEO-electron']['SEE_Error'], label="LEO Electron", capsize=5, capthick=2)
    plt.errorbar(Dict['LEO-trapped-proton']['Shielding'], Dict['LEO-trapped-proton']['SEE_Rate'], 
                 yerr=Dict['LEO-trapped-proton']['SEE_Error'], label="LEO trapped proton", capsize=5, capthick=2)
    plt.errorbar(Dict['LEO-solar-proton']['Shielding'], Dict['LEO-solar-proton']['SEE_Rate'], 
                 yerr=Dict['LEO-solar-proton']['SEE_Error'], label="LEO Solar Proton", capsize=5, capthick=2)
    plt.errorbar(Dict['LEO-cosmic-proton']['Shielding'], Dict['LEO-cosmic-proton']['SEE_Rate'], 
                 yerr=Dict['LEO-cosmic-proton']['SEE_Error'], label="LEO Cosmic Proton", capsize=5, capthick=2)
    plt.errorbar(Dict['LEO-cosmic-iron']['Shielding'], Dict['LEO-cosmic-iron']['SEE_Rate'], 
                 yerr=Dict['LEO-cosmic-iron']['SEE_Error'], label="LEO Cosmic Iron", capsize=5, capthick=2)

    # Combine the LEO trapped proton, Solar Proton and Cosmic data
    Dict['LEO'] = {}
    Dict['LEO']['Shielding'] = Dict['LEO-solar-proton']['Shielding']
    Dict['LEO']['SEE_Rate'] = (# Dict['LEO-electron']['SEE_Rate'] 
                               Dict['LEO-trapped-proton']['SEE_Rate'] 
                               + Dict['LEO-solar-proton']['SEE_Rate'] 
                               + Dict['LEO-cosmic-proton']['SEE_Rate'] 
                               + Dict['LEO-cosmic-iron']['SEE_Rate'] )
    Dict['LEO']['SEE_Error'] = ( np.sqrt(#np.square(Dict['LEO-electron']['SEE_Error'])
                                         np.square(Dict['LEO-trapped-proton']['SEE_Error'])
                                         + np.square(Dict['LEO-solar-proton']['SEE_Error']) 
                                         + np.square(Dict['LEO-cosmic-proton']['SEE_Error']) 
                                         + np.square(Dict['LEO-cosmic-iron']['SEE_Error'])) )

    plt.errorbar(Dict['LEO']['Shielding'], Dict['LEO']['SEE_Rate'], yerr=Dict['LEO']['SEE_Error'],
                    label="LEO total 11-year mean", color=LEOColor, linestyle='-.', **default_errorbar_style)

    ########################### GEO ###################################
    # plt.errorbar(Dict['GEO-electron']['Shielding'], Dict['GEO-electron']['SEE_Rate'],
    #                 yerr=Dict['GEO-electron']['SEE_Error'], label="GEO Electron", capsize=5, capthick=2)
    # plt.errorbar(Dict['GEO-trapped-proton']['Shielding'], Dict['GEO-trapped-proton']['SEE_Rate'],
    #                 yerr=Dict['GEO-trapped-proton']['SEE_Error'], label="GEO trapped proton", capsize=5, capthick=2)
    # plt.errorbar(Dict['GEO-solar-proton']['Shielding'], Dict['GEO-solar-proton']['SEE_Rate'],
    #                 yerr=Dict['GEO-solar-proton']['SEE_Error'], label="GEO Solar Proton", capsize=5, capthick=2)
    # plt.errorbar(Dict['GEO-solar-helium']['Shielding'], Dict['GEO-solar-helium']['SEE_Rate'],
    #                 yerr=Dict['GEO-solar-helium']['SEE_Error'], label="GEO Solar Helium", capsize=5, capthick=2)
    # plt.errorbar(Dict['GEO-solar-oxygen']['Shielding'], Dict['GEO-solar-oxygen']['SEE_Rate'],
    #                 yerr=Dict['GEO-solar-oxygen']['SEE_Error'], label="GEO Solar Oxygen", capsize=5, capthick=2)
    # plt.errorbar(Dict['GEO-solar-iron']['Shielding'], Dict['GEO-solar-iron']['SEE_Rate'],
    #                 yerr=Dict['GEO-solar-iron']['SEE_Error'], label="GEO Solar Iron", capsize=5, capthick=2)
    # plt.errorbar(Dict['GEO-solar-nickel']['Shielding'], Dict['GEO-solar-nickel']['SEE_Rate'],
    #                 yerr=Dict['GEO-solar-nickel']['SEE_Error'], label="GEO Solar Nickel", capsize=5, capthick=2)
    # plt.errorbar(Dict['GEO-cosmic-proton']['Shielding'], Dict['GEO-cosmic-proton']['SEE_Rate'],
    #                 yerr=Dict['GEO-cosmic-proton']['SEE_Error'], label="GEO Cosmic Proton", capsize=5, capthick=2)
    # plt.errorbar(Dict['GEO-cosmic-iron']['Shielding'], Dict['GEO-cosmic-iron']['SEE_Rate'],
    #                 yerr=Dict['GEO-cosmic-iron']['SEE_Error'], label="GEO Cosmic Iron", capsize=5, capthick=2)


    # plt.errorbar(Dict['GEO-SAPPHIRE-Solar-AllHeavyIons']['Shielding'], Dict['GEO-SAPPHIRE-Solar-AllHeavyIons']['SEE_Rate'],
    #                 yerr=Dict['GEO-SAPPHIRE-Solar-AllHeavyIons']['SEE_Error'], label="GEO SAPPHIRE Solar Heavy Ions", color=SAP_PeakColor, capsize=5, capthick=2)
    # plt.errorbar(Dict['GEO-GCR-AllHeavyIons']['Shielding'], Dict['GEO-GCR-AllHeavyIons']['SEE_Rate'],
    #                 yerr=Dict['GEO-GCR-AllHeavyIons']['SEE_Error'], label="GEO GCR Heavy Ions", color='C4', capsize=5, capthick=2)


    # # Combine the GEO Solar Proton, Cosmic and trapped proton data
    # Dict['GEO'] = {}
    # Dict['GEO']['Shielding'] = Dict['GEO-solar-proton']['Shielding']
    # Dict['GEO']['SEE_Rate'] = (#Dict['GEO-electron']['SEE_Rate'] 
    #                            # Dict['GEO-trapped-proton']['SEE_Rate']
    #                         Dict['GEO-solar-proton']['SEE_Rate']
    #                         #+ Dict['GEO-solar-helium']['SEE_Rate']
    #                         #+ Dict['GEO-solar-oxygen']['SEE_Rate']
    #                         #+ Dict['GEO-solar-iron']['SEE_Rate']
    #                         #+ Dict['GEO-solar-nickel']['SEE_Rate']
    #                         #+ Dict['GEO-cosmic-iron']['SEE_Rate'] 
    #                         + Dict['GEO-cosmic-proton']['SEE_Rate']
    #                         + Dict['GEO-SAPPHIRE-Solar-AllHeavyIons']['SEE_Rate']
    #                         + Dict['GEO-GCR-AllHeavyIons']['SEE_Rate'])
    # Dict['GEO']['SEE_Error'] = ( np.sqrt(#np.square(Dict['GEO-electron']['SEE_Error']) 
    #                                     #np.square(Dict['GEO-trapped-proton']['SEE_Error'])
    #                                     np.square(Dict['GEO-solar-proton']['SEE_Error']) 
    #                                     #+ np.square(Dict['GEO-solar-helium']['SEE_Error'])
    #                                     #+ np.square(Dict['GEO-solar-oxygen']['SEE_Error'])
    #                                     #+ np.square(Dict['GEO-solar-iron']['SEE_Error'])
    #                                     #+ np.square(Dict['GEO-solar-nickel']['SEE_Error'])
    #                                     #+ np.square(Dict['GEO-cosmic-iron']['SEE_Error']) 
    #                                     + np.square(Dict['GEO-cosmic-proton']['SEE_Error'])) 
    #                                     + np.square(Dict['GEO-SAPPHIRE-Solar-AllHeavyIons']['SEE_Error']
    #                                     + np.square(Dict['GEO-GCR-AllHeavyIons']['SEE_Error'])))

    # plt.errorbar(Dict['GEO']['Shielding'], Dict['GEO']['SEE_Rate'], yerr=Dict['GEO']['SEE_Error'],
    #                 label="GEO total 11-year mean", color=GEOColor, **default_errorbar_style)


    ########################### VAP ###################################
    # plt.errorbar(Dict['VAP-electron']['Shielding'], Dict['VAP-electron']['SEE_Rate'],
    #                 yerr=Dict['VAP-electron']['SEE_Error'], label="VAP Electron", capsize=5, capthick=2)
    # plt.errorbar(Dict['VAP-trapped-proton']['Shielding'], Dict['VAP-trapped-proton']['SEE_Rate'],
    #                 yerr=Dict['VAP-trapped-proton']['SEE_Error'], label="VAP trapped proton", capsize=5, capthick=2)
    # plt.errorbar(Dict['VAP-solar-proton']['Shielding'], Dict['VAP-solar-proton']['SEE_Rate'],
    #                 yerr=Dict['VAP-solar-proton']['SEE_Error'], label="VAP Solar Proton", capsize=5, capthick=2)
    # plt.errorbar(Dict['VAP-cosmic-proton']['Shielding'], Dict['VAP-cosmic-proton']['SEE_Rate'],
    #                 yerr=Dict['VAP-cosmic-proton']['SEE_Error'], label="VAP Cosmic Proton", capsize=5, capthick=2)
    # plt.errorbar(Dict['VAP-cosmic-iron']['Shielding'], Dict['VAP-cosmic-iron']['SEE_Rate'],
    #                 yerr=Dict['VAP-cosmic-iron']['SEE_Error'], label="VAP Cosmic Iron", capsize=5, capthick=2)

    # # Combine the VAP Solar Proton, Cosmic and trapped proton data
    # Dict['VAP'] = {}
    # Dict['VAP']['Shielding'] = Dict['VAP-solar-proton']['Shielding']
    # Dict['VAP']['SEE_Rate'] = ( Dict['VAP-electron']['SEE_Rate'] 
    #                            + Dict['VAP-trapped-proton']['SEE_Rate'] 
    #                            + Dict['VAP-solar-proton']['SEE_Rate'] 
    #                            + Dict['VAP-cosmic-proton']['SEE_Rate'] 
    #                            + Dict['VAP-cosmic-iron']['SEE_Rate'] )
    # Dict['VAP']['SEE_Error'] = ( np.sqrt(np.square(Dict['VAP-electron']['SEE_Error'])
    #                                      + np.square(Dict['VAP-trapped-proton']['SEE_Error'])
    #                                      + np.square(Dict['VAP-solar-proton']['SEE_Error']) 
    #                                      + np.square(Dict['VAP-cosmic-proton']['SEE_Error']) 
    #                                      + np.square(Dict['VAP-cosmic-iron']['SEE_Error'])) )

    # plt.errorbar(Dict['VAP']['Shielding'], Dict['VAP']['SEE_Rate'], yerr=Dict['VAP']['SEE_Error']
    #                 , label="VAP total 11-year mean", color=VAPColor, linestyle='--', **default_errorbar_style)

    # plt.axhline(1/(8e+3 * 60 * 60 * 24), linestyle='-', label="1 upset per kByte per day", color='black')
    # plt.axhline(1/(8e+6 * 60 * 60 * 24), linestyle='--', label="1 upset per MByte per day", color='black')
    # # plt.axhline(1/(8e+9 * 60 * 60 * 24), linestyle=':', label="1 upset per GByte per day", color='black')

    # if CrossectionName == "Cypress CY62167GE30-45ZXI":
    #     plt.ylim(2e-13, 3e-6)
    #     pass
    # else:
    #     plt.ylim(1e-17, 2e-9)

    plt.yscale("log")
    plt.grid()
    plt.title(CrossectionName + "\nLET Cross Section SEU Rate Estimates")

    plt.xlabel("Aluminium Shielding Thickness")
    plt.ylabel("Single Event Upset Rate [s-1 bit-1]")

    # Manually add legend entries in the desired order
    # handles, labels = plt.gca().get_legend_handles_labels()

    # Create a new order for the legend
    # handles, labels = plt.gca().get_legend_handles_labels()
    # order = [0, 4, 1, 5, 7, 6, 8, 2, 3]  # Adjust this list to reorder as needed
    # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='lower right')
    plt.legend(loc='lower right')

    plt.savefig(Directory + CrossectionName + " Rates.pdf", format='pdf', bbox_inches="tight")
    # plt.show()