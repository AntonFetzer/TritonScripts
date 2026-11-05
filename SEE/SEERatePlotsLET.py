import matplotlib.pyplot as plt
import numpy as np
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Dependencies.SEERates import readRateCSVs
from Dependencies.CarringtonColors import Expected, PlusColor, MinusColor, LEOColor, GEOColor, MEOColor, CREME96Color, CREME96IonColor

Directory = "/scratch/work/fetzera1/GRAS/1Tile/LET_Histograms/Carrington/"

CrossectionNames = ["Cypress CY62167GE30-45ZXI", "NanoXplore NG-Medium"]  # List of cross-section names to process

# Define default errorbar style
default_errorbar_style = {
    'capsize': 4,        # Size of the error bar caps
    'elinewidth': 2,     # Width of the error bar lines
    'capthick': 2,       # Thickness of the error bar caps
}

for CrossectionName in CrossectionNames:


    # Read the CSV file of every run. The AllHeavyIons totals are written by SEERatePlotsAllSpecies.py
    Dict = readRateCSVs([Directory + "*/* " + CrossectionName + " SEE-Rate.csv",
                         Directory + "../Carrington-*-All-Species/* " + CrossectionName + " SEE-Rate.csv"], CrossectionName)

    plt.figure(figsize=(5, 8))
    ########################### Carrington SEP ###################################
    # Shielding = Dict['Carrington-SEP-Expected-Int']['Shielding']

    # ExpectedInt = Dict['Carrington-SEP-Expected-Int']['SEE_Rate']
    # Error = Dict['Carrington-SEP-Expected-Int']['SEE_Error']
    
    # Minus = Dict['Carrington-SEP-Minus2Sigma-Int']['SEE_Rate']
    # Plus = Dict['Carrington-SEP-Plus2Sigma-Int']['SEE_Rate']
    
    Shielding = Dict['Carrington-SEP-Expected-EVT2026']['Shielding']

    ExpectedInt = Dict['Carrington-SEP-Expected-EVT2026']['SEE_Rate']
    Error = Dict['Carrington-SEP-Expected-EVT2026']['SEE_Error']
    
    Minus = Dict['Carrington-SEP-Minus2Sigma-EVT2026']['SEE_Rate']
    Plus = Dict['Carrington-SEP-Plus2Sigma-EVT2026']['SEE_Rate']

    plt.errorbar(Shielding, ExpectedInt, yerr=Error, label="GEO EVT 1-in-150 year peak SEP"
                 , color=Expected, **default_errorbar_style)

    plt.fill_between(Shielding, ExpectedInt, Plus, color=PlusColor, alpha=0.5, label="GEO EVT 1-in-150 year peak SEP +2\u03C3") 
    plt.fill_between(Shielding, ExpectedInt, Minus, color=MinusColor, alpha=0.5, label="GEO EVT 1-in-150 year peak SEP -2\u03C3")


    # EVT2026 = Dict['Carrington-SEP-Expected-EVT2026']['SEE_Rate']
    # plt.errorbar(Shielding, EVT2026, label="SEP EVT2026", color='C4', **default_errorbar_style)    
    
    # EVT2026_Extr = Dict['Carrington-SEP-Expected-EVT2026-10to1000MeV']['SEE_Rate']
    # plt.errorbar(Shielding, EVT2026_Extr, label="SEP EVT2026 10-1000 MeV", color='C5', **default_errorbar_style)    
    
    # EVT2026_Fit = Dict['Carrington-SEP-Expected-EVT2026-Fit']['SEE_Rate']
    # plt.errorbar(Shielding, EVT2026_Fit, label="SEP EVT2026 Fit", color='C6', **default_errorbar_style)    
    
    # EVT2026_Fit_Extr = Dict['Carrington-SEP-Expected-EVT2026-Fit-10to1000MeV']['SEE_Rate']
    # plt.errorbar(Shielding, EVT2026_Fit_Extr, label="SEP EVT2026 Fit 10-1000 MeV", color='C7', **default_errorbar_style)

    # Carrington_SEP_Expected_Int = Dict['Carrington-SEP-Expected-Int']['SEE_Rate']
    # plt.errorbar(Shielding, Carrington_SEP_Expected_Int, label="Old", color='C8', **default_errorbar_style)


    # Other EVT plots. Do not plot for now.
    # plt.errorbar(Dict['Carrington-SEP-Plus2Sigma-Int']['Shielding'], Dict['Carrington-SEP-Plus2Sigma-Int']['SEE_Rate'],
    #             yerr=Dict['Carrington-SEP-Plus2Sigma-Int']['SEE_Error'], label="Carrington SEP EVT +2\u03C3"
    #             , color=PlusColor, **default_errorbar_style)
    # plt.errorbar(Dict['Carrington-SEP-Minus2Sigma-Int']['Shielding'], Dict['Carrington-SEP-Minus2Sigma-Int']['SEE_Rate'],
    #             yerr=Dict['Carrington-SEP-Minus2Sigma-Int']['SEE_Error'], label="Carrington SEP EVT -2\u03C3"
    #             , color=MinusColor, **default_errorbar_style)

    # Carrington Electrons
    # plt.errorbar(Dict['CarringtonElectronINTEGRALPowTabelated']['Shielding'],
    #               Dict['CarringtonElectronINTEGRALPowTabelated']['SEE_Rate'],
    #                 yerr=Dict['CarringtonElectronINTEGRALPowTabelated']['SEE_Error'],
    #                   label="Carrington Electron SEU Rate", color='C3', **default_errorbar_style)



    # # Plot the GEO Solar Proton 5min Peak flux
    plt.errorbar(Dict['GEO-SolarProton-5minPeakFlux']['Shielding'], Dict['GEO-SolarProton-5minPeakFlux']['SEE_Rate'],
                  yerr=Dict['GEO-SolarProton-5minPeakFlux']['SEE_Error'], marker='o', label="GEO CREME96 peak 5 min Proton"
                  , color=CREME96Color, linestyle=':', **default_errorbar_style)

    # plt.errorbar(Dict['GEO-CREME96-Solar-Iron']['Shielding'], Dict['GEO-CREME96-Solar-Iron']['SEE_Rate'],
    #               yerr=Dict['GEO-CREME96-Solar-Iron']['SEE_Error'], marker='o', label="GEO CREME96 peak 5 min SEP Iron"
    #                 , color='C4', linestyle='-', **default_errorbar_style)

    plt.errorbar(Dict['GEO-CREME96-Solar-AllHeavyIons']['Shielding'], Dict['GEO-CREME96-Solar-AllHeavyIons']['SEE_Rate'],
                  yerr=Dict['GEO-CREME96-Solar-AllHeavyIons']['SEE_Error'], marker='o', label="GEO CREME96 peak 5 min\nHeavy Ions GRAS LET Folding"
                    , color=CREME96IonColor, linestyle='-.', **default_errorbar_style)


    ########################### LEO ###################################
    # plt.errorbar(Dict['LEO-electron']['Shielding'], Dict['LEO-electron']['SEE_Rate'],
    #               yerr=Dict['LEO-electron']['SEE_Error'], label="LEO Electron", **default_errorbar_style)
    # plt.errorbar(Dict['LEO-trapped-proton']['Shielding'], Dict['LEO-trapped-proton']['SEE_Rate'], 
    #              yerr=Dict['LEO-trapped-proton']['SEE_Error'], label="LEO trapped proton", **default_errorbar_style)
    # plt.errorbar(Dict['LEO-solar-proton']['Shielding'], Dict['LEO-solar-proton']['SEE_Rate'], 
    #              yerr=Dict['LEO-solar-proton']['SEE_Error'], label="LEO Solar Proton", **default_errorbar_style)
    # plt.errorbar(Dict['LEO-cosmic-proton']['Shielding'], Dict['LEO-cosmic-proton']['SEE_Rate'], 
    #              yerr=Dict['LEO-cosmic-proton']['SEE_Error'], label="LEO Cosmic Proton", **default_errorbar_style)
    # plt.errorbar(Dict['LEO-cosmic-iron']['Shielding'], Dict['LEO-cosmic-iron']['SEE_Rate'], 
    #             yerr=Dict['LEO-cosmic-iron']['SEE_Error'], label="LEO Cosmic Iron", **default_errorbar_style)
    
    # plt.errorbar(Dict['LEO-SAPPHIRE-Solar-AllHeavyIons']['Shielding'], Dict['LEO-SAPPHIRE-Solar-AllHeavyIons']['SEE_Rate'],
    #                 yerr=Dict['LEO-SAPPHIRE-Solar-AllHeavyIons']['SEE_Error'], label="LEO SAPPHIRE Solar Heavy Ions", **default_errorbar_style)
    # plt.errorbar(Dict['LEO-GCR-AllHeavyIons']['Shielding'], Dict['LEO-GCR-AllHeavyIons']['SEE_Rate'],
    #                 yerr=Dict['LEO-GCR-AllHeavyIons']['SEE_Error'], label="LEO GCR Heavy Ions", **default_errorbar_style)

    # Combine the LEO trapped proton, Solar Proton and Cosmic data
    Dict['LEO'] = {}
    Dict['LEO']['Shielding'] = Dict['LEO-solar-proton']['Shielding']

    Dict['LEO']['SEE_Rate'] = (Dict['LEO-trapped-proton']['SEE_Rate'] 
                               + Dict['LEO-solar-proton']['SEE_Rate'] 
                               + Dict['LEO-cosmic-proton']['SEE_Rate'] 
                               + Dict['LEO-SAPPHIRE-Solar-AllHeavyIons']['SEE_Rate']
                               + Dict['LEO-GCR-AllHeavyIons']['SEE_Rate'])
    
    Dict['LEO']['SEE_Error'] = ( np.sqrt(np.square(Dict['LEO-trapped-proton']['SEE_Error'])
                                         + np.square(Dict['LEO-solar-proton']['SEE_Error']) 
                                         + np.square(Dict['LEO-cosmic-proton']['SEE_Error']) 
                                         + np.square(Dict['LEO-SAPPHIRE-Solar-AllHeavyIons']['SEE_Error'])
                                         + np.square(Dict['LEO-GCR-AllHeavyIons']['SEE_Error'])))

    plt.errorbar(Dict['LEO']['Shielding'], Dict['LEO']['SEE_Rate'], yerr=Dict['LEO']['SEE_Error'],
                    label="LEO total 11-year mean", color=LEOColor, linestyle='-.', **default_errorbar_style)


    #################################### MEO ###################################
    # plt.errorbar(Dict['MEO-electron']['Shielding'], Dict['MEO-electron']['SEE_Rate'],
    #                 yerr=Dict['MEO-electron']['SEE_Error'], label="MEO Electron", **default_errorbar_style)
    plt.errorbar(Dict['MEO-trapped-proton']['Shielding'], Dict['MEO-trapped-proton']['SEE_Rate'],
                    yerr=Dict['MEO-trapped-proton']['SEE_Error'], label="MEO trapped proton", **default_errorbar_style)
    # plt.errorbar(Dict['MEO-solar-proton']['Shielding'], Dict['MEO-solar-proton']['SEE_Rate'],
    #                 yerr=Dict['MEO-solar-proton']['SEE_Error'], label="MEO Solar Proton", **default_errorbar_style)
    plt.errorbar(Dict['MEO-cosmic-proton']['Shielding'], Dict['MEO-cosmic-proton']['SEE_Rate'],
                    yerr=Dict['MEO-cosmic-proton']['SEE_Error'], label="MEO Cosmic Proton", **default_errorbar_style)
    plt.errorbar(Dict['MEO-cosmic-iron']['Shielding'], Dict['MEO-cosmic-iron']['SEE_Rate'],
                    yerr=Dict['MEO-cosmic-iron']['SEE_Error'], label="MEO Cosmic Iron", **default_errorbar_style)

    plt.errorbar(Dict['MEO-SAPPHIRE-Solar-AllHeavyIons']['Shielding'], Dict['MEO-SAPPHIRE-Solar-AllHeavyIons']['SEE_Rate'],
                    yerr=Dict['MEO-SAPPHIRE-Solar-AllHeavyIons']['SEE_Error'], label="MEO SAPPHIRE Solar Heavy Ions", **default_errorbar_style)
    plt.errorbar(Dict['MEO-GCR-AllHeavyIons']['Shielding'], Dict['MEO-GCR-AllHeavyIons']['SEE_Rate'],
                    yerr=Dict['MEO-GCR-AllHeavyIons']['SEE_Error'], label="MEO GCR Heavy Ions", **default_errorbar_style)

    # Combine the MEO trapped proton, Solar Proton and Cosmic data
    Dict['MEO'] = {}
    Dict['MEO']['Shielding'] = Dict['MEO-trapped-proton']['Shielding']

    Dict['MEO']['SEE_Rate'] = (Dict['MEO-trapped-proton']['SEE_Rate']
                               + Dict['MEO-solar-proton']['SEE_Rate']
                                 + Dict['MEO-cosmic-proton']['SEE_Rate']
                                 + Dict['MEO-SAPPHIRE-Solar-AllHeavyIons']['SEE_Rate']
                                 + Dict['MEO-GCR-AllHeavyIons']['SEE_Rate'])

    Dict['MEO']['SEE_Error'] = ( np.sqrt(np.square(Dict['MEO-trapped-proton']['SEE_Error'])
                                         + np.square(Dict['MEO-solar-proton']['SEE_Error'])
                                         + np.square(Dict['MEO-cosmic-proton']['SEE_Error'])
                                         + np.square(Dict['MEO-SAPPHIRE-Solar-AllHeavyIons']['SEE_Error'])
                                         + np.square(Dict['MEO-GCR-AllHeavyIons']['SEE_Error'])))

    plt.errorbar(Dict['MEO']['Shielding'], Dict['MEO']['SEE_Rate'], yerr=Dict['MEO']['SEE_Error'],
                    label="MEO total 11-year mean", color=MEOColor, linestyle='-.', **default_errorbar_style)



    ########################### GEO ###################################
    # plt.errorbar(Dict['GEO-electron']['Shielding'], Dict['GEO-electron']['SEE_Rate'],
    #                 yerr=Dict['GEO-electron']['SEE_Error'], label="GEO Electron", **default_errorbar_style)
    # plt.errorbar(Dict['GEO-trapped-proton']['Shielding'], Dict['GEO-trapped-proton']['SEE_Rate'],
    #                 yerr=Dict['GEO-trapped-proton']['SEE_Error'], label="GEO trapped proton", **default_errorbar_style)
    # plt.errorbar(Dict['GEO-solar-proton']['Shielding'], Dict['GEO-solar-proton']['SEE_Rate'],
    #                 yerr=Dict['GEO-solar-proton']['SEE_Error'], label="GEO Solar Proton", **default_errorbar_style)
    # plt.errorbar(Dict['GEO-solar-iron']['Shielding'], Dict['GEO-solar-iron']['SEE_Rate'],
    #                 yerr=Dict['GEO-solar-iron']['SEE_Error'], label="GEO Solar Iron", **default_errorbar_style)
    # plt.errorbar(Dict['GEO-cosmic-proton']['Shielding'], Dict['GEO-cosmic-proton']['SEE_Rate'],
    #                 yerr=Dict['GEO-cosmic-proton']['SEE_Error'], label="GEO Cosmic Proton", **default_errorbar_style)
    # plt.errorbar(Dict['GEO-cosmic-iron']['Shielding'], Dict['GEO-cosmic-iron']['SEE_Rate'],
    #                 yerr=Dict['GEO-cosmic-iron']['SEE_Error'], label="GEO Cosmic Iron", **default_errorbar_style)


    # plt.errorbar(Dict['GEO-SAPPHIRE-Solar-AllHeavyIons']['Shielding'], Dict['GEO-SAPPHIRE-Solar-AllHeavyIons']['SEE_Rate'],
    #                 yerr=Dict['GEO-SAPPHIRE-Solar-AllHeavyIons']['SEE_Error'], label="GEO SAPPHIRE Solar Heavy Ions", **default_errorbar_style)
    # plt.errorbar(Dict['GEO-GCR-AllHeavyIons']['Shielding'], Dict['GEO-GCR-AllHeavyIons']['SEE_Rate'],
    #                 yerr=Dict['GEO-GCR-AllHeavyIons']['SEE_Error'], label="GEO GCR Heavy Ions", **default_errorbar_style)


    # Combine the GEO Solar Proton, Cosmic and trapped proton data
    Dict['GEO'] = {}
    Dict['GEO']['Shielding'] = Dict['GEO-solar-proton']['Shielding']

    Dict['GEO']['SEE_Rate'] = (Dict['GEO-trapped-proton']['SEE_Rate']
                            + Dict['GEO-solar-proton']['SEE_Rate']
                            + Dict['GEO-cosmic-proton']['SEE_Rate']
                            + Dict['GEO-SAPPHIRE-Solar-AllHeavyIons']['SEE_Rate']
                            + Dict['GEO-GCR-AllHeavyIons']['SEE_Rate'])
    
    Dict['GEO']['SEE_Error'] = (np.sqrt(np.square(Dict['GEO-trapped-proton']['SEE_Error'])
                                        + np.square(Dict['GEO-solar-proton']['SEE_Error']) 
                                        + np.square(Dict['GEO-cosmic-proton']['SEE_Error'])
                                        + np.square(Dict['GEO-SAPPHIRE-Solar-AllHeavyIons']['SEE_Error'])
                                        + np.square(Dict['GEO-GCR-AllHeavyIons']['SEE_Error'])))

    plt.errorbar(Dict['GEO']['Shielding'], Dict['GEO']['SEE_Rate'], yerr=Dict['GEO']['SEE_Error'],
                    label="GEO total 11-year mean", color=GEOColor, **default_errorbar_style)


    ########################### VAP ###################################
    # plt.errorbar(Dict['VAP-electron']['Shielding'], Dict['VAP-electron']['SEE_Rate'],
    #                 yerr=Dict['VAP-electron']['SEE_Error'], label="VAP Electron", **default_errorbar_style)
    # plt.errorbar(Dict['VAP-trapped-proton']['Shielding'], Dict['VAP-trapped-proton']['SEE_Rate'],
    #                 yerr=Dict['VAP-trapped-proton']['SEE_Error'], label="VAP trapped proton", **default_errorbar_style)
    # plt.errorbar(Dict['VAP-solar-proton']['Shielding'], Dict['VAP-solar-proton']['SEE_Rate'],
    #                 yerr=Dict['VAP-solar-proton']['SEE_Error'], label="VAP Solar Proton", **default_errorbar_style)
    # plt.errorbar(Dict['VAP-cosmic-proton']['Shielding'], Dict['VAP-cosmic-proton']['SEE_Rate'],
    #                 yerr=Dict['VAP-cosmic-proton']['SEE_Error'], label="VAP Cosmic Proton", **default_errorbar_style)
    # plt.errorbar(Dict['VAP-cosmic-iron']['Shielding'], Dict['VAP-cosmic-iron']['SEE_Rate'],
    #                 yerr=Dict['VAP-cosmic-iron']['SEE_Error'], label="VAP Cosmic Iron", **default_errorbar_style)

    # plt.errorbar(Dict['VAP-SAPPHIRE-Solar-AllHeavyIons']['Shielding'], Dict['VAP-SAPPHIRE-Solar-AllHeavyIons']['SEE_Rate'],
    #                 yerr=Dict['VAP-SAPPHIRE-Solar-AllHeavyIons']['SEE_Error'], label="VAP SAPPHIRE Solar Heavy Ions", **default_errorbar_style)
    # plt.errorbar(Dict['VAP-GCR-AllHeavyIons']['Shielding'], Dict['VAP-GCR-AllHeavyIons']['SEE_Rate'],
    #                 yerr=Dict['VAP-GCR-AllHeavyIons']['SEE_Error'], label="VAP GCR Heavy Ions", **default_errorbar_style)


    # # Combine the VAP Solar Proton, Cosmic and trapped proton data
    # Dict['VAP'] = {}
    # Dict['VAP']['Shielding'] = Dict['VAP-solar-proton']['Shielding']
    # Dict['VAP']['SEE_Rate'] = (# Dict['VAP-electron']['SEE_Rate'] 
    #                             Dict['VAP-trapped-proton']['SEE_Rate'] 
    #                            + Dict['VAP-solar-proton']['SEE_Rate'] 
    #                            + Dict['VAP-cosmic-proton']['SEE_Rate'] 
    #                            #+ Dict['VAP-cosmic-iron']['SEE_Rate']
    #                            + Dict['VAP-SAPPHIRE-Solar-AllHeavyIons']['SEE_Rate']
    #                              + Dict['VAP-GCR-AllHeavyIons']['SEE_Rate'])
    # Dict['VAP']['SEE_Error'] = ( np.sqrt(#np.square(Dict['VAP-electron']['SEE_Error'])
    #                                       np.square(Dict['VAP-trapped-proton']['SEE_Error'])
    #                                      + np.square(Dict['VAP-solar-proton']['SEE_Error']) 
    #                                      + np.square(Dict['VAP-cosmic-proton']['SEE_Error']) 
    #                                      #+ np.square(Dict['VAP-cosmic-iron']['SEE_Error'])
    #                                      + np.square(Dict['VAP-SAPPHIRE-Solar-AllHeavyIons']['SEE_Error'])
    #                                      + np.square(Dict['VAP-GCR-AllHeavyIons']['SEE_Error'])))
    # plt.errorbar(Dict['VAP']['Shielding'], Dict['VAP']['SEE_Rate'], yerr=Dict['VAP']['SEE_Error']
    #                 , label="VAP total 11-year mean", color=VAPColor, linestyle='--', **default_errorbar_style)

    plt.axhline(1/(8e+6 * 60 * 60 * 24), linestyle='--', label="1 upset per MByte per day", color='black')

    if CrossectionName == "Cypress CY62167GE30-45ZXI":
        # plt.ylim(2e-13, 3e-6)
        plt.axhline(1/(8e+3 * 60 * 60 * 24), linestyle='-', label="1 upset per kByte per day", color='black')
        # No CREME96 HUP curve: no defensible single-RPP sensitive volume exists for
        # this part (Serban et al. 2024 require four nested RPPs)
    else:
        # plt.ylim(1e-18, 3e-11)
        plt.axhline(1/(8e+9 * 60 * 60 * 24), linestyle=':', label="1 upset per GByte per day", color='black')

        # CREME96 HUP (Vanderbilt web tool): GEO peak 5 min SEP heavy ions (Z = 2-92),
        # IRPP with the vendor-documented Z = 2 um sensitive depth (RADSAGA 2017)
        CREME96HUPShielding = np.array(['1mm', '2mm', '4mm', '8mm', '16mm'])
        CREME96HUPRate = np.array([2.45e-11, 8.84e-13, 1.89e-13, 1.73e-14, 2.06e-15])
        plt.plot(CREME96HUPShielding, CREME96HUPRate, marker='s', label="GEO CREME96 peak 5 min\nHeavy Ions CREME 96 HUP IRPP"
                    , color='C9', linestyle=':')

    plt.yscale("log")
    plt.grid()
    plt.title(CrossectionName + "\nLET Cross Section SEU Rate Estimates")

    plt.xlabel("Aluminium Shielding Thickness")
    plt.ylabel("Single Event Upset Rate [s-1 bit-1]")

    # Manually add legend entries in the desired order
    # handles, labels = plt.gca().get_legend_handles_labels()

    # Create a new order for the legend
    # handles, labels = plt.gca().get_legend_handles_labels()
    # Order by label: matplotlib groups handles by artist type, so plain indices
    # shift whenever a line/patch/errorbar is added or removed
    # order = [labels.index(l) for l in [
    #     "GEO EVT 1-in-150 year peak SEP +2σ",
    #     "GEO EVT 1-in-150 year peak SEP",
    #     "GEO EVT 1-in-150 year peak SEP -2σ",
    #     "GEO CREME96 peak 5 min Proton",
    #     "GEO CREME96 peak 5 min\nHeavy Ions GRAS LET Folding",
    #     "GEO CREME96 peak 5 min\nHeavy Ions CREME 96 HUP IRPP",
    #     "VAP total 11-year mean",
    #     "GEO total 11-year mean",
    #     "LEO total 11-year mean",
    #     "1 upset per MByte per day",
    #     "1 upset per kByte per day" if CrossectionName == "Cypress CY62167GE30-45ZXI" else "1 upset per GByte per day",
    # ] if l in labels]
    # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='lower right')
    
    plt.legend()
    
    print("Saving " + Directory + CrossectionName + " Rates.pdf")
    plt.savefig(Directory + CrossectionName + " Rates.pdf", format='pdf', bbox_inches="tight")
    # plt.show()