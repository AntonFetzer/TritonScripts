import matplotlib.pyplot as plt
import numpy as np
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Dependencies.SEERates import readRateCSVs
from Dependencies.CarringtonColors import Expected, PlusColor, MinusColor, LEOColor, GEOColor, CREME96Color, SAP_PeakColor

Directory = "/scratch/work/fetzera1/GRAS/1Tile/Fluence_Histograms/CarringtonShielded/"

CrossectionNames = ["Cypress CY62167GE30-45ZXI", "NanoXplore NG-Medium"]  

# Define default errorbar style
default_errorbar_style = {
    'capsize': 4,        # Size of the error bar caps
    'elinewidth': 2,     # Width of the error bar lines
    'capthick': 2,       # Thickness of the error bar caps
}

for CrossectionName in CrossectionNames:

    # Read the CSV file of every run
    Dict = readRateCSVs([Directory + "*/* " + CrossectionName + " SEE-Rate.csv"], CrossectionName)


    plt.figure(figsize=(5, 8))
    ########################### Carrington SEP ###################################
    # The EVT2026 runs replace the -Int runs, which are in the deprecated CarringtonShielded_2xTooHigh
    Shielding = Dict['Carrington-SEP-Expected-EVT2026']['Shielding']
    EVT2026 = Dict['Carrington-SEP-Expected-EVT2026']['SEE_Rate']
    Minus = Dict['Carrington-SEP-Minus2Sigma-EVT2026']['SEE_Rate']
    Plus = Dict['Carrington-SEP-Plus2Sigma-EVT2026']['SEE_Rate']
    Error = Dict['Carrington-SEP-Expected-EVT2026']['SEE_Error']

    plt.errorbar(Shielding, EVT2026, yerr=Error, label="GEO EVT 1-in-150 year peak SEP"
                , color=Expected, **default_errorbar_style)

    plt.fill_between(Shielding, EVT2026, Plus, color=PlusColor, alpha=0.5, label="GEO EVT 1-in-150 year peak SEP +2\u03C3")
    plt.fill_between(Shielding, EVT2026, Minus, color=MinusColor, alpha=0.5, label="GEO EVT 1-in-150 year peak SEP -2\u03C3")



    EVT2026_Extr = Dict['Carrington-SEP-Expected-EVT2026-10to1000MeV']['SEE_Rate']
    EVT2026_Fit = Dict['Carrington-SEP-Expected-EVT2026-Fit']['SEE_Rate']
    EVT2026_Fit_Extr = Dict['Carrington-SEP-Expected-EVT2026-Fit-10to1000MeV']['SEE_Rate']


    plt.errorbar(Shielding, EVT2026_Extr, label="SEP EVT2026 10-1000 MeV"
                 , color='C5', **default_errorbar_style)
    plt.errorbar(Shielding, EVT2026_Fit, label="SEP EVT2026 Fit"
                 , color='C6', **default_errorbar_style)
    plt.errorbar(Shielding, EVT2026_Fit_Extr, label="SEP EVT2026 Fit 10-1000 MeV"
                 , color='C7', **default_errorbar_style)


    # Other EVT plots. Do not plot for now.
    # plt.errorbar(Dict['Carrington-SEP-Plus2Sigma-Int']['Shielding'], Dict['Carrington-SEP-Plus2Sigma-Int']['SEE_Rate'],
    #             yerr=Dict['Carrington-SEP-Plus2Sigma-Int']['SEE_Error'], label="Carrington SEP EVT +2\u03C3"
    #             , color=PlusColor, **default_errorbar_style)
    # plt.errorbar(Dict['Carrington-SEP-Minus2Sigma-Int']['Shielding'], Dict['Carrington-SEP-Minus2Sigma-Int']['SEE_Rate'],
    #             yerr=Dict['Carrington-SEP-Minus2Sigma-Int']['SEE_Error'], label="Carrington SEP EVT -2\u03C3"
    #             , color=MinusColor, **default_errorbar_style)

    # Plot the GEO Solar Proton 5min Peak Flux, once its new one-sided run exists
    if 'GEO-SolarProton-5minPeakFlux' in Dict:
        plt.errorbar(Dict['GEO-SolarProton-5minPeakFlux']['Shielding'], Dict['GEO-SolarProton-5minPeakFlux']['SEE_Rate'],
                    yerr=Dict['GEO-SolarProton-5minPeakFlux']['SEE_Error'], marker='o', label="GEO CREME96 peak 5 min Proton\nGRAS Shieleded Spectrum Folding"
                    , color=CREME96Color, linestyle=':', **default_errorbar_style)


    ########################### LEO ###################################
    # plt.errorbar(Dict['LEO-trapped-proton']['Shielding'], Dict['LEO-trapped-proton']['SEE_Rate'], 
    #              yerr=Dict['LEO-trapped-proton']['SEE_Error'], label="LEO trapped proton", capsize=5, capthick=2)
    # plt.errorbar(Dict['LEO-solar-proton']['Shielding'], Dict['LEO-solar-proton']['SEE_Rate'], 
    #              yerr=Dict['LEO-solar-proton']['SEE_Error'], label="LEO Solar Proton", capsize=5, capthick=2)
    # plt.errorbar(Dict['LEO-cosmic-proton']['Shielding'], Dict['LEO-cosmic-proton']['SEE_Rate'], 
    #              yerr=Dict['LEO-cosmic-proton']['SEE_Error'], label="LEO Cosmic Proton", capsize=5, capthick=2)

    # Combine the LEO trapped proton, Solar Proton and Cosmic data
    # Dict['LEO'] = {}
    # Dict['LEO']['Shielding'] = Dict['LEO-solar-proton']['Shielding']

    # Dict['LEO']['SEE_Rate'] = (  Dict['LEO-trapped-proton']['SEE_Rate'] 
    #                         + Dict['LEO-solar-proton']['SEE_Rate'] 
    #                         + Dict['LEO-cosmic-proton']['SEE_Rate'] )

    # Dict['LEO']['SEE_Error'] = ( np.sqrt(np.square(Dict['LEO-trapped-proton']['SEE_Error'])
    #                                     + np.square(Dict['LEO-solar-proton']['SEE_Error']) 
    #                                     + np.square(Dict['LEO-cosmic-proton']['SEE_Error']) ) )

    # plt.errorbar(Dict['LEO']['Shielding'], Dict['LEO']['SEE_Rate'], yerr=Dict['LEO']['SEE_Error'],
    #                 label="LEO total 11-year mean", color=LEOColor, linestyle='-.', **default_errorbar_style)

    ########################### GEO ###################################
    # plt.errorbar(Dict['GEO-solar-proton']['Shielding'], Dict['GEO-solar-proton']['SEE_Rate'],
                    # yerr=Dict['GEO-solar-proton']['SEE_Error'], label="GEO Solar Proton", capsize=5, capthick=2)
    # plt.errorbar(Dict['GEO-cosmic-proton']['Shielding'], Dict['GEO-cosmic-proton']['SEE_Rate'],
    #                 yerr=Dict['GEO-cosmic-proton']['SEE_Error'], label="GEO Cosmic Proton", capsize=5, capthick=2)

    # # Combine the GEO Solar Proton, Cosmic and trapped proton data
    # Dict['GEO'] = {}
    # Dict['GEO']['Shielding'] = Dict['GEO-solar-proton']['Shielding']

    # Dict['GEO']['SEE_Rate'] = ( Dict['GEO-solar-proton']['SEE_Rate']  
    #                         + Dict['GEO-cosmic-proton']['SEE_Rate'] )

    # Dict['GEO']['SEE_Error'] = ( np.sqrt(np.square(Dict['GEO-solar-proton']['SEE_Error'])  
    #                                     + np.square(Dict['GEO-cosmic-proton']['SEE_Error'])) )

    # plt.errorbar(Dict['GEO']['Shielding'], Dict['GEO']['SEE_Rate'], yerr=Dict['GEO']['SEE_Error'],
    #                 label="GEO total 11-year mean", color=GEOColor, **default_errorbar_style)


    ########################### VAP ###################################
    # plt.errorbar(Dict['VAP-trapped-proton']['Shielding'], Dict['VAP-trapped-proton']['SEE_Rate'],
                    # yerr=Dict['VAP-trapped-proton']['SEE_Error'], label="VAP trapped proton", capsize=5, capthick=2)
    # plt.errorbar(Dict['VAP-solar-proton']['Shielding'], Dict['VAP-solar-proton']['SEE_Rate'],
    #                 yerr=Dict['VAP-solar-proton']['SEE_Error'], label="VAP Solar Proton", capsize=5, capthick=2)
    # plt.errorbar(Dict['VAP-cosmic-proton']['Shielding'], Dict['VAP-cosmic-proton']['SEE_Rate'],
    #                 yerr=Dict['VAP-cosmic-proton']['SEE_Error'], label="VAP Cosmic Proton", capsize=5, capthick=2)

    # Combine the VAP Solar Proton, Cosmic and trapped proton data
    # Dict['VAP'] = {}
    # Dict['VAP']['Shielding'] = Dict['VAP-solar-proton']['Shielding']

    # Dict['VAP']['SEE_Rate'] = (+ Dict['VAP-trapped-proton']['SEE_Rate'] 
    #                         + Dict['VAP-solar-proton']['SEE_Rate'] 
    #                         + Dict['VAP-cosmic-proton']['SEE_Rate'] )

    # Dict['VAP']['SEE_Error'] = ( np.sqrt(np.square(Dict['VAP-trapped-proton']['SEE_Error'])
    #                                     + np.square(Dict['VAP-solar-proton']['SEE_Error']) 
    #                                     + np.square(Dict['VAP-cosmic-proton']['SEE_Error']) ) )

    # plt.errorbar(Dict['VAP']['Shielding'], Dict['VAP']['SEE_Rate'], yerr=Dict['VAP']['SEE_Error']
    #                 , label="VAP total 11-year mean", color=VAPColor, linestyle='--', **default_errorbar_style)

    plt.axhline(1/(8e+6 * 60 * 60 * 24), linestyle='--', label="1 upset per MByte per day", color='black')

    if CrossectionName == "Cypress CY62167GE30-45ZXI":
        plt.ylim(2e-13, 3e-6)
        plt.axhline(1/(8e+3 * 60 * 60 * 24), linestyle='-', label="1 upset per kByte per day", color='black')
    else:
        plt.ylim(1e-18, 3e-11)
        plt.axhline(1/(8e+9 * 60 * 60 * 24), linestyle=':', label="1 upset per GByte per day", color='black')

        # CREME96 PUP (Vanderbilt web tool): GEO peak 5 min SEP protons folded with the
        # vendor proton-energy Weibull; complete for the NG-Medium since its proton
        # direct-ionization contribution is zero (no PUP curve for the Cypress, whose
        # PDI-dominated response PUP structurally cannot represent)
        CREME96PUPShielding = np.array(['1mm', '2mm', '4mm', '8mm', '16mm'])
        CREME96PUPRate = np.array([1.98e-11, 1.56e-11, 1.28e-11, 8.10e-12, 4.53e-12])
        plt.plot(CREME96PUPShielding, CREME96PUPRate, marker='s', label="GEO CREME96 peak 5 min Proton\nCREME96 PUP"
                    , color=CREME96Color, linestyle=':')

    plt.yscale("log")
    plt.grid()
    plt.title(CrossectionName + "\nProton Energy Cross Section SEU Rate Estimates")

    plt.xlabel("Aluminium Shielding Thickness")
    plt.ylabel("Single Event Upset Rate [s-1 bit-1]")

    # Create a new order for the legend
    # handles, labels = plt.gca().get_legend_handles_labels()
    # Order by label: matplotlib groups handles by artist type, so plain indices
    # shift whenever a line/patch/errorbar is added or removed
    # order = [labels.index(l) for l in [
    #     "GEO EVT 1-in-150 year peak SEP +2σ",
    #     "GEO EVT 1-in-150 year peak SEP",
    #     "GEO EVT 1-in-150 year peak SEP -2σ",
    #     "GEO CREME96 peak 5 min Proton\nGRAS Shieleded Spectrum Folding",
    #     "GEO CREME96 peak 5 min Proton\nCREME96 PUP",
    #     "VAP total 11-year mean",
    #     "GEO total 11-year mean",
    #     "LEO total 11-year mean",
    #     "1 upset per MByte per day",
    #     "1 upset per kByte per day" if CrossectionName == "Cypress CY62167GE30-45ZXI" else "1 upset per GByte per day",
    # ] if l in labels]
    # plt.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='lower right')
    plt.legend(loc='lower right')

    plt.savefig(Directory + CrossectionName + " Rates.pdf", format='pdf', bbox_inches="tight")
    # plt.show()