import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Read.ReadSD2Q import readSDQ2
from Dependencies.TotalDose import totalDose
from Dependencies.MergeTotalDose import mergeTotalDose



# Define default errorbar style
default_errorbar_style = {
    'capsize': 1.5,     # Size of the error bar caps
    'elinewidth': 1.5,    # Width of the error bar lines
    'capthick': 2,      # Thickness of the error bar caps
    'linestyle':''      # No line connecting the points
}

Expected = 'blue' # Blue
PlusColor = 'C1'  # Orange
MinusColor = 'C2' # Green
LEOColor = 'C8'   # Yellow
MEOColor = 'C9'   # Turquoise
VAPColor = 'C3'   # Red
GEOColor = 'C7'   # Grey


# Aluminium Shielding Thickness [mm]
# 101 Tiles from 0 to 10 mm
x = np.linspace(0, 10, num=101, dtype=float, endpoint=True)
# print(x)
####### kRad line #########
CriticalDose1 = [1 for i in x]
CriticalDose10 = [10 for i in x]

Path = "/scratch/work/fetzera1/GRAS/2D/Shielding_Curves/Carrington/"
# The carrington fluxes are per second, but we are interested in the event fluence which corresponds to 1.23 d
PerStoEvent = 1.23 * 24 * 60 * 60

for Species in ['Electrons', 'Protons']: 

    plt.figure(0, [4, 6])
    plt.plot(x, CriticalDose1, color='k', linewidth=2, label='1 krad')
    plt.plot(x, CriticalDose10, '--', color='k', linewidth=2, label='10 krad')

    ##################################### ELECTRON EVT Carrington Spectra ########################################
    if Species == 'Electrons':
        #### Carrington Electron fluxes ######
        # EVT_Electron_NoExt = totalDose(Path + "Carrington-Electron-Expected-EVT2026-Arb/Res")
        # plt.errorbar(x, EVT_Electron_NoExt['dose'] * PerStoEvent,
        #              yerr=EVT_Electron_NoExt['error'] * PerStoEvent, label='EVT electron\nevent dose', 
        #              color=Expected, **default_errorbar_style)
        
        EVT_Electron_Ext = totalDose(Path + "Carrington-Electron-Expected-EVT2026-Arb-130keVto1000MeV-Production/Res")
        plt.errorbar(x, EVT_Electron_Ext['dose'] * PerStoEvent, 
                    yerr=EVT_Electron_Ext['error'] * PerStoEvent, label='EVT electron\nevent dose', 
                    color=Expected, **default_errorbar_style)
        
        # EVT_Electron_ExtDense = totalDose(Path + "Carrington-Electron-Expected-EVT2026-ArbDense-130keVto1000MeV/Res")
        # plt.errorbar(x, EVT_Electron_ExtDense['dose'] * PerStoEvent, 
        #             yerr=EVT_Electron_ExtDense['error'] * PerStoEvent, label='EVT electron\nevent dose',
        #             color=Expected, **default_errorbar_style)
        
        # Old Electron spectra
        # CarringtonElectronINTEGRALPowTabelated = totalDose(Path + "CarringtonElectronINTEGRALPowTabelated/Res")
        # plt.errorbar(x, CarringtonElectronINTEGRALPowTabelated['dose'] * PerStoEvent, 
        #             yerr=CarringtonElectronINTEGRALPowTabelated['error'] * PerStoEvent, label='EVT electron\nevent dose', 
        #             color=Expected, **default_errorbar_style)
        
        # CarringtonElectronDiffPowTabelated = totalDose(Path + "CarringtonElectronDiffPowTabelated/Res")
        # plt.errorbar(x, CarringtonElectronDiffPowTabelated['dose'] * PerStoEvent,
        #             yerr=CarringtonElectronDiffPowTabelated['error'] * PerStoEvent, label='Carrington electron Diff Pow Tabelated',
        #             color=PlusColor, **default_errorbar_style)
        
        # CarringtonElectronDiffPow = totalDose(Path + "CarringtonElectronDiffPow/Res")
        # plt.errorbar(x, CarringtonElectronDiffPow['dose'] * PerStoEvent,
        #             yerr=CarringtonElectronDiffPow['error'] * PerStoEvent, label='Carrington electron Diff Pow',
        #             color=MinusColor, **default_errorbar_style)   



    ##################################### Proton EVT Carrington Spectra ########################################
    if Species == 'Protons':
        #### Carrington SEP fluxes ######
        Ex_EVT2026 = totalDose(Path + "Carrington-SEP-Expected-EVT2026/Res")
        plt.errorbar(x, Ex_EVT2026['dose'] * PerStoEvent, Ex_EVT2026['error'] * PerStoEvent, label='EVT SEP event dose',
                    color=Expected, **default_errorbar_style)

        # Ex_EVT2026_ext = totalDose(Path + "Carrington-SEP-Expected-EVT2026-10to1000MeV/Res")
        # plt.errorbar(x, Ex_EVT2026_ext['dose'] * PerStoEvent, Ex_EVT2026_ext['error'] * PerStoEvent, label='EVT SEP event dose',
        #             color=Expected, **default_errorbar_style)

        # Ex_EVT2026_Fit = totalDose(Path + "Carrington-SEP-Expected-EVT2026-Fit/Res")
        # plt.errorbar(x, Ex_EVT2026_Fit['dose'] * PerStoEvent, Ex_EVT2026_Fit['error'] * PerStoEvent, label='EVT SEP event dose',
        #             color=Expected, **default_errorbar_style)

        # Ex_EVT2026_Fit_ext = totalDose(Path + "Carrington-SEP-Expected-EVT2026-Fit-10to1000MeV/Res")
        # plt.errorbar(x, Ex_EVT2026_Fit_ext['dose'] * PerStoEvent, Ex_EVT2026_Fit_ext['error'] * PerStoEvent, label='EVT SEP event dose',
        #             color=Expected, **default_errorbar_style)

        # Ex = totalDose(Path + "Carrington-SEP-Expected-Int/Res")
        # plt.errorbar(x, Ex['dose'] * PerStoEvent, Ex['error'] * PerStoEvent, label='Carrington-SEP-Expected-Int', 
        #             color='C5', **default_errorbar_style)

        # Minus = totalDose(Path + "Carrington-SEP-Minus2Sigma-Int/Res")
        # plt.fill_between(x, Ex['dose'] * PerStoEvent, Minus['dose'] * PerStoEvent, color=MinusColor, alpha=0.5, 
        #                 label='EVT SEP − 2\u03C3')

        Minus_EVT2026 = totalDose(Path + "Carrington-SEP-Minus2Sigma-EVT2026/Res")
        plt.fill_between(x, Ex_EVT2026['dose'] * PerStoEvent, Minus_EVT2026['dose'] * PerStoEvent, color=MinusColor, alpha=0.5, 
                        label='EVT SEP − 2\u03C3')
        
        # Plus = totalDose(Path + "Carrington-SEP-Plus2Sigma-Int/Res")
        # plt.fill_between(x, Plus['dose'] * PerStoEvent, Ex['dose'] * PerStoEvent, color=PlusColor, alpha=0.5, 
        #                 label='EVT SEP + 2\u03C3')
        
        Plus_EVT2026 = totalDose(Path + "Carrington-SEP-Plus2Sigma-EVT2026/Res")
        plt.fill_between(x, Plus_EVT2026['dose'] * PerStoEvent, Ex_EVT2026['dose'] * PerStoEvent, color=PlusColor, alpha=0.5, 
                        label='EVT SEP + 2\u03C3')

        # Carrington_SEP_Expected_Diff = totalDose(Path + "Carrington-SEP-Expected-Diff/Res")
        # plt.errorbar(x, Carrington_SEP_Expected_Diff['dose'] * PerStoEvent, 
        #              Carrington_SEP_Expected_Diff['error'] * PerStoEvent,
        #             label='Carrington SEP Expected Diff', color=Expected, **default_errorbar_style)

        # Carrington_SEP_Expected_Int_With0 = totalDose(Path + "Carrington-SEP-Expected-Int-With0/Res")
        # plt.errorbar(x, Carrington_SEP_Expected_Int_With0['dose'] * PerStoEvent, 
        #             Carrington_SEP_Expected_Int_With0['error'] * PerStoEvent, 
        #             label='Carrington SEP Expected Int With 0', color='C8', **default_errorbar_style)
        
        # MinusWith0 = totalDose(Path + "Carrington-SEP-Minus2Sigma-Int-With0/Res")
        # plt.errorbar(x, MinusWith0['dose'] * PerStoEvent, MinusWith0['error'] * PerStoEvent,
        #             label='Carrington SEP − 2\u03C3 With 0', color=MinusColor, **default_errorbar_style)
        
        # PlusWith0 = totalDose(Path + "Carrington-SEP-Plus2Sigma-Int-With0/Res")
        # plt.errorbar(x, PlusWith0['dose'] * PerStoEvent, PlusWith0['error'] * PerStoEvent,
        #             label='Carrington SEP + 2\u03C3 With 0', color=PlusColor, **default_errorbar_style)

        ################################# Other extreme Proton Spectra ########################################
    
        ## CREME96 Peak Solar Proton Flux ##
        # solar_proton_5minPeakFlux = totalDose(Path + "Solar-Proton-5minPeakFlux/Res")
        # The flux is in 1/s cm2, therefore we need to multiply by to get the event dose in krad
        # plt.errorbar(x, solar_proton_5minPeakFlux['dose'] * 60 * 60 * 24 * 1.23, 
        #             yerr=solar_proton_5minPeakFlux['error'] * 60 * 60 * 24 * 1.23, color=GEOColor, 
        #             label="1 month peak flux on GEO", **default_errorbar_style)
        # solar_proton_5minPeakFlux_New = totalDose(Path + "Solar-Proton-5minPeakFlux_New/Res")
        # plt.errorbar(x, solar_proton_5minPeakFlux_New['dose'] * 60 * 60 * 24 * 1.23, 
        #             yerr=solar_proton_5minPeakFlux_New['error'] * 60 * 60 * 24 * 1.23, color='C9', 
        #             label="1 month peak flux on GEO", **default_errorbar_style)



    ################################## Nominal Mission Spectra ########################################
    
    #### LEO ######

    LEO_electron = totalDose(Path + "LEO-electron/Res")
    # plt.errorbar(x, LEO_electron['dose']*100/11, yerr=LEO_electron['error']*100/11, label="LEO Electrons 100 year dose", 
    #             color=LEOColor, **default_errorbar_style)

    LEO_trapped_proton = totalDose(Path + "LEO-trapped-proton/Res")
    # plt.errorbar(x, LEO_trapped_proton['dose'], yerr=LEO_trapped_proton['error'], label="LEO Trapped Protons", 
    #             color='C1', **default_errorbar_style)

    LEO_solar_proton = totalDose(Path + "LEO-solar-proton/Res") # Insigngificant on LEO
    # plt.errorbar(x, LEO_solar_proton['dose'], yerr=LEO_solar_proton['error'], label="LEO Solar Protons", 
    #             color='C2', **default_errorbar_style)

    LEO_cosmic_proton = totalDose(Path + "LEO-cosmic-proton/Res") # Negiligible
    # plt.errorbar(x, LEO_cosmic_proton['dose'], yerr=LEO_cosmic_proton['error'], label="LEO Cosmic Protons", 
    #             color='C3', **default_errorbar_style)

    LEO_cosmic_iron = totalDose(Path + "LEO-cosmic-iron/Res") # Negiligible
    # plt.errorbar(x, LEO_cosmic_iron['dose'], yerr=LEO_cosmic_iron['error'], label="LEO Cosmic Iron", 
    #             color='C4', **default_errorbar_style)

    LEO_Total = mergeTotalDose([LEO_electron, LEO_trapped_proton, LEO_solar_proton, LEO_cosmic_proton, LEO_cosmic_iron])
    plt.errorbar(x, LEO_Total['dose'] * 100/11, yerr=LEO_Total['error'] * 100/11, color=LEOColor, 
                label="LEO 100 year dose", **default_errorbar_style)

    # LEO_SD = readSDQ2("/scratch/work/fetzera1/Spectra/Carrington/LEO/spenvis_sqo.txt")
    # plt.plot(LEO_SD['Thickness'], ( LEO_SD['Electrons'] + LEO_SD['Bremsstrahlung'] ), 
    #         label="SHIELDOSE-2Q LEO Electrons + Bremsstrahlung", color='C0')
    # plt.plot(LEO_SD['Thickness'], LEO_SD['Trapped Protons'], label="SHIELDOSE-2Q LEO Trapped Protons", color='C1')
    # plt.plot(LEO_SD['Thickness'], LEO_SD['Total'] * 100/11, label="SHIELDOSE-2Q LEO Total Dose", color=LEOColor)


    #### MEO ######
    # MEO_electron = totalDose(Path + "MEO-electron/Res")
    # plt.errorbar(x, MEO_electron['dose'], yerr=MEO_electron['error'], label="MEO Electrons", 
    #             color='C0', **default_errorbar_style)
    
    # MEO_trapped_proton = totalDose(Path + "MEO-trapped-proton/Res")  # Effectively no trapped protons on MEO
    # plt.errorbar(x, MEO_trapped_proton['dose'], yerr=MEO_trapped_proton['error'], label="MEO Trapped Protons", 
    #             color='C1', **default_errorbar_style)

    # MEO_solar_proton = totalDose(Path + "MEO-solar-proton/Res")
    # plt.errorbar(x, MEO_solar_proton['dose'], yerr=MEO_solar_proton['error'], label="MEO Solar Protons", 
    #             color='C2', **default_errorbar_style)

    # MEO_cosmic_proton = totalDose(Path + "MEO-cosmic-proton/Res")
    # plt.errorbar(x, MEO_cosmic_proton['dose'], yerr=MEO_cosmic_proton['error'], label="MEO Cosmic Protons", 
    #             color='C3', **default_errorbar_style)

    # MEO_cosmic_iron = totalDose(Path + "MEO-cosmic-iron/Res")
    # plt.errorbar(x, MEO_cosmic_iron['dose'], yerr=MEO_cosmic_iron['error'], label="MEO Cosmic Iron", 
    #             color='C4', **default_errorbar_style)

    # MEO_Total = mergeTotalDose([MEO_electron, MEO_trapped_proton, MEO_solar_proton, MEO_cosmic_proton, MEO_cosmic_iron])
    # plt.errorbar(x, MEO_Total['dose'], yerr=MEO_Total['error'], color=MEOColor, 
    #             label="MEO Total", **default_errorbar_style)

    # MEO_SD = readSDQ2("/scratch/work/fetzera1/Spectra/Carrington/MEO/spenvis_sqo.txt")
    # plt.plot(MEO_SD['Thickness'], ( MEO_SD['Electrons'] + MEO_SD['Bremsstrahlung'] ), 
    #         label="SHIELDOSE-2Q MEO Electrons + Bremsstrahlung", color='C0')
    # plt.plot(MEO_SD['Thickness'], MEO_SD['Trapped Protons'], label="SHIELDOSE-2Q MEO Trapped Protons", color='C1')
    # plt.plot(MEO_SD['Thickness'], MEO_SD['Total'], label="SHIELDOSE-2Q MEO Total Dose", color=MEOColor)

    if Species == 'Electrons':
        #### VAP ######
        VAP_electron = totalDose(Path + "VAP-electron/Res")
        # plt.errorbar(x, VAP_electron['dose']/(11*12), yerr=VAP_electron['error']/(11*12), label="VAP Electrons 1 month dose", 
        #             color=VAPColor, **default_errorbar_style)
        
        VAP_trapped_proton = totalDose(Path + "VAP-trapped-proton/Res")
        # plt.errorbar(x, VAP_trapped_proton['dose'], yerr=VAP_trapped_proton['error'], label="VAP Trapped Protons", 
        #             color='C1', **default_errorbar_style)
        
        VAP_trapped_proton_integral = totalDose(Path + "VAP-trapped-proton-integral/Res")
        # plt.errorbar(x, VAP_trapped_proton_integral['dose'], yerr=VAP_trapped_proton_integral['error'], label="VAP Trapped Protons Integral", 
        #             color='black', **default_errorbar_style)
        
        VAP_solar_proton = totalDose(Path + "VAP-solar-proton/Res")
        # plt.errorbar(x, VAP_solar_proton['dose'], yerr=VAP_solar_proton['error'], label="VAP Solar Protons", 
        #             color='C2', **default_errorbar_style)
        
        VAP_cosmic_proton = totalDose(Path + "VAP-cosmic-proton/Res")
        # plt.errorbar(x, VAP_cosmic_proton['dose'], yerr=VAP_cosmic_proton['error'], label="VAP Cosmic Protons", 
        #             color='C3', **default_errorbar_style)
        
        VAP_cosmic_iron = totalDose(Path + "VAP-cosmic-iron/Res")
        # plt.errorbar(x, VAP_cosmic_iron['dose'], yerr=VAP_cosmic_iron['error'], label="VAP Cosmic Iron", 
        #             color='C4', **default_errorbar_style)


        VAP_Total = mergeTotalDose([VAP_electron, VAP_trapped_proton, VAP_solar_proton, VAP_cosmic_proton, VAP_cosmic_iron])
        plt.errorbar(x, VAP_Total['dose']/(11*12), yerr=VAP_Total['error']/(11*12), 
                    label="VAP 1 month dose\n(GRAS)", color=VAPColor, **default_errorbar_style)

        if Species == 'Electrons':
            VAP_SD = readSDQ2("/scratch/work/fetzera1/Spectra/Carrington/VAP/spenvis_sqo.txt")
            # plt.plot(VAP_SD['Thickness'], ( VAP_SD['Electrons'] + VAP_SD['Bremsstrahlung'] ), label="SHIELDOSE-2Q VAP Electrons + Bremsstrahlung", color='C0')
            # plt.plot(VAP_SD['Thickness'], VAP_SD['Trapped Protons'], label="SHIELDOSE-2Q VAP Trapped Protons", color='C1')
            plt.plot(VAP_SD['Thickness'], VAP_SD['Total'] /(11*12), '.', label="VAP 1 month dose\n(SHIELDOSE-2Q)", color=VAPColor) 
            
    
    # #### GEO ######
    GEO_electron = totalDose(Path + "GEO-electron/Res") # Similar but slightly less than VAP-electron, therefore not interesting for shielding
    # plt.errorbar(x, GEO_electron['dose']/11, yerr=GEO_electron['error']/11, label="GEO Electrons 1 year dose", 
    #              color=GEOColor, **default_errorbar_style)    
    GEO_trapped_proton = totalDose(Path + "GEO-trapped-proton/Res") # The trapped proton spectrum on GEO ends at 8MeV, therefore not intetrsting for shielding
    # plt.errorbar(x, GEO_trapped_proton['dose'], yerr=GEO_trapped_proton['error'], label="GEO Trapped Protons", 
    #              color='C1', **default_errorbar_style)

    GEO_solar_proton = totalDose(Path + "GEO-solar-proton/Res")
    # plt.errorbar(x, GEO_solar_proton['dose'], yerr=GEO_solar_proton['error'], label="GEO Solar Protons", 
    #              color='C2', **default_errorbar_style)

    GEO_cosmic_proton = totalDose(Path + "GEO-cosmic-proton/Res") # Cosmics have super low fluxes --> negligible TID
    # plt.errorbar(x, GEO_cosmic_proton['dose'], yerr=GEO_cosmic_proton['error'], label="GEO Cosmic Protons", 
    #              color='C3', **default_errorbar_style)

    GEO_cosmic_iron = totalDose(Path + "GEO-cosmic-iron/Res")     # Cosmics have super low fluxes --> negligible TID
    # plt.errorbar(x, GEO_cosmic_iron['dose'], yerr=GEO_cosmic_iron['error'], label="GEO Cosmic Iron", 
    #              color='C4', **default_errorbar_style)

    GEO_Total = mergeTotalDose([GEO_electron, GEO_trapped_proton, GEO_solar_proton, GEO_cosmic_proton, GEO_cosmic_iron])
    plt.errorbar(x, GEO_Total['dose']/11, yerr=GEO_Total['error']/11, color=GEOColor, label="GEO 1 year dose"
                , **default_errorbar_style)

    # GEO_SD = readSDQ2("/scratch/work/fetzera1/Spectra/Carrington/GEO/spenvis_sqo.txt")
    # plt.plot(GEO_SD['Thickness'], ( GEO_SD['Electrons'] + GEO_SD['Bremsstrahlung'] ), label="SHIELDOSE-2Q GEO Electrons + Bremsstrahlung", color='C0')
    # plt.plot(GEO_SD['Thickness'], GEO_SD['Trapped Protons'], label="SHIELDOSE-2Q GEO Trapped Protons", color='C1')
    # plt.plot(GEO_SD['Thickness'], GEO_SD['Total'] * 10/11, label="SHIELDOSE-2Q GEO Total Dose", color=GEOColor) 
    
    

    ######################## Plot formatting ############################

    #plt.title("Total Ionising Dose Behind Shielding")
    plt.xlabel("Aluminium Shielding Thickness [mm]")
    plt.ylabel("Total Ionising Dose [krad]")
    plt.yscale("log")
    # plt.xscale("log")
    # Increase the number of ticks on the y-axis
    # plt.yticks(np.logspace(-4, 4, num=9))

    plt.xlim(0.5, 10)
    plt.ylim(1e-1, 2e+2)
    plt.grid(which='both')

    ################### Re-ordering the legend labels ############################
    '''
    # Store the handles and labels when plotting
    handles, labels = plt.gca().get_legend_handles_labels()

    # Reorder the handles and labels to your desired order
    # Example: Suppose you want to reorder them based on some criteria or the original plotting order
    if Species == 'Electrons':
        order = [0, 1, 3, 4, 6, 5, 2]  # ELECTRON specify the desired order of indices
    if Species == 'Protons':
        order = [0, 1, 6, 5, 2, 4, 3]  # PROTON Plot specify the desired order of indices 

    # Apply the reordering
    handles = [handles[i] for i in order]
    labels = [labels[i] for i in order]

    # Then pass them to plt.legend()
    plt.legend(handles, labels)
    '''
    # External legend
    plt.legend()

    if Species == 'Electrons':
        plt.savefig(Path + "CompareCarringtonElectron.pdf", format='pdf', bbox_inches="tight")
    if Species == 'Protons':
        plt.savefig(Path + "CompareCarringtonProton.pdf", format='pdf', bbox_inches="tight")
    
    plt.close()
    # plt.show()
