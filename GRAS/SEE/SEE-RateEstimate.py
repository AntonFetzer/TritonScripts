import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from GRAS.Dependencies.TotalLETHistos import totalGRASLETHistos
from uncertainties import ufloat

'''
The functional form of the Weibull is:
F(x) = A (1- exp{-[(x-x0)/W] ** s}) where
    x = effective LET in MeV-cm2/milligram;
    F(x) = SEE cross-section in square-microns/bit;
    A0 = limiting or plateau cross-section;
    x0 = onset parameter, such that F(x) = 0 for x < x0;
    W = width parameter;
    s = a dimensionless exponent.
https://creme.isde.vanderbilt.edu/CREME-MC/help/weibull
'''


ThickList = ["0mm", "2mm", "4mm", "6mm"]

FitParams = "SEFI"

if FitParams == "SEU":
    CrossectionName = "Hercules SEU"
    A = 1.25E-11
    x0 = 0.95535
    W = 13.23393071
    s = 1
elif FitParams == "SET":
    CrossectionName = "Hercules SET"
    A = 4.02E-07
    x0 = 3.54024
    W = 18.31970977
    s = 1
elif FitParams == "SEFI":
    CrossectionName = "Hercules SEFI"
    A = 2.48E-06
    x0 = 0.95337
    W = 1.90630695
    s = 1


def f(LET):
    """
    Calculate the rate estimate based on the Linear Energy Transfer (LET) values.

    Parameters:
    LET (numpy.ndarray): Array of LET values.

    Returns:
    numpy.ndarray: Array of rate estimates.
    """
    
    result = np.zeros_like(LET)
    # Remove all values of LET that are less than x0
    mask = LET > x0
    # A * (1 - np.exp(-((x - x0) / W) ** s))
    result[mask] = A * (1 - np.exp(-((LET[mask] - x0) / W) ** s))

    return result


Paths = [
    "/l/triton_work/LET/Foresail1-Hercules/FS1-SolarProtons/",
    "/l/triton_work/LET/Foresail1-Hercules/FS1-CosmicProt/",
    "/l/triton_work/LET/Foresail1-Hercules/FS1-CosmicHe/",
    "/l/triton_work/LET/Foresail1-Hercules/FS1-CosmicFe/",
    "/l/triton_work/LET/Foresail1-Hercules/FS1-CosmicO/",
    "/l/triton_work/LET/Foresail1-Hercules/FS1-TrappedProtons/",
    "/l/triton_work/LET/Foresail1-Hercules/FS1-TrappedElectrons/"
    ]

DataName = [
    "FS1-Solar Protons",
    "FS1-Cosmic Protons",
    "FS1-Cosmic Helium",
    "FS1-Cosmic Iron",
    "FS1-Cosmic Oxygen",
    "FS1-Trapped Protons",
    "FS1-Trapped Electrons"
    ]

for Thick in ThickList:

    for P, Path in enumerate(Paths):
        
        path = Path + Thick + "/Res/"
        ## ----------------------------------- LET Read-in -----------------------------------------------------------
        # Only works if all input files have the same number of particle!!!!!
        LETHist = totalGRASLETHistos(path)

        NumberEntriesLETHist = np.sum(LETHist['entries'])
        TotalLET = np.sum(LETHist['mean'] * LETHist['value'])
        TotalLETError = 0
        for i in range(len(LETHist['value'])):
            TotalLETError += LETHist['error'][i] ** 2
        TotalLETError = np.sqrt(TotalLETError)
        TotalLETU = ufloat(TotalLET, TotalLETError)

        ### LET Histogram ###############
        fig, ax1 = plt.subplots(1)
        plt.bar(LETHist['lower'], LETHist['value'], width=LETHist['upper'] - LETHist['lower'], align='center', alpha=0.3)
        plt.errorbar(LETHist['mean'], LETHist['value'], LETHist['error'], fmt=' ', capsize=5, elinewidth=1,
                     capthick=1, label="LET Histogram")
        plt.plot([], [], label="SEE cross-section", color='C1')
        plt.yscale("log")
        plt.xscale("log")
        plt.grid(which='both')
        plt.title(DataName[P] + "  &  " + CrossectionName + " \nTotal LET = " + str(TotalLETU) + " MeV cm2 mg-1")
        plt.xlabel("LET [MeV cm2 mg-1]")
        ax1.legend(loc='center right')
        ax1.set_ylabel("Rate per LET bin [cm-2 s-1]", color='C0')
        ax1.tick_params(axis='y', colors='C0')

        ax2 = ax1.twinx()
        plt.plot(LETHist['lower'], f(LETHist['lower']), color='C1')
        ax2.set_ylabel("Cross Section [cm2 bit-1]", color='C1')
        plt.yscale("log")
        ax2.tick_params(axis='y', colors='C1')

        plt.savefig(path + "../" + DataName[P] + " " + CrossectionName + " LET.png", dpi=300, bbox_inches='tight')
        #plt.savefig("/l/triton_work/CARRINGTON/Histograms/" + DataName[P] + "/" + DataName[P] + " " + CrossectionName + " " + Thick + " LET-Hist.pdf",
        #    format='pdf', bbox_inches="tight")
        #plt.close('all')
        #plt.show()


        ### SEE rate ###############

        SEEHist = LETHist

        SEEHist['value'] = LETHist['value'] * f(LETHist['mean'])
        SEEHist['error'] = LETHist['error'] * f(LETHist['mean'])

        ################## Calculating total SEE Rate #####################################

        SEERate = np.sum(SEEHist['value'])
        SEERateError = 0
        for i in range(len(SEEHist['value'])):
            SEERateError += SEEHist['error'][i] ** 2
        SEERateError = np.sqrt(SEERateError)
        SEERateU = ufloat(SEERate, SEERateError)

        print("The total SEE rate is:", SEERateU, " s-1 bit-1 ")
        # print("The total number of Errors is:", SEERateU, " during the 3 year mission ")
        print("or:", SEERateU * 8e+9, " s-1 Gbyte-1 ")

        fig, ax1 = plt.subplots(1)
        plt.bar(SEEHist['lower'], SEEHist['value'], width=SEEHist['upper'] - SEEHist['lower'], align='edge',
                alpha=0.3)
        plt.errorbar(SEEHist['mean'], SEEHist['value'], SEEHist['error'], fmt=' ', capsize=5, elinewidth=1,
                     capthick=1, label="SEE Histogram")
        plt.plot([], [], label="SEE cross-section", color='C1')
        plt.yscale("log")
        plt.xscale("log")
        plt.grid(which='both')
        plt.title(DataName[P] + " + " + Thick + "Al  +  " + CrossectionName + " \nTotal SEERate = " + str(SEERateU) + " s-1 bit-1")
        #plt.title(DataName[P] + " + " + Thick + "Al  +  " + CrossectionName + " \nTotal Number of Errors = " + str(SEERateU) + " during the 3 year mission")
        plt.xlabel("LET [MeV cm2 mg-1]")
        ax1.legend(loc='center right')
        ax1.set_ylabel("SEE Rate per LET bin [s-1 bit-1]", color='C0')
        #ax1.set_ylabel("Errors per LET bin", color='C0')
        ax1.tick_params(axis='y', colors='C0')

        ax2 = ax1.twinx()
        plt.plot(LETHist['lower'], f(LETHist['lower']), color='C1')
        ax2.set_ylabel("Cross Section [cm2 bit-1]", color='C1')
        #ax2.set_ylabel("Cross Section of the receiver [cm2]", color='C1')
        plt.yscale("log")
        ax2.tick_params(axis='y', colors='C1')

        plt.savefig(path + "../" + DataName[P] + " " + CrossectionName + " SEE.png", dpi=300, bbox_inches='tight')
        #plt.savefig("/l/triton_work/CARRINGTON/Histograms/" + DataName[P] + "/" + DataName[P] + " " + CrossectionName + " " + Thick + " SEE-Hist.pdf",
        #    format='pdf', bbox_inches="tight")
        plt.close('all')
        #plt.show()


        CSVFile = open("/l/triton_work/LET/Foresail1-Hercules/SEERates.csv", 'a')
        List = (DataName[P], Thick, CrossectionName, SEERate, SEERateError)

        String = ', '.join(map(str, List))
        print(String)
        CSVFile.writelines(String + "\n")
        CSVFile.close()

