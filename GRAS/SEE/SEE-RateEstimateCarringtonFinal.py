import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from GRAS.Dependencies.TotalLETHistos import totalLETHistos
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


Correctable = 0

ThickList = ["0mm", "1mm", "2mm", "4mm", "8mm", "16mm"]

if Correctable:
    CrossectionName = "LSRAM correctable SBU"
    L0 = 0.4
    W = 18
    S = 0.98
    A0 = 4.01e-9
else:
    CrossectionName = "LSRAM uncorrectable SBU"
    L0 = 0.4
    W = 1
    S = 0.4
    A0 = 5.5e-14


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
    mask = LET > L0
    # A * (1 - np.exp(-((x - x0) / W) ** s))
    result[mask] = A0 * (1 - np.exp(-((LET[mask] - L0) / W) ** S))

    return result

# List of folders to read in
# Carrington-SEP-Expected-Int-With0 Carrington-SEP-Plus2Sigma-Int-With0 GEO-AP9-mission GEO-SolarProton-mission  ISS-SolarProton-mission
# Template VAB-AP9-mission CarringtonElectronINTEGRALPowTabelated Carrington-SEP-Minus2Sigma-Int-With0 GEO-SolarProton-5minPeakFlux 
# ISS-AP9-mission VAB-AE9-mission  VAB-SolarProton-mission

Paths = [
    "/l/triton_work/LET_Histograms/Carrington/Carrington-SEP-Expected-Int-With0/",
    "/l/triton_work/LET_Histograms/Carrington/Carrington-SEP-Plus2Sigma-Int-With0/",
    "/l/triton_work/LET_Histograms/Carrington/Carrington-SEP-Minus2Sigma-Int-With0/",

    "/l/triton_work/LET_Histograms/Carrington/CarringtonElectronINTEGRALPowTabelated/",

    "/l/triton_work/LET_Histograms/Carrington/ISS-AP9-mission/",
    "/l/triton_work/LET_Histograms/Carrington/ISS-SolarProton-mission/",

    "/l/triton_work/LET_Histograms/Carrington/VAB-AP9-mission/",
    "/l/triton_work/LET_Histograms/Carrington/VAB-AE9-mission/",
    "/l/triton_work/LET_Histograms/Carrington/VAB-SolarProton-mission/",

    "/l/triton_work/LET_Histograms/Carrington/GEO-AP9-mission/",
    "/l/triton_work/LET_Histograms/Carrington/GEO-SolarProton-mission/",
    "/l/triton_work/LET_Histograms/Carrington/GEO-SolarProton-5minPeakFlux/",
    ]

DataName = [
    "Carrington SEP EVT",
    "Carrington SEP +2 Sigma",
    "Carrington SEP -2 Sigma",

    "Carrington Electron",

    "ISS AP9",
    "ISS Solar Proton",

    "VAB AP9",
    "VAB AE9",
    "VAB Solar Proton",

    "GEO AP9",
    "GEO Solar Proton",
    "GEO Solar Proton 5min Peak Flux",
    ]

for Thick in ThickList:

    for P, Path in enumerate(Paths):
        
        path = Path + Thick + "/Res/"
        ## ----------------------------------- LET Read-in -----------------------------------------------------------
        # Only works if all input files have the same number of particle!!!!!
        LETHist, _ = totalLETHistos(path)

        NumberEntriesLETHist = np.sum(LETHist['entries'])
        TotalLET = np.sum(LETHist['mean'] * LETHist['value'])
        TotalLETError = np.sum(np.square(LETHist['error']))
        TotalLETError = np.sqrt(TotalLETError)
        TotalLETU = ufloat(TotalLET, TotalLETError)

        ### LET Histogram ###############
        fig, ax1 = plt.subplots(1)
        plt.bar(LETHist['lower'], LETHist['value'], width=LETHist['upper'] - LETHist['lower'], align='edge', alpha=0.3)
        plt.errorbar(LETHist['mean'], LETHist['value'], LETHist['error'], fmt=' ', capsize=5, elinewidth=1,
                     capthick=1, label="LET Histogram")
        plt.plot([], [], label="SEE cross-section", color='C1')
        plt.yscale("log")
        plt.xscale("log")
        plt.grid(which='both')
        plt.title(DataName[P] + " " + Thick + "Al " + CrossectionName + "\nTotal LET = " + str(TotalLETU) + " MeV cm2 mg-1")
        plt.xlabel("LET [MeV cm2 mg-1]")
        ax1.legend(loc='center right')
        ax1.set_ylabel("Rate per LET bin [cm-2 s-1]", color='C0')
        ax1.tick_params(axis='y', colors='C0')

        ax2 = ax1.twinx()
        plt.plot(LETHist['lower'], f(LETHist['lower']), color='C1')
        ax2.set_ylabel("Cross Section [cm2 bit-1]", color='C1')
        plt.yscale("log")
        ax2.tick_params(axis='y', colors='C1')

        #plt.savefig(path + "../" + DataName[P] + " " + CrossectionName + " LET.png", dpi=300, bbox_inches='tight')
        plt.savefig(path + "../" + DataName[P] + " " + DataName[P] + " " + CrossectionName + " " + Thick + " LET-Hist.pdf",
            format='pdf', bbox_inches="tight")
        plt.close('all')
        #plt.show()


        ### SEE rate ###############

        SEEHist = LETHist
        # Scale the value and error by the cross section fuction
        SEEHist['value'] = LETHist['value'] * f(LETHist['mean'])
        SEEHist['error'] = LETHist['error'] * f(LETHist['mean'])

        ################## Calculating total SEE Rate #####################################

        SEERate = np.sum(SEEHist['value'])
        SEERateError = np.sum(np.square(SEEHist['error']))
        SEERateError = np.sqrt(SEERateError)
        SEERateU = ufloat(SEERate, SEERateError)

        # Calculate number of entries contributing to the total SEE rate
        # Only entries with LET > 0.4 MeV cm2 mg-1 are considered
        EntriesContributingToSEE = np.sum(LETHist['entries'][LETHist['mean'] > L0])

        print("The total SEE rate is:", SEERateU, " s-1 bit-1 ")
        # print("The total number of Errors is:", SEERateU, " during the 3 year mission ")
        print("or:", SEERateU * 8e+9, " s-1 Gbyte-1 ")

        fig, ax1 = plt.subplots(1)
        plt.bar(SEEHist['lower'], SEEHist['value'], width=SEEHist['upper'] - SEEHist['lower'], align='edge', alpha=0.3)
        plt.errorbar(SEEHist['mean'], SEEHist['value'], SEEHist['error'], fmt=' ', capsize=5, elinewidth=1,
                     capthick=1, label="SEE Histogram")
        plt.plot([], [], label="SEE cross-section", color='C1')
        plt.yscale("log")
        plt.xscale("log")
        plt.grid(which='both')
        plt.title(DataName[P] + " " + Thick + "Al " + CrossectionName + "\nTotal SEERate = " + str(SEERateU) + " s-1 bit-1")
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

        # plt.savefig(path + "../" + DataName[P] + " " + CrossectionName + " SEE.png", dpi=300, bbox_inches='tight')
        plt.savefig(path + "../" + DataName[P] + " " + DataName[P] + " " + CrossectionName + " " + Thick + " SEE-Hist.pdf",
            format='pdf', bbox_inches="tight")
        plt.close('all')
        #plt.show()

        
        CSVFile = open("/l/triton_work/LET_Histograms/Carrington/SEERates.csv", 'a')
        List = (DataName[P], Thick, CrossectionName, SEERate, SEERateError, EntriesContributingToSEE)

        String = ', '.join(map(str, List))
        print(String)
        CSVFile.writelines(String + "\n")
        CSVFile.close()
        

