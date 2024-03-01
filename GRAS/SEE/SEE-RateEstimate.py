import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from GRAS.Dependencies.TotalGRASLETHistos import totalGRASLETHistos
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


ThickList = ["0mm", "2mm", "4mm"]

FitParams = "SEU"

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
    result[mask] = A * (1 - np.exp(-((LET[mask] - x0) / W) ** s))

    return result


Paths = [
    "/l/triton_work/LET/Foresail1-Hercules/FS1-SolarProtons/"
    ]

DataName = [
    "FS1-SolarProtons"
    ]

lowerID = 0
upperID = 1
meanID = 2
valueID = 3
errorID = 4
entriesID = 5

for Thick in ThickList:

    for P, Path in enumerate(Paths):
        
        path = Path + Thick + "/Res/"
        ## ----------------------------------- LET Read-in -----------------------------------------------------------
        # Only works if all input files have the same number of particle!!!!!
        LETHist, EffHist = totalGRASLETHistos(path, "")

        C = 2330  # to convert from MeV/cm to Mev cm2 mg-1 in silicon with silicion density being 2.33 g/cm3

        LETHist[:, lowerID] = LETHist[:, lowerID] / C
        LETHist[:, upperID] = LETHist[:, upperID] / C
        LETHist[:, meanID] = LETHist[:, meanID] / C

        NumberEntriesLETHist = sum(LETHist[:, entriesID])
        TotalLET = sum(LETHist[:, meanID] * LETHist[:, valueID])
        TotalLETError = 0
        for i in range(len(LETHist[:, valueID])):
            TotalLETError += LETHist[i, errorID] ** 2
        TotalLETError = np.sqrt(TotalLETError)
        TotalLETU = ufloat(TotalLET, TotalLETError)

        ### LET Histogram ###############
        fig, ax1 = plt.subplots(1)
        plt.bar(LETHist[:, lowerID], LETHist[:, valueID], width=LETHist[:, upperID] - LETHist[:, lowerID], align='center', alpha=0.3)
        plt.errorbar(LETHist[:, meanID], LETHist[:, valueID], LETHist[:, errorID], fmt=' ', capsize=5, elinewidth=1,
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
        plt.plot(LETHist[:, lowerID], f(LETHist[:, lowerID]), color='C1')
        ax2.set_ylabel("Cross Section [cm2 bit-1]", color='C1')
        plt.yscale("log")
        ax2.tick_params(axis='y', colors='C1')

        #plt.savefig("/l/triton_work/CARRINGTON/Histograms/" + DataName[P] + "/" + DataName[P] + " " + CrossectionName + " " + Thick + " LET-Hist.pdf",
        #    format='pdf', bbox_inches="tight")
        #plt.close('all')
        plt.show()


        ### SEE rate ###############

        SEEHist = LETHist

        SEEHist[:, valueID] = LETHist[:, valueID] * f(LETHist[:, meanID])
        SEEHist[:, errorID] = LETHist[:, errorID] * f(LETHist[:, meanID])

        ################## Calculating total SEE Rate #####################################

        SEERate = np.sum(SEEHist[:, valueID])
        SEERateError = 0
        for i in range(len(SEEHist[:, valueID])):
            SEERateError += SEEHist[i, errorID] ** 2
        SEERateError = np.sqrt(SEERateError)
        SEERateU = ufloat(SEERate, SEERateError)

        print("The total SEE rate is:", SEERateU, " s-1 bit-1 ")
        # print("The total number of Errors is:", SEERateU, " during the 3 year mission ")
        print("or:", SEERateU * 8e+9, " s-1 Gbyte-1 ")

        fig, ax1 = plt.subplots(1)
        plt.bar(SEEHist[:, lowerID], SEEHist[:, valueID], width=SEEHist[:, upperID] - SEEHist[:, lowerID], align='edge',
                alpha=0.3)
        plt.errorbar(SEEHist[:, meanID], SEEHist[:, valueID], SEEHist[:, errorID], fmt=' ', capsize=5, elinewidth=1,
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
        plt.plot(LETHist[:, lowerID], f(LETHist[:, lowerID]), color='C1')
        ax2.set_ylabel("Cross Section [cm2 bit-1]", color='C1')
        #ax2.set_ylabel("Cross Section of the receiver [cm2]", color='C1')
        plt.yscale("log")
        ax2.tick_params(axis='y', colors='C1')

        #plt.savefig("/l/triton_work/CARRINGTON/Histograms/" + DataName[P] + "/" + DataName[P] + " " + CrossectionName + " " + Thick + " SEE-Hist.pdf",
        #    format='pdf', bbox_inches="tight")
        # plt.close('all')
        plt.show()


        CSVFile = open("/l/triton_work/CARRINGTON/SEERates.csv", 'a')
        List = (DataName[P], Thick, CrossectionName, SEERate, SEERateError)

        String = ', '.join(map(str, List))
        print(String)
        CSVFile.writelines(String + "\n")
        CSVFile.close()

