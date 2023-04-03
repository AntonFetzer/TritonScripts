import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from GRAS.Dependencies.TotalGRASLETHistos import totalGRASLETHistos
from uncertainties import ufloat

'''
The functional form of the Weibull is:
F(x) = A (1- exp{-[(x-x0)/W]s}) where
    x = effective LET in MeV-cm2/milligram;
    F(x) = SEE cross-section in square-microns/bit;
    A = limiting or plateau cross-section;
    x0 = onset parameter, such that F(x) = 0 for x < x0;
    W = width parameter;
    s = a dimensionless exponent.
https://creme.isde.vanderbilt.edu/CREME-MC/help/weibull
'''

Correctable = 1
ThickList = ["0mm", "1mm", "2mm", "4mm", "8mm", "16mm", ]

for Thick in ThickList:

    Paths = [#"/home/anton/Desktop/triton_work/LET/Carrington-SEP-Plus2Sigma-Int-With0/" + Thick + "/Res/",
         #"/home/anton/Desktop/triton_work/LET/Carrington-SEP-Expected-Int-With0/" + Thick + "/Res/",
         #"/home/anton/Desktop/triton_work/LET/Carrington-SEP-Minus2Sigma-Int-With0/" + Thick + "/Res/",
         #"/home/anton/Desktop/triton_work/LET/SEP2003-INTEGRAL-FluxBasedOnFluenceDividedBy24h/" + Thick + "/Res/",
         "/home/anton/Desktop/triton_work/LET/LETAP910MeV/" + Thick + "/Res/",
         "/home/anton/Desktop/triton_work/LET/ISS-LEO-Proton10MeV/" + Thick + "/Res/",
        "/home/anton/Desktop/triton_work/LET/LunarCosmic-H-Flux/" + Thick + "/Res/",
        "/home/anton/Desktop/triton_work/LET/LunarSEP10MeVFlux/" + Thick + "/Res/"]

    DataName = [#"Carrington SEP +2 Sigma",
            #"Carrington SEP EVT",
            #"Carrington SEP -2 Sigma",
            #"2003 SPE",
            "AP9 GTO trapped protons",
            "AP9 LEO trapped protons",
        "Cosmic Protons",
        "Solar Protons"
    ]

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

    LET = sp.Symbol('LET')

    f = A0 * (1 - sp.exp(-(LET / W) * S))

    f = sp.lambdify(LET, f, "numpy")

    '''
    Evals = np.geomspace(1, 80, num=100)
    
    plt.figure(0)
    
    plt.plot(Evals, f(Evals), label="Weibul")
    plt.legend()
    plt.grid(which='both')
    plt.yscale("log")
    plt.xscale("log")
    plt.title("LET crossection")
    plt.xlabel("LET [MeV cm2 mg-1]")
    plt.ylabel("Cross Section cm2/bit")
    '''

    for P, path in enumerate(Paths):
        ## ----------------------------------- LET Read-in -----------------------------------------------------------

        # Only works if all input files have the same number of particle!!!!!
        LETHist, EffHist = totalGRASLETHistos(path, "")

        lowerID = 0
        upperID = 1
        meanID = 2
        valueID = 3
        errorID = 4
        entriesID = 5

        C = 2330  # to convert from MeV/cm to Mev cm2 mg-1

        LETHist[:, lowerID] = LETHist[:, lowerID] / C
        LETHist[:, upperID] = LETHist[:, upperID] / C
        LETHist[:, meanID] = LETHist[:, meanID] / C

        NumberEntriesLETHist = sum(LETHist[:, entriesID])
        TotalLET = sum(LETHist[:, meanID] * LETHist[:, valueID])
        TotalLETError = 0
        for i in range(len(LETHist[:, valueID])):
            TotalLETError += LETHist[i, errorID] ** 2
        TotalLETError = np.sqrt(TotalLETError)
        TotalLET = ufloat(TotalLET, TotalLETError)

        '''
        ### LET Histogram ###############
        fig, ax1 = plt.subplots(1)
        # plt.figure(1)
        plt.bar(LETHist[:, lowerID], LETHist[:, valueID], width=LETHist[:, upperID] - LETHist[:, lowerID], align='edge', alpha=0.3)
        plt.errorbar(LETHist[:, meanID], LETHist[:, valueID], LETHist[:, errorID], fmt=' ', capsize=5, elinewidth=1, capthick=1, label="LET Histogram")
        plt.plot([], [], label="SEE cross-section", color='C1')
        plt.yscale("log")
        plt.xscale("log")
        plt.grid(which='both')
        plt.title(DataName[P] + "  &  " + CrossectionName + " \nTotal LET = " + str(TotalLET) + " cm-2 s-1")
        plt.xlabel("LET [MeV cm2 mg-1]")
        ax1.legend(loc='center right')
        ax1.set_ylabel("Rate per LET bin [cm-2 s-1]", color='C0')
        ax1.tick_params(axis='y', colors='C0')
        
        ax2 = ax1.twinx()
        plt.plot(LETHist[:, lowerID], f(LETHist[:, lowerID]), color='C1')
        ax2.set_ylabel("Cross Section [cm2 bit-1]", color='C1')
        plt.yscale("log")
        ax2.tick_params(axis='y', colors='C1')
        
        plt.savefig("/home/anton/Desktop/TritonPlots/Carrington/SEE/LET-Rate" + DataName[P] + " " + CrossectionName + ".pdf", format='pdf',
                    bbox_inches="tight")
        '''

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

        print("The total SEE rate is:", SEERate, " s-1 bit-1 ")
        print("or:", SEERateU * 8e+6, " s-1 Mbyte-1 ")

        fig, ax1 = plt.subplots(1)
        plt.bar(SEEHist[:, lowerID], SEEHist[:, valueID], width=SEEHist[:, upperID] - SEEHist[:, lowerID], align='edge',
                alpha=0.3)
        plt.errorbar(SEEHist[:, meanID], SEEHist[:, valueID], SEEHist[:, errorID], fmt=' ', capsize=5, elinewidth=1,
                     capthick=1, label="SEE Histogram")
        plt.plot([], [], label="SEE cross-section", color='C1')
        plt.yscale("log")
        plt.xscale("log")
        plt.grid(which='both')
        plt.title(DataName[P] + " + " + Thick + "Al  +  " + CrossectionName + " \nTotal SEERate = " + f"{SEERate:.2}" + " s-1 bit-1")
        plt.xlabel("LET [MeV cm2 mg-1]")
        ax1.legend(loc='center right')
        ax1.set_ylabel("SEE Rate per LET bin [s-1 bit-1]", color='C0')
        ax1.tick_params(axis='y', colors='C0')

        ax2 = ax1.twinx()
        plt.plot(LETHist[:, lowerID], f(LETHist[:, lowerID]), color='C1')
        ax2.set_ylabel("Cross Section [cm2 bit-1]", color='C1')
        plt.yscale("log")
        ax2.tick_params(axis='y', colors='C1')

        plt.savefig(
            "/home/anton/Desktop/TritonPlots/Luna/SEE/SEE-Rate" + DataName[P] + Thick + CrossectionName + ".svg",
            format='svg', bbox_inches="tight")


        CSVFile = open("/home/anton/Desktop/TritonPlots/Luna/SEE/SEERates.csv", 'a')
        List = (DataName[P], Thick, CrossectionName, SEERate, SEERateError)

        String = ', '.join(map(str, List))
        print(String)
        CSVFile.writelines(String + "\n")
        CSVFile.close()

        # plt.show()
