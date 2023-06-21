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
    A0 = limiting or plateau cross-section;
    x0 = onset parameter, such that F(x) = 0 for x < x0;
    W = width parameter;
    s = a dimensionless exponent.
https://creme.isde.vanderbilt.edu/CREME-MC/help/weibull
'''

HardSoft = 1
ThickList = ["0mm", "1mm", "2mm", "4mm", "8mm", "16mm"]

# NO Cutoff
# if HardSoft:
#     CrossectionName = "Hard reset cross section"
#     A = 9.014836795252201e-10
#     x0 = 1e-10
#     W = 0.02399012436035362
#     s = 1.0
# else:
#     CrossectionName = "Soft reset cross section"
#     A = 2.58987341772152e-09
#     x0 = 1e-10
#     W = 0.02399012436035362
#     s = 1.0

# With Cutoff
if HardSoft:
    CrossectionName = "Hard reset cross section"
    A = 1.00148367952522e-09
    x0 = 1e-10
    W = 0.002399012436035362
else:
    CrossectionName = "Soft reset cross section"
    A = 2.68987341772152e-09
    x0 = 1e-10
    W = 0.002399012436035362


LET = sp.Symbol('LET')
f_expr = A * (1 - sp.exp(-(LET - x0 / W)))
f_lamb = sp.lambdify(LET, f_expr, "numpy")


def f(array):
    return np.where(f_lamb(array) < 0, 0, f_lamb(array))


for Thick in ThickList:

    Paths = [
        "/home/anton/Desktop/triton_work/LET/Chess/Chess1-Proton-Mission-AP8/" + Thick + "/Res/",
        "/home/anton/Desktop/triton_work/LET/Chess/Chess1-SolarProton-Mission-Sapphire/" + Thick + "/Res/",
        "/home/anton/Desktop/triton_work/LET/Chess/Chess1-Electron-Mission-AE8/" + Thick + "/Res/",
        "/home/anton/Desktop/triton_work/LET/Chess/Chess1-CosmicProton-Mission-ISO/" + Thick + "/Res/",
        "/home/anton/Desktop/triton_work/LET/Chess/Chess1-CosmicFe-Mission-ISO/" + Thick + "/Res/"
    ]

    DataName = ["AP8", "SolarP", "AE8", "CosmicP", "CosmicFe"]

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
        # plt.figure(1)
        plt.bar(LETHist[:, lowerID], LETHist[:, valueID], width=LETHist[:, upperID] - LETHist[:, lowerID], align='edge',
                alpha=0.3)
        plt.errorbar(LETHist[:, meanID], LETHist[:, valueID], LETHist[:, errorID], fmt=' ', capsize=5, elinewidth=1,
                     capthick=1, label="LET Histogram")
        plt.plot([], [], label="SEE cross-section", color='C1')
        plt.yscale("log")
        plt.xscale("log")
        plt.grid(which='both')
        # plt.title(DataName + "  &  " + CrossectionName + " \nTotal LET = " + str(TotalLET) + " cm-2 s-1")
        plt.title(DataName[P] + "  &  " + CrossectionName + " \nTotal LET = " + str(TotalLETU) + " MeV cm2 mg-1")
        plt.xlabel("LET [MeV cm2 mg-1]")
        ax1.legend(loc='center right')
        # ax1.set_ylabel("Rate per LET bin [cm-2 s-1]", color='C0')
        ax1.set_ylabel("Total LET per LET bin [MeV cm2 mg-1]", color='C0')
        ax1.tick_params(axis='y', colors='C0')

        ax2 = ax1.twinx()
        plt.plot(LETHist[:, lowerID], f(LETHist[:, lowerID]), color='C1')
        # ax2.set_ylabel("Cross Section [cm2 bit-1]", color='C1')
        ax2.set_ylabel("Cross Section of the receiver [cm2]", color='C1')
        plt.yscale("log")
        ax2.tick_params(axis='y', colors='C1')
        plt.savefig("/home/anton/Desktop/triton_work/Chess1 GNSS SEE Analysis/Histograms/" + DataName[P] + "/" + DataName[P] + " " + CrossectionName + " " + Thick + " LET-Hist.pdf",
            format='pdf', bbox_inches="tight")
        # plt.show()
        plt.close('all')

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

        # print("The total SEE rate is:", SEERate, " s-1 bit-1 ")
        print("The total number of Errors is:", SEERateU, " during the 3 year mission ")
        # print("or:", SEERateU * 8e+6, " s-1 Mbyte-1 ")

        fig, ax1 = plt.subplots(1)
        plt.bar(SEEHist[:, lowerID], SEEHist[:, valueID], width=SEEHist[:, upperID] - SEEHist[:, lowerID], align='edge',
                alpha=0.3)
        plt.errorbar(SEEHist[:, meanID], SEEHist[:, valueID], SEEHist[:, errorID], fmt=' ', capsize=5, elinewidth=1,
                     capthick=1, label="SEE Histogram")
        plt.plot([], [], label="SEE cross-section", color='C1')
        plt.yscale("log")
        plt.xscale("log")
        plt.grid(which='both')
        # plt.title(DataName + " + " + Thick + "Al  +  " + CrossectionName + " \nTotal SEERate = " + f"{SEERate:.2}" + " s-1 bit-1")
        plt.title(
            DataName[P] + " + " + Thick + "Al  +  " + CrossectionName + " \nTotal Number of Errors = " + str(SEERateU) + " during the 3 year mission")
        plt.xlabel("LET [MeV cm2 mg-1]")
        ax1.legend(loc='center right')
        # ax1.set_ylabel("SEE Rate per LET bin [s-1 bit-1]", color='C0')
        ax1.set_ylabel("Errors per LET bin", color='C0')
        ax1.tick_params(axis='y', colors='C0')

        ax2 = ax1.twinx()
        plt.plot(LETHist[:, lowerID], f(LETHist[:, lowerID]), color='C1')
        # ax2.set_ylabel("Cross Section [cm2 bit-1]", color='C1')
        ax2.set_ylabel("Cross Section of the receiver [cm2]", color='C1')
        plt.yscale("log")
        ax2.tick_params(axis='y', colors='C1')

        plt.savefig("/home/anton/Desktop/triton_work/Chess1 GNSS SEE Analysis/Histograms/" + DataName[P] + "/" + DataName[P] + " " + CrossectionName + " " + Thick + " SEE-Hist.pdf",
            format='pdf', bbox_inches="tight")
        plt.close('all')

        CSVFile = open("/home/anton/Desktop/triton_work/Chess1 GNSS SEE Analysis/SEERates.csv", 'a')
        List = (DataName[P], Thick, CrossectionName, SEERate, SEERateError)

        String = ', '.join(map(str, List))
        print(String)
        CSVFile.writelines(String + "\n")
        CSVFile.close()

        # plt.show()
