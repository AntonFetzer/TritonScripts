import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from GRAS.Dependencies.TotalLETHistos import totalLETHistos
from uncertainties import ufloat
import os
import natsort

'''
The functional form of the Weibull is:
F(x) = A (1- exp{-[(x-x0)/W] ** s}) where
    x = effective LET in MeV-cm2/milligram;
    F(x) = SEU cross-section in square-microns/bit;
    A0 = limiting or plateau cross-section;
    x0 = onset parameter, such that F(x) = 0 for x < x0;
    W = width parameter;
    s = a dimensionless exponent.
https://creme.isde.vanderbilt.edu/CREME-MC/help/weibull
'''

Directory = "/l/triton_work/LET_Histograms/Carrington/"
#Correctable = 0

# if Correctable:
#     CrossectionName = "LSRAM correctable SEU"
#     L0 = 0.4
#     W = 18
#     S = 0.98
#     A0 = 4.01e-9
# else:
#     CrossectionName = "LSRAM uncorrectable SEU" # The uncorrectable parameters are only estimated upper bounds. This means the actual cross section is unknown.
#     L0 = 0.4
#     W = 1
#     S = 0.4
#     A0 = 5.5e-14
# else:   # nanoXplore https://nanoxplore-wiki.atlassian.net/wiki/spaces/NAN/pages/46497810/NG-MEDIUM+Radiative+Test#Weibull-fitting
CrossectionName = "nanoXplore SEU"  
L0 = 0.11 
W = 36
S = 4.4
A0 = 5.2E-09


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


# Open the CSV file to write the results to
# if Correctable:
#     CSVFile = open("/l/triton_work/LET_Histograms/Carrington/SEERatesCorrectable.csv", 'w')
# if not Correctable:
#     CSVFile = open("/l/triton_work/LET_Histograms/Carrington/SEERatesUncorrectable.csv", 'w')
CSVFile = open("/l/triton_work/LET_Histograms/Carrington/SEERatesnanoXplore.csv", 'w')
# Write the header to the file
header = "Data,Shielding,Crossection,SEE_Rate,SEE_Error,EntriesContributingToSEE"
CSVFile.write(header + "\n")


# Generate list of all folder names in the directory
FolderList = [f for f in os.listdir(Directory) if os.path.isdir(os.path.join(Directory, f))]

for F, Folder in enumerate(FolderList):

    if Folder == "NoSEEs":
        continue

    # Generate list of all subfolder names in the Folder
    SubFolderList = [f for f in os.listdir(Directory + Folder) if os.path.isdir(os.path.join(Directory + Folder, f))]
    
    # Sort the SubFolderList
    SubFolderList = natsort.natsorted(SubFolderList)

    for SubFolder in SubFolderList:
        # Print a dot for each folder processed as a progress indicator
        print(".", end='', flush=True)
        
        path = Directory + Folder + "/" + SubFolder + "/Res/"
        ## ----------------------------------- LET Read-in -----------------------------------------------------------
        # Only works if all input files have the same number of particle!!!!!
        try:
            LETHist, _ = totalLETHistos(path)
        except:
            print("Error in", path)
            continue

        # Convert Fluences to Fluxes for comparison.
        # Path names sets containing Carrington-SEP and CarringtonElectron are fluxes in [cm-2 s-1]
        # Datasets containing "mission" are 30 day fluences in [cm-2] (per 30 days)
        if "mission" in Folder:
            LETHist['value'] = LETHist['value'] / ( 30 * 24 * 3600 ) 
            LETHist['error'] = LETHist['error'] / ( 30 * 24 * 3600 )


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
        plt.plot([], [], label="SEU Cross Section", color='C1')
        plt.yscale("log")
        plt.xscale("log")
        plt.grid(which='major')
        plt.title(Folder + " " + CrossectionName + " " + SubFolder + " Al\nTotal LET = " + str(TotalLETU) + " MeV cm2 mg-1")
        plt.xlabel("LET [MeV cm2 mg-1]")
        ax1.legend(loc='lower left')
        ax1.set_ylabel("Rate per LET bin [cm-2 s-1]", color='C0')
        ax1.tick_params(axis='y', colors='C0')

        ax2 = ax1.twinx()
        plt.plot(LETHist['lower'], f(LETHist['lower']), color='C1')
        ax2.set_ylabel(CrossectionName + " Cross Section [cm2 bit-1]", color='C1')
        plt.yscale("log")
        ax2.tick_params(axis='y', colors='C1')

        plt.savefig(path + "../" + Folder + " " + CrossectionName + " " + SubFolder + " LET-Hist.pdf", format='pdf', bbox_inches="tight")
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

        #print("The total SEE rate is:", SEERateU, " s-1 bit-1 ")
        #print("or:", SEERateU * 8e+9, " s-1 Gbyte-1 ")

        # Write the results to the CSV file
        List = (Folder, SubFolder, CrossectionName, SEERate, SEERateError, EntriesContributingToSEE)

        String = ','.join(map(str, List))
        #print(String)
        CSVFile.writelines(String + "\n")

        if SEERate == 0:
            print("No SEEs in", path)
            continue

        fig, ax1 = plt.subplots(1)
        plt.bar(SEEHist['lower'], SEEHist['value'], width=SEEHist['upper'] - SEEHist['lower'], align='edge', alpha=0.3)
        plt.errorbar(SEEHist['mean'], SEEHist['value'], SEEHist['error'], fmt=' ', capsize=5, elinewidth=1,
                     capthick=1, label="SEU Rate Histogram")
        plt.plot([], [], label="SEU Cross Section", color='C1')
        plt.yscale("log")
        plt.xscale("log")
        plt.grid(which='major')
        plt.title(Folder + " " + CrossectionName + " " + SubFolder + " Al\nTotal SEU Rate = " + str(SEERateU) + " s-1 bit-1")
        plt.xlabel("LET [MeV cm2 mg-1]")
        ax1.legend(loc='lower left')
        ax1.set_ylabel("SEU Rate per LET bin [s-1 bit-1]", color='C0')
        ax1.tick_params(axis='y', colors='C0')

        ax2 = ax1.twinx()
        plt.plot(LETHist['lower'], f(LETHist['lower']), color='C1')
        ax2.set_ylabel(CrossectionName + " Cross Section [cm2 bit-1]", color='C1')
        plt.yscale("log")
        ax2.tick_params(axis='y', colors='C1')

        plt.savefig(path + "../" + Folder + " " + CrossectionName + " " + SubFolder + " SEE-Hist.pdf", format='pdf', bbox_inches="tight")
        plt.close('all')
        #plt.show()
        
CSVFile.close()
