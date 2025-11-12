import os

import numpy as np
import matplotlib.pyplot as plt
from Dependencies.TotalLETHistos import totalLETHistos
from uncertainties import ufloat
import os
import natsort

Expected = 'blue' # Blue
PlusColor = 'C1'  # Orange
MinusColor = 'C2' # Green
LEOColor = 'C8'   # Yellow
MEOColor = 'C9'   # Turquoise
VAPColor = 'C3'   # Red
GEOColor = 'C7'   # Grey

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
# Directory = "/scratch/work/fetzera1/LET_Histograms/Carrington/"
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
CrossectionName = "NanoXplore SEU"  
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
CSVFile = open(Directory + "/SEERatesnanoXplore.csv", 'w')
# Write the header to the file
header = "Data,Shielding,Crossection,SEE_Rate,SEE_Error,Relative_SEE_Error,Entries_Contributing_To_SEE"
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
        
        path = Directory + Folder + "/" + SubFolder + "/Res/"
        ## ----------------------------------- LET Read-in -----------------------------------------------------------
        LETHist, _ = totalLETHistos(path)

        # Check if the LET histogram is None (no files found)
        if LETHist is None:
            print("No LET histogram found in", path, "-> skipping")
            continue

        # Ensure the histogram arrays returned have consistent binning/length.
        # The script requires 'lower','upper','mean','value','error','entries' to be the same length.
        required_keys = ['lower', 'upper', 'mean', 'value', 'error', 'entries']
        lengths = [len(LETHist[k]) for k in required_keys]
        if not all(l == lengths[0] for l in lengths):
            print("Inconsistent histogram lengths in", path, "-> skipping")
            continue

        # Normalise from 11 year fluence to flux
        if "-electron" in Folder or "-solar-proton" in Folder or "-trapped-proton" in Folder or "-cosmic-proton" in Folder or "-cosmic-iron" in Folder:
            NormalisationFactor = 4015 * 24 * 3600  # seconds in 11 years
            LETHist['value'] = LETHist['value'] / NormalisationFactor
            LETHist['error'] = LETHist['error'] / NormalisationFactor

        NumberEntriesLETHist = np.sum(LETHist['entries'])
        TotalLET = np.sum(LETHist['value'] * LETHist['mean'])
        TotalLETError = np.square(LETHist['error'] * LETHist['mean'])
        TotalLETError = np.sqrt(np.sum(TotalLETError))
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
        plt.title(Folder + " LET Histogram behind " + SubFolder + " Al\nTotal LET rate = " + str(TotalLETU) + " MeV cm2 mg-1 s-1\nTotal Entries = " + str(int(NumberEntriesLETHist)))
        plt.xlabel("LET [MeV cm2 mg-1]")
        ax1.legend(loc='center right')
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

        SEERateError = np.square(SEEHist['error'])
        SEERateError = np.sqrt(np.sum(SEERateError))

        SEERateU = ufloat(SEERate, SEERateError)
        if SEERate != 0:
            RelativeError = SEERateError / SEERate
        else:
            RelativeError = np.nan

        # Calculate number of entries contributing to the total SEE rate
        EntriesContributingToSEE = 0  # Initialize counter
        for i in range(len(SEEHist['value'])):
            if SEEHist['value'][i] > 0:
                EntriesContributingToSEE += SEEHist['entries'][i]

        print("The total SEE rate is:", SEERateU, " s-1 bit-1 with ", int(EntriesContributingToSEE), " entries contributing to the SEE rate.")
        #print("or:", SEERateU * 8e+9, " s-1 Gbyte-1 ")

        # Write the results to the CSV file
        List = (Folder, SubFolder, CrossectionName, SEERate, SEERateError, RelativeError, EntriesContributingToSEE)

        String = ','.join(map(str, List))
        #print(String)
        CSVFile.writelines(String + "\n")

        if SEERate == 0:
            print("No SEEs in", path)
            continue

        ### SEE Histogram ###############

        fig, ax1 = plt.subplots(1)
        plt.bar(SEEHist['lower'], SEEHist['value'], width=SEEHist['upper'] - SEEHist['lower'], align='edge', alpha=0.3)
        plt.errorbar(SEEHist['mean'], SEEHist['value'], SEEHist['error'], fmt=' ', capsize=5, elinewidth=1,
                     capthick=1, label="SEU Rate Histogram")
        plt.plot([], [], label="SEU Cross Section", color='C1')
        plt.yscale("log")
        plt.xscale("log")
        plt.grid(which='major')
        plt.title(Folder + " SEU Rate Histogram behind " + SubFolder + " Al\nTotal SEU Rate = " + str(SEERateU) + " s-1 bit-1")
        plt.xlabel("LET [MeV cm2 mg-1]")
        ax1.legend(loc='center left')
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
