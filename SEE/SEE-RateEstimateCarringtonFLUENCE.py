import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from Dependencies.TotalFluenceHistograms import totalFluenceHistos
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
    x = effective Fluence in MeV-cm2/milligram;
    F(x) = SEU cross-section in square-microns/bit;
    A0 = limiting or plateau cross-section;
    x0 = onset parameter, such that F(x) = 0 for x < x0;
    W = width parameter;
    s = a dimensionless exponent.
https://creme.isde.vanderbilt.edu/CREME-MC/help/weibull
'''

Directory = "/l/triton_work/Fluence_Histograms/CarringtonShielded/"

# nanoXplore https://nanoxplore-wiki.atlassian.net/wiki/spaces/NAN/pages/46497810/NG-MEDIUM+Radiative+Test#Weibull-fitting
CrossectionName = "NanoXplore Proton SEU"  
L0 = 29.999
W = 29.68
S = 502E-3
A0 = 4.85E-016 # From NanoXplore wiki in cm2 bit-1
#A0 = 1.52E-15 # From CERN test in cm2 bit-1


def f(Fluence):
    """
    Calculate the rate estimate based on the Linear Energy Transfer (Fluence) values.

    Parameters:
    Fluence (numpy.ndarray): Array of Fluence values.

    Returns:
    numpy.ndarray: Array of rate estimates.
    """
    
    result = np.zeros_like(Fluence)
    # Remove all values of Fluence that are less than x0
    mask = Fluence > L0
    # A * (1 - np.exp(-((x - x0) / W) ** s))
    result[mask] = A0 * (1 - np.exp(-((Fluence[mask] - L0) / W) ** S))

    return result


# Open the CSV file to write the results to
CSVFile = open(Directory + "/SEERatesnanoXplorePROTON.csv", 'w')
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
        ## ----------------------------------- Fluence Read-in -----------------------------------------------------------
        ElectronHist, ProtonHist = totalFluenceHistos(path)

        # Normalise from 11 year fluence to flux
        if "-solar-proton" in Folder or "-trapped-proton" in Folder or "-cosmic-proton" in Folder:
            NormalisationFactor = 4015 * 24 * 3600  # seconds in 11 years
            ProtonHist['value'] = ProtonHist['value'] / NormalisationFactor
            ProtonHist['error'] = ProtonHist['error'] / NormalisationFactor

        NumberEntriesProtonHist = np.sum(ProtonHist['entries'])
        TotalFluence = np.sum(ProtonHist['value'] * ProtonHist['mean'])
        TotalFluenceError = np.square(ProtonHist['error'] * ProtonHist['mean'])
        TotalFluenceError = np.sqrt(np.sum(TotalFluenceError))
        TotalFluenceU = ufloat(TotalFluence, TotalFluenceError)

        ### Fluence Histogram ###############
        fig, ax1 = plt.subplots(1)
        plt.bar(ProtonHist['lower'], ProtonHist['value'], width=ProtonHist['upper'] - ProtonHist['lower'], align='edge', alpha=0.3) # pyright: ignore[reportCallIssue, reportArgumentType]
        plt.errorbar(ProtonHist['mean'], ProtonHist['value'], ProtonHist['error'], fmt=' ', capsize=5, elinewidth=1, # pyright: ignore[reportArgumentType] # pyright: ignore[reportCallIssue] # pyright: ignore[reportCallIssue] # type: ignore
                     capthick=1, label="Fluence Histogram")
        plt.plot([], [], label="SEU Cross Section", color='C1')
        plt.yscale("log")
        plt.xscale("log")
        plt.grid(which='major')
        plt.title(Folder + " " + CrossectionName + " " + SubFolder + " Al\nTotal Fluence = " + str(TotalFluenceU) + " cm-2 s-1")
        plt.xlabel("Energy [MeV]")
        ax1.legend(loc='center right')
        ax1.set_ylabel("Fluence per Energy bin [cm-2 s-1]", color='C0')
        ax1.tick_params(axis='y', colors='C0')

        ax2 = ax1.twinx()
        plt.plot(ProtonHist['lower'], f(ProtonHist['lower']), color='C1')
        ax2.set_ylabel(CrossectionName + " Cross Section [cm2 bit-1]", color='C1')
        plt.yscale("log")
        ax2.tick_params(axis='y', colors='C1')

        plt.savefig(path + "../" + Folder + " " + CrossectionName + " " + SubFolder + " Fluence-Hist.pdf", format='pdf', bbox_inches="tight")
        plt.close('all')
        #plt.show()


        ### SEE rate ###############

        SEEHist = ProtonHist
        # Scale the value and error by the cross section fuction
        SEEHist['value'] = ProtonHist['value'] * f(ProtonHist['mean'])
        SEEHist['error'] = ProtonHist['error'] * f(ProtonHist['mean'])

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

        # Write the results to the CSV file
        List = (Folder, SubFolder, CrossectionName, SEERate, SEERateError, RelativeError, EntriesContributingToSEE)

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
        plt.xlabel("Energy [MeV]")
        ax1.legend(loc='center right')
        ax1.set_ylabel("SEU Rate per Energy bin [s-1 bit-1]", color='C0')
        ax1.tick_params(axis='y', colors='C0')

        ax2 = ax1.twinx()
        plt.plot(ProtonHist['lower'], f(ProtonHist['lower']), color='C1')
        ax2.set_ylabel(CrossectionName + " Cross Section [cm2 bit-1]", color='C1')
        plt.yscale("log")
        ax2.tick_params(axis='y', colors='C1')

        plt.savefig(path + "../" + Folder + " " + CrossectionName + " " + SubFolder + " SEE-Hist.pdf", format='pdf', bbox_inches="tight")
        plt.close('all')
        #plt.show()
        
CSVFile.close()
