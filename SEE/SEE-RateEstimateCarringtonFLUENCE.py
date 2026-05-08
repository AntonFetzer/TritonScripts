import os
import numpy as np
import matplotlib.pyplot as plt
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

# CrossectionName = "NanoXplore Proton SEU" 
CrossectionName = "Cypress CY62167GE30-45ZXI Proton SEU"

# nanoXplore https://nanoxplore-wiki.atlassian.net/wiki/spaces/NAN/pages/46497810/NG-MEDIUM+Radiative+Test#Weibull-fitting
if CrossectionName == "NanoXplore Proton SEU":
    L0 = 29.999
    W = 29.68
    S = 502E-3
    A0 = 4.85E-016 # From NanoXplore wiki in cm2 bit-1
    #A0 = 1.52E-15 # From CERN test in cm2 bit-1

    def f(Energy):
        """
        Calculate the rate estimate based on the Linear Energy Transfer (Fluence) values.

        Parameters:
        Energy (numpy.ndarray): Array of Energy values.

        Returns:
        numpy.ndarray: Array of rate estimates.
        """
        
        result = np.zeros_like(Energy)
        # Remove all values of Energy that are less than x0
        mask = Energy > L0
        # A * (1 - np.exp(-((x - x0) / W) ** s))
        result[mask] = A0 * (1 - np.exp(-((Energy[mask] - L0) / W) ** S))

        return result


if CrossectionName == "Cypress CY62167GE30-45ZXI Proton SEU":

    E_pts = np.array([
            0.6018484958761148,
            0.7017038286703828,
            0.7990318713533592,
            0.9045035097830326,
            1.0,
            1.1055789050944362,
            1.502875750260605,
            2.007077315996211,
            2.51188643150958,
            2.9986313485755667,
            4.004645731836126,
            5.011872336272722,
            40.04645731836127,
            80.37633606969334,
            123.68233972929845,
            164.20433388871635,
            183.69721161885397,
            # 1000.0,
        ], dtype=float)

    S_pts = np.array([
        4.0343096658867766e-11,
        1.3043213867190094e-10,
        1.068676240164437e-9,
        8.952061908580012e-10,
        1.2757607330635113e-9,
        4.2169650342858224e-10,
        5.7493054776844584e-11,
        4.923882631706752e-12,
        1.664002008471059e-12,
        6.281696115744896e-13,
        1.5229733635477114e-13,
        6.713099386829804e-14,
        1.0e-13,
        9.152473108773893e-14,
        8.193350588635765e-14,
        7.838460410979093e-14,
        7.666822074546222e-14,
        # 7.666822074546223e-14,
    ], dtype=float)

    def f(E_MeV: float | np.ndarray) -> float | np.ndarray:
        """
        Cypress SRAM proton SEU cross section (cm^2/bit) vs proton energy (MeV),
        using ONLY interpolation through the measured points from Coronetti et al. (2021) Figure 3.  

        Behavior:
        - Log-log piecewise linear interpolation between measured (E, sigma) points.
        - For E below the lowest point: Zero cross section (no extrapolation below the lowest measured energy, since that would be very uncertain and likely non-physical given the steep drop-off at low energy).
        - For E above the highest point: keep cross section constant at the highest-point value.

        Input:
        E_MeV: scalar or array-like proton energy in MeV

        Output:
        sigma in cm^2/bit (scalar if scalar input; numpy array otherwise)
        """
        # --- LOW-ENERGY PROTON (LEP) EMPIRICAL DATA ---
        # Extracted from Coronetti et al. (2021) Figure 3.  
        

        # Log-log interpolation
        logE_pts = np.log10(E_pts)
        logS_pts = np.log10(S_pts)
        logS = np.interp(np.log10(E_MeV), logE_pts, logS_pts)
        S = 10.0 ** logS

        # All S values below the lowest measured energy are set to zero (no extrapolation below the lowest measured energy, since that would be very uncertain and likely non-physical given the steep drop-off at low energy).
        S = np.where(E_MeV < E_pts[0], 0.0, S)

        return S



# Open the CSV file to write the results to
CSVFile = open(Directory + "/SEERates_" + CrossectionName + ".csv", 'w')
# Write the header to the file
header = "Data,Shielding,Crossection,SEE_Rate,SEE_Error,Relative_SEE_Error,Entries_Contributing_To_SEE"
CSVFile.write(header + "\n")


# Generate list of all folder names in the directory
FolderList = [f for f in os.listdir(Directory) if os.path.isdir(os.path.join(Directory, f))]

for F, Folder in enumerate(FolderList):

    # if Folder != "Carrington-SEP-Expected-Int":
    #     continue

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

        BinWidths = ProtonHist['upper'] - ProtonHist['lower']

        ### Fluence Histogram ###############
        fig, ax1 = plt.subplots(1)
        plt.bar(ProtonHist['lower'], ProtonHist['value'], width=BinWidths, align='edge', alpha=0.3) # pyright: ignore[reportCallIssue, reportArgumentType]
        plt.errorbar(ProtonHist['mean'], ProtonHist['value'], ProtonHist['error'], fmt=' ', capsize=2, elinewidth=1, # pyright: ignore[reportArgumentType] # pyright: ignore[reportCallIssue] # pyright: ignore[reportCallIssue] # type: ignore
                     capthick=1, label="Fluence Histogram")
        plt.plot([], [], label="SEU Cross Section", color='C1')
        plt.yscale("log")
        plt.xscale("log")
        plt.xlim(5e-1, 1.5e2)
        plt.grid(which='major')
        plt.title(Folder + " " + CrossectionName + " " + SubFolder + " Al\nTotal Fluence = " + str(TotalFluenceU) + " cm-2 s-1")
        plt.xlabel("Energy [MeV]")
        ax1.legend(loc='upper left')
        ax1.set_ylabel("Fluence per Energy bin [cm-2 s-1]", color='C0')
        ax1.tick_params(axis='y', colors='C0')

        ax2 = ax1.twinx()
        # plt.plot(E_pts, S_pts, '.-', color='C1')
        plt.plot(ProtonHist['lower'], f(ProtonHist['lower']), color='C1')
        ax2.set_ylabel(CrossectionName + " Cross Section [cm2 bit-1]", color='C1')
        plt.yscale("log")
        ax2.tick_params(axis='y', colors='C1')

        plt.savefig(path + "../" + Folder + " " + CrossectionName + " " + SubFolder + " Fluence-Hist.pdf", format='pdf', bbox_inches="tight")
        plt.close('all')
        #plt.show()


        ### SEE rate ###############

        SEEHist = ProtonHist

        # Some mean bin energies are zero.
        # Must replace with the lower bin edge to avoid zero-cross-section issue with the Weibull function.
        for i in range(len(SEEHist['mean'])):
            if SEEHist['mean'][i] == 0:
                SEEHist['mean'][i] = SEEHist['lower'][i] 

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
        plt.bar(SEEHist['lower'], SEEHist['value'], width=BinWidths, align='edge', alpha=0.3)
        plt.errorbar(SEEHist['mean'], SEEHist['value'], SEEHist['error'], fmt=' ', capsize=2, elinewidth=1,
                     capthick=1, label="SEU Rate Histogram")
        plt.plot([], [], label="SEU Cross Section", color='C1')
        plt.yscale("log")
        plt.xscale("log")
        plt.grid(which='major')
        plt.title(Folder + " " + CrossectionName + " " + SubFolder + " Al\nTotal SEU Rate = " + str(SEERateU) + " s-1 bit-1")
        plt.xlabel("Energy [MeV]")
        ax1.legend()
        ax1.set_ylabel("SEU Rate per Energy bin [s-1 bit-1]", color='C0')
        ax1.tick_params(axis='y', colors='C0')

        ax2 = ax1.twinx()
        # plt.plot(E_pts, S_pts, '.-', color='C1')
        plt.plot(SEEHist['lower'], f(SEEHist['lower']), color='C1')
        ax2.set_ylabel(CrossectionName + " Cross Section [cm2 bit-1]", color='C1')
        plt.yscale("log")
        plt.xlim(5e-1, 1.5e2)
        ax2.tick_params(axis='y', colors='C1')

        plt.savefig(path + "../" + Folder + " " + CrossectionName + " " + SubFolder + " SEE-Hist.pdf", format='pdf', bbox_inches="tight")
        plt.close('all')
        #plt.show()
        
CSVFile.close()
