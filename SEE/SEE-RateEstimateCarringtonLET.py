import os
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Dependencies.TotalLETHistos import totalLETHistos
from Dependencies.HistogramPlots import histogramWithCrossSection
from Dependencies.SEERates import (LET_CROSS_SECTIONS, runFolders, thicknessFolders, usableHistogram,
                                   normaliseToRate, crossSection, foldCrossSection, seeRate, writeRateCSV)
from uncertainties import ufloat

Directory = "/scratch/work/fetzera1/GRAS/1Tile/LET_Histograms/Carrington/"

# Process all run folders in the directory, or only the ones given as command line arguments (all in one directory)
Directory, FolderList = runFolders(Directory, sys.argv[1:])

for parameters in LET_CROSS_SECTIONS:
    CrossectionName = parameters["name"]

    for F, Folder in enumerate(FolderList):

        # One line per thickness in the CSV of this run
        Rows = []

        for SubFolder in thicknessFolders(Directory, Folder):

            path = Directory + Folder + "/" + SubFolder + "/Res/"
            ## ----------------------------------- LET Read-in -----------------------------------------------------------
            LETHist, _ = totalLETHistos(path)

            if not usableHistogram(LETHist, path):
                continue

            # Normalise from 11 year fluence to flux
            LETHist = normaliseToRate(LETHist, Folder)

            NumberEntriesLETHist = np.sum(LETHist['entries'])
            TotalLET = np.sum(LETHist['value'] * LETHist['mean'])
            TotalLETError = np.square(LETHist['error'] * LETHist['mean'])
            TotalLETError = np.sqrt(np.sum(TotalLETError))
            TotalLETU = ufloat(TotalLET, TotalLETError)

            Sigma = crossSection(LETHist['lower'], parameters)

            ### LET Histogram ###############
            histogramWithCrossSection(
                LETHist, Sigma, path + "../" + Folder + " " + CrossectionName + " " + SubFolder + " LET-Hist.pdf",
                title=Folder + " LET histogram behind " + SubFolder + " Al\nTotal LET rate = " + str(TotalLETU) + " MeV mg-1 s-1\nTotal entries = " + str(int(NumberEntriesLETHist)),
                xlabel="LET [MeV cm2 mg-1]", ylabel="Rate per LET bin [cm-2 s-1]", label="LET histogram",
                sigmaLabel="SEU cross section", sigmaAxisLabel=CrossectionName + " cross section [cm2 bit-1]", legendLoc='center right')


            ### SEE rate ###############

            SEEHist = foldCrossSection(LETHist, parameters)

            SEERate, SEERateError, RelativeError, EntriesContributingToSEE = seeRate(SEEHist)
            SEERateU = ufloat(SEERate, SEERateError)

            print("The total SEE rate is:", SEERateU, " s-1 bit-1 with ", int(EntriesContributingToSEE), " entries contributing to the SEE rate.")

            Rows.append((SubFolder, SEERate, SEERateError, RelativeError, EntriesContributingToSEE))

            if SEERate == 0:
                print("No SEEs in", path)
                continue

            ### SEE Histogram ###############
            histogramWithCrossSection(
                SEEHist, Sigma, path + "../" + Folder + " " + CrossectionName + " " + SubFolder + " SEE-Hist.pdf",
                title=Folder + " SEU rate histogram behind " + SubFolder + " Al\nTotal SEU rate = " + str(SEERateU) + " s-1 bit-1\nEntries contributing to SEU = " + str(int(EntriesContributingToSEE)),
                xlabel="LET [MeV cm2 mg-1]", ylabel="SEU rate per LET bin [s-1 bit-1]", label="SEU rate histogram",
                sigmaLabel="SEU cross section", sigmaAxisLabel=CrossectionName + " cross section [cm2 bit-1]", legendLoc='center left')

        # Write the results of all thicknesses to a CSV file in the run folder
        writeRateCSV(Directory + Folder + "/" + Folder + " " + CrossectionName + " SEE-Rate.csv", Folder, CrossectionName, Rows)
