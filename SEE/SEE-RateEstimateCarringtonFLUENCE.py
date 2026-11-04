import os
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Dependencies.TotalFluenceHistograms import totalFluenceHistos
from Dependencies.HistogramPlots import histogramWithCrossSection
from Dependencies.SEERates import (PROTON_CROSS_SECTIONS, runFolders, thicknessFolders, usableHistogram,
                                   normaliseToRate, crossSection, foldCrossSection, seeRate, writeRateCSV)
from uncertainties import ufloat

Expected = 'blue' # Blue
PlusColor = 'C1'  # Orange
MinusColor = 'C2' # Green
LEOColor = 'C8'   # Yellow
MEOColor = 'C9'   # Turquoise
VAPColor = 'C3'   # Red
GEOColor = 'C7'   # Grey

Directory = "/scratch/work/fetzera1/GRAS/1Tile/Fluence_Histograms/CarringtonShielded/"

# Process all run folders in the directory, or only the ones given as command line arguments (all in one directory)
Directory, FolderList = runFolders(Directory, sys.argv[1:])

for parameters in PROTON_CROSS_SECTIONS:
    CrossectionName = parameters["name"]

    for F, Folder in enumerate(FolderList):

        # One line per thickness in the CSV of this run
        Rows = []

        for SubFolder in thicknessFolders(Directory, Folder):

            path = Directory + Folder + "/" + SubFolder + "/Res/"
            ## ----------------------------------- Fluence Read-in -----------------------------------------------------------
            ElectronHist, ProtonHist = totalFluenceHistos(path)

            if not usableHistogram(ProtonHist, path):
                continue

            # Normalise from 11 year fluence to flux
            ProtonHist = normaliseToRate(ProtonHist, Folder)

            # Each bin value is the content of the bin, so their sum is the total flux
            TotalFlux = np.sum(ProtonHist['value'])
            TotalFluxError = np.sqrt(np.sum(np.square(ProtonHist['error'])))
            TotalFluxU = ufloat(TotalFlux, TotalFluxError)

            Sigma = crossSection(ProtonHist['lower'], parameters)

            ### Flux Histogram ###############
            histogramWithCrossSection(
                ProtonHist, Sigma, path + "../" + Folder + " " + CrossectionName + " " + SubFolder + " Fluence-Hist.pdf",
                title=Folder + " " + CrossectionName + " " + SubFolder + " Al\nTotal Flux = " + str(TotalFluxU) + " cm-2 s-1",
                xlabel="Energy [MeV]", ylabel="Flux per Energy bin [cm-2 s-1]", label="Flux Histogram",
                sigmaLabel="SEU Cross Section", sigmaAxisLabel=CrossectionName + " Cross Section [cm2 bit-1]", legendLoc='upper left', capsize=2)


            ### SEE rate ###############

            SEEHist = foldCrossSection(ProtonHist, parameters)

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
                title=Folder + " " + CrossectionName + " " + SubFolder + " Al\nTotal SEU Rate = " + str(SEERateU) + " s-1 bit-1",
                xlabel="Energy [MeV]", ylabel="SEU Rate per Energy bin [s-1 bit-1]", label="SEU Rate Histogram",
                sigmaLabel="SEU Cross Section", sigmaAxisLabel=CrossectionName + " Cross Section [cm2 bit-1]", capsize=2)

        # Write the results of all thicknesses to a CSV file in the run folder
        writeRateCSV(Directory + Folder + "/" + Folder + " " + CrossectionName + " SEE-Rate.csv", Folder, CrossectionName, Rows)
