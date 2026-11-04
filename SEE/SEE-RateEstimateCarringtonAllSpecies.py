"""Estimate SEE rates for the Carrington all-species datasets.

Adapted from SEE-RateEstimateCarringtonLET.py.  By default, all folders with
All-Species in their name in LET_Histograms are analysed.  Alternative dataset
directories can be supplied as command line arguments.

Every species folder is a run and gets one SEE-Rate CSV per cross section.  No
histogram plots are made: with about 2900 thicknesses they would take hours.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Dependencies.TotalLETHistos import totalLETHistos
from Dependencies.SEERates import (LET_CROSS_SECTIONS, runFolders, thicknessFolders, usableHistogram,
                                   normaliseToRate, foldCrossSection, seeRate, writeRateCSV)
from uncertainties import ufloat

BaseDirectory = "/scratch/work/fetzera1/GRAS/1Tile/LET_Histograms/"

# Generate list of all folders with All-Species in their name
SetList = [BaseDirectory + f_ for f_ in sorted(os.listdir(BaseDirectory)) if "All-Species" in f_ and os.path.isdir(BaseDirectory + f_)]

# Alternatively process only the folders given as command line arguments
if len(sys.argv) > 1:
    SetList = [os.path.abspath(p) for p in sys.argv[1:]]

for parameters in LET_CROSS_SECTIONS:
    CrossectionName = parameters["name"]

    for SetFolder in SetList:

        Directory, FolderList = runFolders(SetFolder + "/", [])

        for F, Folder in enumerate(FolderList):

            # One line per thickness in the CSV of this run
            Rows = []

            for SubFolder in thicknessFolders(Directory, Folder):

                path = Directory + Folder + "/" + SubFolder + "/Res/"
                ## ----------------------------------- LET Read-in -----------------------------------------------------------
                LETHist, _ = totalLETHistos(path)

                if not usableHistogram(LETHist, path):
                    continue

                # Normalise from 11 year fluence to flux where needed, then scale by the cross section
                SEEHist = foldCrossSection(normaliseToRate(LETHist, Folder), parameters)

                SEERate, SEERateError, RelativeError, EntriesContributingToSEE = seeRate(SEEHist)

                print("The total SEE rate is:", ufloat(SEERate, SEERateError), " s-1 bit-1 with ", int(EntriesContributingToSEE), " entries contributing to the SEE rate.")

                Rows.append((SubFolder, SEERate, SEERateError, RelativeError, EntriesContributingToSEE))

            # Write the results of all thicknesses to a CSV file in the run folder
            writeRateCSV(Directory + Folder + "/" + Folder + " " + CrossectionName + " SEE-Rate.csv", Folder, CrossectionName, Rows)
