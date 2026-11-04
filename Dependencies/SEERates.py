"""Shared steps of the SEU rate scripts in SEE/.

Fold a LET or proton-energy histogram with an SEU cross section, compute the
total rate, write it next to the data and read it back for the rate plots.
A run folder holds one CSV per cross section, with one data line per shielding
thickness: <run>/<run> <cross section> SEE-Rate.csv.
"""
import csv
import glob
import os

import natsort
import numpy as np

# All cross sections are in cm2 bit-1. LET_CROSS_SECTIONS take LET in MeV cm2 mg-1,
# PROTON_CROSS_SECTIONS take proton energy in MeV.
#
# The Weibull form is A0 (1 - exp(-((x - L0) / W) ** S)), zero at and below the onset L0.
# It is the form of https://creme.isde.vanderbilt.edu/CREME-MC/help/weibull, but the
# parameters here are not in CREME-MC's input units (square microns per bit for LET).
#    A0 = limiting or plateau cross section;
#    L0 = onset parameter;
#    W = width parameter;
#    S = a dimensionless exponent.

# nanoXplore https://nanoxplore-wiki.atlassian.net/wiki/spaces/NAN/pages/46497810/NG-MEDIUM+Radiative+Test#Weibull-fitting
# Cypress: Garcia Alia et al. (2025), Table I
LET_CROSS_SECTIONS = (
    {
        "name": "NanoXplore SEU",
        "L0": 0.11,
        "W": 36.0,
        "S": 4.4,
        "A0": 5.2e-9,
    },
    {
        "name": "Cypress CY62167GE30-45ZXI",
        "L0": 0.1,
        "W": 70.0,
        "S": 1.2,
        "A0": 2.6e-7,
    },
)

# nanoXplore https://nanoxplore-wiki.atlassian.net/wiki/spaces/NAN/pages/46497810/NG-MEDIUM+Radiative+Test#Weibull-fitting
# Cypress: measured points extracted from Coronetti et al. (2021) Figure 3
PROTON_CROSS_SECTIONS = (
    {
        "name": "NanoXplore NG-Medium",
        "L0": 29.999,
        "W": 29.68,
        "S": 502E-3,
        "A0": 4.85E-016, # From NanoXplore wiki in cm2 bit-1
        #"A0": 1.52E-15, # From CERN test in cm2 bit-1
    },
    {
        "name": "Cypress CY62167GE30-45ZXI",
        "E_pts": np.array([
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
        ], dtype=float),
        "S_pts": np.array([
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
        ], dtype=float),
    },
)

MISSION_SECONDS = 4015 * 24 * 3600  # seconds in 11 years, the duration of the mission-fluence macros

# Runs simulated with mission-fluence macros (cm-2 per mission), matched in their folder names.
# All other runs use flux macros (cm-2 s-1) and are already rates. The GEO SAPPHIRE all-species
# campaign was simulated with the old mission-fluence macros. All macros generated since
# 2026-07-13, including the LEO SAPPHIRE set, are flux-normalised.
MISSION_FLUENCE_RUNS = ("-electron", "-solar-proton", "-trapped-proton", "-cosmic-proton", "-cosmic-iron",
                        "-solar-helium", "-solar-oxygen", "-solar-iron", "-solar-nickel", "GEO-SAPPHIRE")

RATE_CSV_HEADER = "Data,Shielding,Crossection,SEE_Rate,SEE_Error,Relative_SEE_Error,Entries_Contributing_To_SEE"


def runFolders(Directory, Arguments, Skip=("OldResults", "NoSEEs")):
    """Return (Directory, FolderList): all run folders in Directory, or only the ones given as arguments.

    Args:
        Directory (str): default directory holding the run folders, ending in '/'.
        Arguments (list): run folders from the command line, all in one directory. Empty uses the default.
        Skip (tuple): folder names in Directory that are not runs.

    Returns:
        tuple: the directory actually used, ending in '/', and the names of its run folders.
    """
    if Arguments:
        Paths = [os.path.abspath(p) for p in Arguments]
        Directory = os.path.dirname(Paths[0]) + "/"
        if any(os.path.dirname(p) + "/" != Directory for p in Paths):
            raise SystemExit("The run folders must all be in one directory")
        return Directory, [os.path.basename(p) for p in Paths]
    return Directory, sorted(f for f in os.listdir(Directory) if os.path.isdir(os.path.join(Directory, f)) and f not in Skip)


def thicknessFolders(Directory, Folder):
    """Return the shielding subfolders of a run in natural order (2mm before 16mm)."""
    Path = os.path.join(Directory, Folder)
    return natsort.natsorted(f for f in os.listdir(Path) if os.path.isdir(os.path.join(Path, f)))


def usableHistogram(Hist, Path):
    """Return True if a reader gave a histogram with consistent arrays; otherwise say why it is skipped."""
    if Hist is None:
        print("No histogram found in", Path, "-> skipping")
        return False
    Lengths = {len(Hist[Key]) for Key in ('lower', 'upper', 'mean', 'value', 'error', 'entries')}
    if len(Lengths) != 1:
        print("Inconsistent histogram lengths in", Path, "-> skipping")
        return False
    return True


def normaliseToRate(Hist, Folder):
    """Return a new histogram per second: runs with mission-fluence macros are divided by the mission duration."""
    Seconds = MISSION_SECONDS if any(Tag in Folder for Tag in MISSION_FLUENCE_RUNS) else 1.0
    RateHist = dict(Hist)
    RateHist['value'] = Hist['value'] / Seconds
    RateHist['error'] = Hist['error'] / Seconds
    return RateHist


def crossSection(X, Parameters):
    """Return the SEU cross section in cm2 bit-1 at X, LET in MeV cm2 mg-1 or proton energy in MeV.

    Parameters with E_pts and S_pts are measured points, interpolated log-log (zero below the
    lowest point, since extrapolating there would be very uncertain, and constant above the
    highest). All others are the Weibull parameters L0, W, S and A0.
    """
    if "E_pts" in Parameters:
        S = 10.0 ** np.interp(np.log10(X), np.log10(Parameters["E_pts"]), np.log10(Parameters["S_pts"]))
        return np.where(X < Parameters["E_pts"][0], 0.0, S)
    Result = np.zeros_like(X)
    Mask = X > Parameters["L0"]
    Result[Mask] = Parameters["A0"] * (1 - np.exp(-((X[Mask] - Parameters["L0"]) / Parameters["W"]) ** Parameters["S"]))
    return Result


def foldCrossSection(Hist, Parameters):
    """Return a new histogram of SEU rate per bin: each bin scaled by the cross section at its mean.

    Empty bins have a zero mean, which uses the lower edge instead to avoid log10(0) in the
    measured points. A zero-mean bin with a non-zero value would keep its cross section too.
    """
    Mean = np.where(Hist['mean'] == 0, Hist['lower'], Hist['mean'])
    Sigma = crossSection(Mean, Parameters)
    SEEHist = dict(Hist)
    SEEHist['mean'] = Mean
    SEEHist['value'] = Hist['value'] * Sigma
    SEEHist['error'] = Hist['error'] * Sigma
    return SEEHist


def seeRate(SEEHist):
    """Return the total rate, its error in quadrature, the relative error and the entries in contributing bins."""
    Rate = np.sum(SEEHist['value'])
    Error = np.sqrt(np.sum(np.square(SEEHist['error'])))
    Relative = Error / Rate if Rate != 0 else np.nan
    return Rate, Error, Relative, np.sum(SEEHist['entries'][SEEHist['value'] > 0])


def writeRateCSV(Path, Folder, CrossectionName, Rows):
    """Write the rates of one run and cross section, with one line per thickness.

    Args:
        Path (str): the CSV to write. Nothing is written if there are no rows.
        Folder (str): the run, written to the Data column.
        CrossectionName (str): the cross section, written to the Crossection column.
        Rows (list): tuples (thickness, rate, error, relative error, contributing entries).
    """
    if not Rows:
        print("No results for", Folder, "-> no CSV written")
        return
    with open(Path, 'w') as File:
        File.write(RATE_CSV_HEADER + "\n")
        for Row in Rows:
            File.write(','.join(map(str, (Folder, Row[0], CrossectionName) + tuple(Row[1:]))) + "\n")


def readRateCSVs(Patterns, CrossectionName):
    """Collect SEE-Rate CSVs into {dataset: {'Shielding': ..., 'SEE_Rate': ..., ...}} of numpy arrays.

    Args:
        Patterns (list): glob patterns of the CSVs. A file keeps its thicknesses in the order it has them.
        CrossectionName (str): the cross section every row must belong to.

    Returns:
        dict: per dataset, the arrays 'Shielding', 'SEE_Rate', 'SEE_Error' and 'Entries_Contributing_To_SEE'.

    Raises:
        ValueError: a row belongs to another cross section, or a dataset has the same thickness twice.
    """
    Dict = {}
    for Pattern in Patterns:
        for Path in glob.glob(Pattern):
            with open(Path, mode='r') as File:
                for Row in csv.DictReader(File):
                    if Row['Crossection'] != CrossectionName:
                        raise ValueError("Expected cross section " + CrossectionName + " in " + Path)
                    Dataset = Dict.setdefault(Row['Data'], {'Shielding': [], 'SEE_Rate': [], 'SEE_Error': [], 'Entries_Contributing_To_SEE': []})
                    if Row['Shielding'] in Dataset['Shielding']:
                        raise ValueError(Row['Data'] + " has " + Row['Shielding'] + " twice, again in " + Path)
                    Dataset['Shielding'].append(Row['Shielding'])
                    Dataset['SEE_Rate'].append(float(Row['SEE_Rate']))
                    Dataset['SEE_Error'].append(float(Row['SEE_Error']))
                    Dataset['Entries_Contributing_To_SEE'].append(float(Row['Entries_Contributing_To_SEE']))
    return {Name: {Key: np.array(Values) for Key, Values in Dataset.items()} for Name, Dataset in Dict.items()}
