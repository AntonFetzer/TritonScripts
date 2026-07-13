import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
from Dependencies.TotalLETHistos import totalLETHistos

'''
Compares two GRAS LET histogram simulation results with each other using the
same algorithm as SEE-RateEstimateCarringtonFinal.py:
merge the per-job csv histograms, normalise the 11 year fluence to flux,
compute the total LET rate and fold the spectrum with a Weibull SEU cross
section to get the total SEE rate.

Here: CREME96 solar iron at GEO sampled directly from the spectrum vs
sampled flat in log energy with biasing weights (FlatLog).
'''

Directory = os.path.expanduser("~/triton_work/GRAS/LET_Histograms/Carrington-CREME96-All-Species/")

DirectFolder = "GEO-CREME96-Solar-Iron"
FlatLogFolder = "GEO-CREME96-Solar-Iron-FlatLog"
SubFolder = "1mm"

DirectColor = 'C0'
FlatLogColor = 'C1'

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
Crossections = [
    # Cypress SRAM
    {'Name': "Cypress CY62167GE30-45ZXI", 'L0': 0.1, 'W': 70, 'S': 1.2, 'A0': 2.6E-07},
    # nanoXplore https://nanoxplore-wiki.atlassian.net/wiki/spaces/NAN/pages/46497810/NG-MEDIUM+Radiative+Test#Weibull-fitting
    {'Name': "NanoXplore SEU", 'L0': 0.11, 'W': 36, 'S': 4.4, 'A0': 5.2E-09},
]


def f(LET, X):
    """Weibull SEU cross section with parameters X, see https://creme.isde.vanderbilt.edu/CREME-MC/help/weibull"""
    result = np.zeros_like(LET)
    mask = LET > X['L0']
    result[mask] = X['A0'] * (1 - np.exp(-((LET[mask] - X['L0']) / X['W']) ** X['S']))
    return result


def readAndNormalise(Folder):
    """Merges the csv histograms of one folder and normalises fluence to flux as in SEE-RateEstimateCarringtonFinal.py."""
    path = Directory + Folder + "/" + SubFolder + "/Res/"

    LETHist, _ = totalLETHistos(path)

    if LETHist is None:
        sys.exit("No LET histogram found in " + path)

    required_keys = ['lower', 'upper', 'mean', 'value', 'error', 'entries']
    lengths = [len(LETHist[k]) for k in required_keys]
    if not all(l == lengths[0] for l in lengths):
        sys.exit("Inconsistent histogram lengths in " + path)

    # Normalise from 11 year fluence to flux
    if "-solar-iron" in Folder.lower():
        NormalisationFactor = 4015 * 24 * 3600  # seconds in 11 years
        LETHist['value'] = LETHist['value'] / NormalisationFactor
        LETHist['error'] = LETHist['error'] / NormalisationFactor

    Result = {'hist': LETHist}
    Result['entries'] = np.sum(LETHist['entries'])

    Result['TotalLET'] = np.sum(LETHist['value'] * LETHist['mean'])
    Result['TotalLETError'] = np.sqrt(np.sum(np.square(LETHist['error'] * LETHist['mean'])))

    return Result


def seeRate(Result, X):
    """Folds the LET spectrum with the Weibull cross section X and returns the SEE rate results."""
    LETHist = Result['hist']
    SEE = {}
    SEE['value'] = LETHist['value'] * f(LETHist['mean'], X)
    SEE['error'] = LETHist['error'] * f(LETHist['mean'], X)

    SEE['SEERate'] = np.sum(SEE['value'])
    SEE['SEERateError'] = np.sqrt(np.sum(np.square(SEE['error'])))
    SEE['RelativeError'] = SEE['SEERateError'] / SEE['SEERate'] if SEE['SEERate'] != 0 else np.nan

    SEE['EntriesContributingToSEE'] = int(np.sum(LETHist['entries'][SEE['value'] > 0]))

    return SEE


Direct = readAndNormalise(DirectFolder)
FlatLog = readAndNormalise(FlatLogFolder)

DirectHist = Direct['hist']
FlatLogHist = FlatLog['hist']

### Print and save the comparison table ###############
CSVFile = open(Directory + "Compare_FlatLog_vs_Direct_" + SubFolder + ".csv", 'w')
header = "Crossection,Quantity,Direct,FlatLog,Ratio_FlatLog/Direct"
print("\n" + header.replace(",", "\t"))
CSVFile.write(header + "\n")

Rows = [
    ("", "Total_entries", Direct['entries'], FlatLog['entries']),
    ("", "Total_LET_rate_MeVcm2mg-1s-1", Direct['TotalLET'], FlatLog['TotalLET']),
    ("", "Total_LET_rate_error", Direct['TotalLETError'], FlatLog['TotalLETError']),
]

SEEResults = {}
for X in Crossections:
    DirectSEE = seeRate(Direct, X)
    FlatLogSEE = seeRate(FlatLog, X)
    SEEResults[X['Name']] = (DirectSEE, FlatLogSEE)

    CombinedSigma = np.sqrt(DirectSEE['SEERateError'] ** 2 + FlatLogSEE['SEERateError'] ** 2)
    Deviation = (FlatLogSEE['SEERate'] - DirectSEE['SEERate']) / CombinedSigma if CombinedSigma != 0 else np.nan

    Rows += [
        (X['Name'], "SEE_rate_s-1bit-1", DirectSEE['SEERate'], FlatLogSEE['SEERate']),
        (X['Name'], "SEE_rate_error", DirectSEE['SEERateError'], FlatLogSEE['SEERateError']),
        (X['Name'], "Relative_SEE_error", DirectSEE['RelativeError'], FlatLogSEE['RelativeError']),
        (X['Name'], "Entries_contributing_to_SEE", DirectSEE['EntriesContributingToSEE'], FlatLogSEE['EntriesContributingToSEE']),
        (X['Name'], "SEE_rate_deviation_sigma", Deviation, np.nan),
    ]

for Crossection, Name, D, F in Rows:
    if Name == "SEE_rate_deviation_sigma":
        print(f"{Crossection}\t{Name}\t{D:.2f}")
        CSVFile.write(f"{Crossection},{Name},{D},,\n")
        continue
    Ratio = F / D if D != 0 else np.nan
    print(f"{Crossection}\t{Name}\t{D:.4e}\t{F:.4e}\t{Ratio:.4f}")
    CSVFile.write(f"{Crossection},{Name},{D},{F},{Ratio}\n")
CSVFile.close()

### Comparison plot ###############
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(8, 11),
                                    gridspec_kw={'height_ratios': [3, 2, 2]})

# LET spectra
for Hist, Label, Color in ((DirectHist, "Direct spectrum sampling", DirectColor),
                           (FlatLogHist, "Flat in log biased sampling", FlatLogColor)):
    ax1.bar(Hist['lower'], Hist['value'], width=Hist['upper'] - Hist['lower'],
            align='edge', alpha=0.3, color=Color)
    ax1.errorbar(Hist['mean'], Hist['value'], Hist['error'], fmt=' ', capsize=3,
                 elinewidth=1, capthick=1, color=Color, label=Label)
ax1.set_yscale("log")
ax1.set_xscale("log")
ax1.grid(which='major')
ax1.set_ylabel("Rate per LET bin [cm-2 s-1]")
ax1.legend(loc='upper right')

TitleLines = ["CREME96 solar iron at GEO behind " + SubFolder + " Al",
              f"Total LET rate: Direct = {Direct['TotalLET']:.3e} ± {Direct['TotalLETError']:.1e}, "
              + f"FlatLog = {FlatLog['TotalLET']:.3e} ± {FlatLog['TotalLETError']:.1e} MeV cm2 mg-1 s-1"]
for X in Crossections:
    DirectSEE, FlatLogSEE = SEEResults[X['Name']]
    TitleLines.append(f"SEE rate ({X['Name']}): Direct = {DirectSEE['SEERate']:.3e} ± {DirectSEE['SEERateError']:.1e}, "
                      + f"FlatLog = {FlatLogSEE['SEERate']:.3e} ± {FlatLogSEE['SEERateError']:.1e} s-1 bit-1")
ax1.set_title("\n".join(TitleLines), fontsize=9)

# Bin by bin ratio FlatLog / Direct with propagated error
mask = (DirectHist['value'] > 0) & (FlatLogHist['value'] > 0)
Ratio = FlatLogHist['value'][mask] / DirectHist['value'][mask]
RatioError = Ratio * np.sqrt((FlatLogHist['error'][mask] / FlatLogHist['value'][mask]) ** 2
                             + (DirectHist['error'][mask] / DirectHist['value'][mask]) ** 2)
BinMean = np.sqrt(DirectHist['lower'][mask] * DirectHist['upper'][mask])  # geometric bin centre
ax2.errorbar(BinMean, Ratio, RatioError, fmt='.', capsize=3, elinewidth=1,
             capthick=1, color='C2', label="FlatLog / Direct")
ax2.axhline(1, color='grey', linestyle='--')
ax2.set_xscale("log")
ax2.grid(which='major')
ax2.set_ylabel("Ratio FlatLog / Direct")
ax2.legend(loc='upper left')

# Normalised Weibull cross sections on the ratio panel x range for orientation
ax2b = ax2.twinx()
for X, Style in zip(Crossections, ('-', '--')):
    ax2b.plot(DirectHist['lower'], f(DirectHist['lower'], X) / X['A0'], color='C3',
              linestyle=Style, alpha=0.6, label=X['Name'] + " (normalised)")
ax2b.set_ylabel("Weibull cross section / A0", color='C3')
ax2b.tick_params(axis='y', colors='C3')
ax2b.legend(loc='center left', fontsize=8)

# Entries per bin: shows where the biasing puts the statistics
for Hist, Label, Color in ((DirectHist, "Direct", DirectColor), (FlatLogHist, "FlatLog", FlatLogColor)):
    ax3.step(np.append(Hist['lower'], Hist['upper'][-1]), np.append(Hist['entries'], Hist['entries'][-1]),
             where='post', color=Color, label=Label)
ax3.set_yscale("log")
ax3.set_xscale("log")
ax3.grid(which='major')
ax3.set_xlabel("LET [MeV cm2 mg-1]")
ax3.set_ylabel("Entries per LET bin")
ax3.legend(loc='upper right')

plt.savefig(Directory + "Compare_FlatLog_vs_Direct_" + SubFolder + ".pdf", format='pdf', bbox_inches="tight")
plt.savefig(Directory + "Compare_FlatLog_vs_Direct_" + SubFolder + ".png", format='png', dpi=150, bbox_inches="tight")
plt.close('all')
print("\nSaved comparison plot and csv to " + Directory)
