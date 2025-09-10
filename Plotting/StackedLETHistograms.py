import matplotlib.pyplot as plt
from GRAS.Dependencies.TotalLETHistos import totalLETHistos
import os
from natsort import natsorted
import numpy as np
from uncertainties import ufloat

Path = "/l/triton_work/LET_Histograms/Mono/Protons16mmAl-200micronSi/"

Title = Path.split("/")[-2]

Folders = [f for f in os.listdir(Path) if f.endswith('MeV')]  # if f.endswith('mm')
Folders = natsorted(Folders, reverse=True)  # Sort the folders in natural order

ThickList = []

for folder in Folders:
    ThickList.append(int(folder.split("MeV")[0]))

NumberOfThicknesses = len(ThickList)

LET_Hist_List = []
Eff_Hist_List = []

for folder in Folders:
    try:
        LET, Eff = totalLETHistos(Path + folder + "/Res/")
        # print("Number of Entries in LET Hist = " + f"{np.sum(LET['entries']):.3g}")

        LET_Hist_List.append(LET)
        Eff_Hist_List.append(Eff)
    except:
        print(f"Error in {folder}")
        

Colours = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15']

### LET by Entries ###############
plt.figure(0)
for i, LET in enumerate(LET_Hist_List):
    plt.bar(LET['lower'], LET['entries'], width=LET['upper'] - LET['lower'], align='edge', alpha=0.5, color=Colours[i], label=f"{ThickList[i]} MeV")

plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("LET Histograms " + Title)
plt.xlabel("LET [MeV cm2 mg-1]")
plt.ylabel("Number of entries per LET bin")
plt.legend()

plt.savefig(Path + "/LET-Entries.pdf", format='pdf', bbox_inches="tight")

### LET by Values with Error Bars ###############
plt.figure(1)
for i, LET in enumerate(LET_Hist_List):
    plt.bar(LET['lower'], LET['value'], width=LET['upper'] - LET['lower'], align='edge', alpha=0.3, color=Colours[i])
    plt.errorbar(LET['mean'], LET['value'], LET['error'], fmt=' ', capsize=5, elinewidth=1, capthick=1, label=f"{ThickList[i]} MeV", color=Colours[i])
    # Calculate total LET by values for LETHist
    TotalLETbyValues = np.sum(LET['mean'] * LET['value'])
    print(f"Total LET by Values for {ThickList[i]} MeV = {TotalLETbyValues:.2e} [MeV cm2 mg-1]")

plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("LET Histograms " + Title)
plt.xlabel("LET [MeV cm2 mg-1]")
plt.ylabel("Rate per LET bin [cm-2 s-1]")
plt.legend()

plt.savefig(Path + "/LET-Values.pdf", format='pdf', bbox_inches="tight")
plt.close(1)
