import matplotlib.pyplot as plt
from GRAS.Dependencies.TotalLETHistos import totalLETHistos
import os
from natsort import natsorted
import numpy as np
from uncertainties import ufloat

Directory = "/l/triton_work/LET_Histograms/Carrington/"
# Go through all sub folders of the directory
for SubDir in os.listdir(Directory):
    # Check if the sub folder is not a directory
    if not os.path.isdir(Directory + SubDir):
        # Skip this entry
        continue

    Path = Directory + SubDir + "/"

# Path = "/l/triton_work/LET_Histograms/Carrington/VAB-CosmicProton-mission/"

    Title = Path.split("/")[-2]


    Folders = [f for f in os.listdir(Path) if f.endswith('mm')]
    Folders = natsorted(Folders)

    ThickList = []

    for folder in Folders:
        ThickList.append(int(folder.split("m")[0]))

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
            

    Colours = ['C3', 'C1', 'C8', 'C2', 'C9', 'C0', 'C7']
    #          red, orange,yellow,green,cyan,blue, grey
    # Colours = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
    #           blue,orange,green,red,violet,brown

    # C0    blue
    # C1    Orange
    # C2    Green
    # C3    red
    # C4    purple
    # C5    brown
    # C6    pink
    # C7    grey
    # C8    yellow
    # C9    cyan

    ### LET by Entries ###############
    # plt.figure(0)
    # for i, LET in enumerate(LET_Hist_List):
    #     plt.bar(LET['lower'], LET['entries'], width=LET['upper'] - LET['lower'], align='edge', alpha=0.5, color=Colours[i], label=f"{ThickList[i]} mm")

    # plt.yscale("log")
    # plt.xscale("log")
    # plt.grid()
    # plt.title("LET Histograms " + Title)
    # plt.xlabel("LET [MeV cm2 mg-1]")
    # plt.ylabel("Number of entries per LET bin")
    # plt.legend()

    # plt.savefig(Path + "/LET-Entries.pdf", format='pdf', bbox_inches="tight")

    ### LET by Values with Error Bars ###############
    plt.figure(1)
    for i, LET in enumerate(LET_Hist_List):
        plt.bar(LET['lower'], LET['value'], width=LET['upper'] - LET['lower'], align='edge', alpha=0.3, color=Colours[i])
        plt.errorbar(LET['mean'], LET['value'], LET['error'], fmt=' ', capsize=5, elinewidth=1, capthick=1, label=f"{ThickList[i]} mm", color=Colours[i])
        # Calculate total LET by values for LETHist
        TotalLETbyValues = np.sum(LET['mean'] * LET['value'])
        print(f"Total LET by Values for {ThickList[i]} mm = {TotalLETbyValues:.2e} [MeV cm2 mg-1]")

    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title("LET Histograms " + Title)
    plt.xlabel("LET [MeV cm2 mg-1]")
    plt.ylabel("Rate per LET bin [cm-2 s-1]")
    plt.legend()

    plt.savefig(Path + "/LET-Values.pdf", format='pdf', bbox_inches="tight")
    plt.close(1)

    """ ### Eff by Entries ###############
    plt.figure(2)
    for i, Eff in enumerate(Eff_Hist_List):
        plt.bar(Eff['lower'], Eff['entries'], width=Eff['upper'] - Eff['lower'], align='edge', alpha=0.5, color=Colours[i], label=f"{ThickList[i]} mm")


    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title(f"EffLET Histograms" + Title)
    plt.xlabel("EffLET [MeV cm2 mg-1]")
    plt.ylabel("Number of entries per EffLET bin")
    plt.legend()

    plt.savefig(Path + "/EffLET-Entries.pdf", format='pdf', bbox_inches="tight")

    ### Eff by Values with Error Bars ###############
    plt.figure(3)
    for i, Eff in enumerate(Eff_Hist_List):
        plt.bar(Eff['lower'], Eff['value'], width=Eff['upper'] - Eff['lower'], align='edge', alpha=0.3, color=Colours[i])
        plt.errorbar(Eff['mean'], Eff['value'], Eff['error'], fmt=' ', capsize=5, elinewidth=1, capthick=1, label=f"{ThickList[i]} mm", color=Colours[i])
        # Calculate total EffLET by values for EffHist
        TotalEffLETbyValues = np.sum(Eff['mean'] * Eff['value'])
        print(f"Total EffLET by Values for {ThickList[i]} mm = {TotalEffLETbyValues:.2e} [MeV cm2 mg-1]")

    plt.yscale("log")
    plt.xscale("log")
    plt.grid()
    plt.title(f"EffLET Histograms" + Title)
    plt.xlabel("EffLET [MeV cm2 mg-1]")
    plt.ylabel("Rate per EffLET bin [cm-2 s-1]")
    plt.legend()

    plt.savefig(Path + "/EffLET-Values.pdf", format='pdf', bbox_inches="tight") """
