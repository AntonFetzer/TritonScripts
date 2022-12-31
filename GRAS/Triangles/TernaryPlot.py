import matplotlib.pyplot as plt
import mpltern  # noqa: F401
import numpy as np
import os
from GRAS.Dependencies.TotalKRadGras import totalkRadGras
from matplotlib import cm
import matplotlib as mpl
from uncertainties import ufloat
from csv import writer


### ------------------------------------ Do not touch !!!! It works but I don't know how and why !!!!!!

def plotTernary(path, particle, materials, mat):

    Fig = 0
    Data = []

    if "Prot" in particle:
        Data = totalkRadGras(path, particle)
        Fig = 0
    elif "Elec" in particle:
        Data = totalkRadGras(path, particle)
        Fig = 1
    if "Total" in particle:
        Fig = 2
        Prot = totalkRadGras(path, "Prot")
        Elec = totalkRadGras(path, "Elec")
        Data = Prot + Elec
        Data[1] = np.sqrt(Prot[1]**2 + Elec[1]**2)

    Data = Data / 2 ############## This is needed becasue GRAS wrongly normalises the surface area of the particle source

    N = 30
    Offset = int((N + 1) * N / 2)

    print("The average Dose is", np.mean(Data[0]))
    print("The average Dose of the main triangle is", np.mean(Data[0][:Offset]))
    print("The average Dose of the upside down triangle is", np.mean(Data[0][Offset:]))
    print("The difference in Dose is", 100 * (np.mean(Data[0][Offset:]) - np.mean(Data[0][:Offset])) / np.mean(Data[0]), "%")

    ColorData = Data[0]
    Min = np.min(ColorData)
    ColorData = ColorData - Min
    Max = np.max(ColorData)
    ColorData = ColorData / Max

    Max = np.max(Data[0])
    Min = np.min(Data[0])
    MinID = np.argmin(Data[0])
    MinErr = Data[1][MinID]

    Colors = cm.turbo(ColorData)
    print(len(ColorData))

    fig = plt.figure(Fig)
    ax = fig.add_subplot(projection='ternary', ternary_scale=100)

    for x in range(N):
        for y in range(N - x):
            n = N - 1
            r = x / n
            g = y / n
            b = 1 - r - g

            if b < 0:
                print(x, y, r, g, b)
                b = 0

            ID = int(x * N + y - (x * (x - 1) / 2))

            ax.fill([r + 1 / n, r, r], [b, b + 1 / n, b], [g, g, g + 1 / n], color=Colors[ID], linewidth=0)

            if ID == MinID:
                ax.fill([r + 1 / n, r, r], [b, b + 1 / n, b], [g, g, g + 1 / n], color="white", linewidth=0)
                rMin = r
                bMin = b
                gMin = g


    for x in range(N - 1):
        for y in range(N - x - 1):
            n = N - 1
            r = x / n + 1 / (3 * n)
            g = y / n + 1 / (3 * n)
            b = 1 - r - g

            h = 1 / (3 * n)  # Don't ask why
            H = 2 / (3 * n)  # This works. I don't know why !

            if b < 0:
                print(x, y, r, g, b)
                b = 0

            ID = int(Offset + x * (N - 1) + y - (x * (x - 1) / 2))

            ##          top axis                Left                    Right
            ax.fill([r - h, r + H, r + H], [b + H, b - h, b + H], [g + H, g + H, g - h], color=Colors[ID], linewidth=0)

            if ID == MinID:
                ax.fill([r - h, r + H, r + H], [b + H, b - h, b + H], [g + H, g + H, g - h], color="white", linewidth=0)
                rMin = r
                bMin = b
                gMin = g


    ax.set_tlabel("Layer 1 [%]: " + Materials[0] + " [%]")
    ax.set_rlabel("Layer 2 [%]: " + Materials[1] + " [%]")
    ax.set_llabel("Layer 3 [%]: " + Materials[2] + " [%]")
    ax.taxis.set_label_position('tick1')
    ax.laxis.set_label_position('tick1')
    ax.raxis.set_label_position('tick1')
    ax.tick_params(labelrotation='horizontal')
    ax.set_title(particle + " ionising dose behind 1.5 g/cm2 of three layer shielding", pad=20) # \n of " + materials[0] + " on " + materials[1] + " on " + materials[2])

    cax = ax.inset_axes([1.02, 0.1, 0.05, 0.9], transform=ax.transAxes)
    norm = mpl.colors.Normalize(vmin=Min, vmax=Max)
    #print(norm)
    Map = mpl.cm.ScalarMappable(norm=norm, cmap=cm.turbo)
    colorbar = fig.colorbar(Map, cax=cax, label='Some Units')
    colorbar.set_label('Ionizing Dose per Month [krad]', rotation=270, va='baseline')

    fig.text(0.04, 0.67, "Minimum dose per month:\n" + str(ufloat(Min, MinErr)) + " kRad.\nAt composition:\n"
             + str(round(rMin*100)) + "% " + materials[0] + " on \n"
             + str(round(gMin*100)) + "% " + materials[1] + " on \n"
             + str(round(bMin*100)) + "% " + materials[2], fontsize='large')

    #plt.show()
    plt.savefig(path + "../Plot/" + particle + "TernaryPlot.eps", format='eps')
    plt.close(Fig)

    if "Total" in particle:
        CSVFile = open("/home/anton/Desktop/triton_work/3MatTriangles/3MatSum.csv", 'a')
        List = (str(round(rMin*100)), "\% " + mat[0] + " &",
                str(round(gMin * 100)), "\% " + mat[1] + " &",
                str(round(bMin * 100)), "\% " + mat[2] + " & \\num{", ufloat(Min, MinErr), "} \\\\")

        String = ', '.join(map(str, List))
        print(String)
        CSVFile.writelines(String + "\n")
        CSVFile.close()



if __name__ == "__main__":

    Path = "/home/anton/Desktop/triton_work/3MatTriangles/ToProcess/"

    Folders = [f for f in os.listdir(Path) if "-" in f]

    print(Folders)

    for F in Folders:

        if not F:
            break

        Mat = F.split("-")

        print(Mat)
        Materials = ["", "", ""]

        for i in range(3):
            if "Al" in Mat[i]:
                Materials[i] = "Aluminium"
            elif "PE" in Mat[i]:
                Materials[i] = "Polyethylene"
            elif "W" in Mat[i]:
                Materials[i] = "Tungsten"
            elif "Pb" in Mat[i]:
                Materials[i] = "Lead"
            elif "FR4" in Mat[i]:
                Materials[i] = "FR4"

        print(Materials)
        #plotTernary(Path, "Proton", Materials, Mat)
        #plotTernary(Path, "Electron", Materials, Mat)
        plotTernary(Path + F + "/Res/", "Total", Materials, Mat)
