import matplotlib.pyplot as plt
import mpltern  # noqa: F401
import numpy as np
from GRAS.Dependencies.TotalKRadGras import totalkRadGras
from matplotlib import cm
import matplotlib as mpl
from uncertainties import ufloat
from GRAS.Triangles.CombinedColormaps import create_average_colormap
import sigfig


### ------------------------------------ Do not touch !!!! It works but I don't know how and why !!!!!!

def plotTernary(path, particle, materials, mat):
    Fig = 0
    Data = []

    if "Elec" in particle:
        Data = totalkRadGras(path, particle)
        Fig = 0
        ColMap = cm.viridis
    elif "Prot" in particle:
        Data = totalkRadGras(path, particle)
        Fig = 1
        ColMap = cm.plasma
    if "Total" in particle:
        Fig = 2
        Prot = totalkRadGras(path, "Prot")
        Elec = totalkRadGras(path, "Elec")
        Data = Prot + Elec
        Data[1] = np.sqrt(Prot[1] ** 2 + Elec[1] ** 2)
        ColMap = create_average_colormap(cm.viridis, cm.plasma)
        Prot = Prot / 2
        Elec = Elec / 2

    Data = Data / 2  ############## This is needed becasue GRAS wrongly normalises the surface area of the particle source

    N = 30
    Offset = int((N + 1) * N / 2)

    print("The average Dose is", np.mean(Data[0]))
    print("The average Dose of the main triangle is", np.mean(Data[0][:Offset]))
    print("The average Dose of the upside down triangle is", np.mean(Data[0][Offset:]))
    print("The difference in Dose is", 100 * (np.mean(Data[0][Offset:]) - np.mean(Data[0][:Offset])) / np.mean(Data[0]),
          "%")

    ColorData = Data[0]
    Min = np.min(ColorData)
    ColorData = ColorData - (Min*0.999)
    Max = np.max(ColorData)
    ColorData = ColorData / Max
    ColorData = np.log(ColorData)
    Min = np.min(ColorData)
    ColorData = ColorData - Min
    Max = np.max(ColorData)
    ColorData = ColorData / Max

    # Define the function to reverse the transformation
    def reverse_transformation(norm_value, Min, Max):
        # Reverse the steps in the normalization
        value = norm_value * Max  # Reverse division by Max
        value = value + Min  # Reverse subtraction of Min
        value = np.exp(value)  # Reverse log
        value = value * Max  # Reverse division by Max
        value = value + (Min * 0.999)  # Reverse subtraction of Min*0.999
        return value

    Max = np.max(Data[0])
    Min = np.min(Data[0])
    MinID = np.argmin(Data[0])
    MinErr = Data[1][MinID]

    Colors = ColMap(ColorData)
    print(len(ColorData))

    fig = plt.figure(Fig)
    ax = fig.add_subplot(projection='ternary', ternary_sum=100)

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

    ax.set_tlabel("Layer 1 [%]: " + materials[0] + " [%]")
    ax.set_rlabel("Layer 2 [%]: " + materials[1] + " [%]")
    ax.set_llabel("Layer 3 [%]: " + materials[2] + " [%]")
    ax.taxis.set_label_position('tick1')
    ax.laxis.set_label_position('tick1')
    ax.raxis.set_label_position('tick1')
    ax.tick_params(labelrotation='horizontal')
    #ax.set_title(particle + " ionising dose behind 1.5 g/cm2 of three layer shielding",
    #             pad=20)  # \n of " + materials[0] + " on " + materials[1] + " on " + materials[2])

    cbar_ticks = np.geomspace(Min, Max, num=8)  # Adjust the number of ticks as needed
    cbar_ticks = [sigfig.round(num, sigfigs=3) for num in cbar_ticks]

    cax = ax.inset_axes([1, 0.1, 0.05, 0.95], transform=ax.transAxes)
    norm = mpl.colors.Normalize(vmin=np.log(Min), vmax=np.log(Max))
    # print(norm)
    Map = mpl.cm.ScalarMappable(norm=norm, cmap=ColMap)
    colorbar = fig.colorbar(Map, cax=cax, ticks=np.log(cbar_ticks))
    colorbar.set_ticklabels(cbar_ticks)
    colorbar.set_label('Ionizing Dose per Month [krad]', rotation=270, va='baseline')

    fig.text(0.1, 0.70, "Minimum dose per month:\n" + str(ufloat(Min, MinErr)) + " kRad.\nAt composition:\n"
             + str(round(rMin * 100)) + "% " + materials[0] + " on \n"
             + str(round(gMin * 100)) + "% " + materials[1] + " on \n"
             + str(round(bMin * 100)) + "% " + materials[2], fontsize='large')

    plt.show()
    # plt.savefig(path + "../Plot/" + particle + "TernaryPlot.pdf", format='pdf')
    # plt.close(Fig)
    '''
    if "Total" in particle:
        CSVFile = open("/l/triton_work/3MatTriangles/3MatSum.csv", 'a')
        List = (str(round(rMin * 100)), "\% " + mat[0] + " &",
                str(round(gMin * 100)), "\% " + mat[1] + " &",
                str(round(bMin * 100)), "\% " + mat[2] + " & \\num{", ufloat(Min, MinErr), "} \\\\")

        String = ', '.join(map(str, List))
        print(String)
        CSVFile.writelines(String + "\n")
        CSVFile.close()
    '''

    print(ufloat(Min, MinErr))

    if "Total" in particle:
        # Extract minimum and its error for Total Dose
        total_dose_min = Min
        total_dose_min_err = MinErr

        # Extract minimum and its error for Electron Dose
        elec_dose_min = Elec[0][MinID]
        elec_dose_min_err = Elec[1][MinID]

        # Extract minimum and its error for Proton Dose
        proton_dose_min = Prot[0][MinID]
        proton_dose_min_err = Prot[1][MinID]

        CSVFilePath = path + "../" + "MinDoses.csv"

        # Convert to string representation with uncertainties
        total_dose_str = str(ufloat(total_dose_min, total_dose_min_err)).split('+/-')
        elec_dose_str = str(ufloat(elec_dose_min, elec_dose_min_err)).split('+/-')
        proton_dose_str = str(ufloat(proton_dose_min, proton_dose_min_err)).split('+/-')

        # Construct the line to write to the CSV
        line = (f"{materials[0]},"
                f"{materials[1]},"
                f"{materials[2]},"
                f"{round(rMin * 100)},"
                f"{round(gMin * 100)},"
                f"{100 - round(rMin * 100) - round(gMin * 100)},"
                f"{elec_dose_str[0]},{elec_dose_str[1]},"
                f"{proton_dose_str[0]},{proton_dose_str[1]},"
                f"{total_dose_str[0]},{total_dose_str[1]}\n")

        # Check if the file exists and write header if it doesn't
        try:
            with open(CSVFilePath, 'r') as file:
                pass
        except FileNotFoundError:
            with open(CSVFilePath, 'w') as file:
                file.write(
                    "Material A, Material B, Material C, % A, % B, % C, Electron Dose, Electron Error, Proton Dose, Proton Error, Total Dose, Total Error\n")

        # Write data to the CSV file
        with open(CSVFilePath, 'a') as file:
            file.write(line)


if __name__ == "__main__":
    Path = "/l/triton_work/3MatTriangles/Al-Ti-Ta/Res/"

    Materials = ["Aluminium", "Titanium", "Tantalum"]
    Mat = ["Al", "Ti", "Ta"]

    # plotTernary(Path, "Electron", Materials, Mat)
    # plotTernary(Path, "Proton", Materials, Mat)
    plotTernary(Path, "Total", Materials, Mat)
