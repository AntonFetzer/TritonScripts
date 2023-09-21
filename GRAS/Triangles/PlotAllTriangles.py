import matplotlib.pyplot as plt
import mpltern  # noqa: F401
import numpy as np
import os
from GRAS.Dependencies.TotalKRadGras import totalkRadGras
from matplotlib import cm
import matplotlib as mpl
from uncertainties import ufloat
from csv import writer
from GRAS.Triangles.TernaryPlot import plotTernary

Path = "/l/triton_work/3MatTriangles/DONE/"

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
    plotTernary(Path + F + "/Res/", "Proton", Materials, Mat)
    plotTernary(Path + F + "/Res/", "Electron", Materials, Mat)
    plotTernary(Path + F + "/Res/", "Total", Materials, Mat)
