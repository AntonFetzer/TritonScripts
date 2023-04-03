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

Path = "/home/anton/Desktop/triton_work/3MatTriangles/PE-Al-Pb/Res/"

Materials = ["Polyethylene", "Aluminium", "Lead"]
Mat = ["PW", "Al", "Pb"]

plotTernary(Path, "Proton", Materials, Mat)
#plotTernary(Path + F + "/Res/", "Electron", Materials, Mat)
#plotTernary(Path + F + "/Res/", "Total", Materials, Mat)
