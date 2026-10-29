from Dependencies.TotalDose import totalDose
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, MultipleLocator
from os import path
from Dependencies.MergeTotalDose import mergeTotalDose
from Read.ReadSD2Q import readSDQ2
import csv
import sys


Path = "/l/triton_work/Shielding_Curves/Mono/"
res_suffix = "/Res/"

Labels = [
    "32MeVProtons",
    "64MeVProtons",
    "128MeVProtons"
]

YLabel = "Ionising Dose [krad/1E11 protons]"
ReplicaError = False  # True for biased GPS sources: error from the replica scatter
MinReplicas = 5
ShieldoseFile = None  # SHIELDOSE-2Q spenvis_sqo.txt: compare <electron label> <proton label> with it
Months = None  # Divide all doses by this number of months, e.g. 132 for an 11-year fluence
Title = None  # Plot title, e.g. --title=TEXT; the default depends on the plot
ShieldoseName = "SHIELDOSE-2Q"  # Legend name, e.g. --shieldose-name="SHIELDOSE-2Q AE9/AP9" when the models differ
XLim = None  # (min, max) in mm, e.g. --xlim=0,6
YLim = None  # (min, max) in krad, e.g. --ylim=1e-2,1e3
Output = None  # PDF name in <path>, e.g. --output=ShieldingCurves_Physics.pdf; default depends on the plot
FigSize = None  # (width, height) in inches, e.g. --figsize=4,6
Names = None  # Legend names, one per label, separated by ';', e.g. --names="EM3;EM4"; default: the labels
LegendLoc = "best"  # Matplotlib legend location, e.g. --legend-loc="upper right"
ErrorbarStyles = {
    'default': {'capsize': 5, 'linestyle': ''},
    'compact': {'capsize': 1.5, 'elinewidth': 1.5, 'capthick': 2, 'linestyle': ''},  # --errorbar-style=compact
}
ErrorbarStyle = ErrorbarStyles['default']

# Optional: python3 -m Plotting.CompareShieldingCurves <path> <label> ... [--replica-error] [--min-replicas=N]
#   [--ylabel=TEXT] [--title=TEXT] [--output=NAME.pdf] [--figsize=W,H] [--names="A;B"] [--errorbar-style=compact]
#   [--legend-loc="upper right"]
#   [--xlim=MIN,MAX --ylim=MIN,MAX] [--shieldose=spenvis_sqo.txt --shieldose-name=TEXT --months=N]
Args = [a for a in sys.argv[1:] if not a.startswith("--")]
if Args:
    Path = Args[0]
    Labels = Args[1:] or Labels
ReplicaError = ReplicaError or "--replica-error" in sys.argv[1:]
for a in sys.argv[1:]:
    if a.startswith("--ylabel="):
        YLabel = a.split("=", 1)[1]
    elif a.startswith("--min-replicas="):
        MinReplicas = int(a.split("=", 1)[1])
    elif a.startswith("--shieldose="):
        ShieldoseFile = a.split("=", 1)[1]
    elif a.startswith("--shieldose-name="):
        ShieldoseName = a.split("=", 1)[1]
    elif a.startswith("--months="):
        Months = float(a.split("=", 1)[1])
    elif a.startswith("--title="):
        Title = a.split("=", 1)[1]
    elif a.startswith("--xlim="):
        XLim = tuple(float(v) for v in a.split("=", 1)[1].split(","))
    elif a.startswith("--ylim="):
        YLim = tuple(float(v) for v in a.split("=", 1)[1].split(","))
    elif a.startswith("--output="):
        Output = a.split("=", 1)[1]
    elif a.startswith("--figsize="):
        FigSize = tuple(float(v) for v in a.split("=", 1)[1].split(","))
    elif a.startswith("--legend-loc="):
        LegendLoc = a.split("=", 1)[1]
    elif a.startswith("--names="):
        Names = a.split("=", 1)[1].split(";")
    elif a.startswith("--errorbar-style="):
        ErrorbarStyle = ErrorbarStyles[a.split("=", 1)[1]]

if Names is not None and len(Names) != len(Labels):
    sys.exit("ERROR !!! --names needs one name per label")

Colours = ['C0', 'C1', 'C2', 'C8', 'C3', 'C9', 'C7', 'k', 'C4', 'C5', 'C6']

# Initialise the Shielding Curves dictionary to store the shielding curves
SC = {}
for label in Labels:
    SC[label] = totalDose(path.join(Path, label) + res_suffix, replica_error=ReplicaError, min_replicas=MinReplicas)


NumTiles = np.shape(SC[Labels[0]]['dose'])[0]

x = np.linspace(0, 10, num=NumTiles, endpoint=True)

if ShieldoseFile is None:
    ## Plot shielding curves with error bars
    plt.figure(1, FigSize)

    # print("Shielding curves: ", SC)

    for i, label in enumerate(Labels):
       plt.errorbar(x, SC[label]['dose'], SC[label]['error'], fmt='', label=Names[i] if Names else label,
                    color=Colours[i], **ErrorbarStyle)

    ####### Plot 10kRad line #########
    CriticalDose = [1 for i in x]
    plt.plot(x, CriticalDose, color='k', linewidth=2, label='1 krad')
    CriticalDose = [10 for i in x]
    plt.plot(x, CriticalDose, '--', color='k', linewidth=2, label='10 krad')

    plt.grid()
    plt.yscale("log")
    # plt.xscale("log")
    plt.title(Title or "Ionising Dose behind shielding of varying thickness")
    plt.xlabel("Aluminium Thickness [mm]")
    plt.ylabel(YLabel)
    plt.legend(loc=LegendLoc)
    if XLim:
        plt.xlim(XLim)
    if YLim:
        plt.ylim(YLim)

    plt.savefig(path.join(Path, Output or "ShieldingCurves.pdf"), format='pdf', bbox_inches="tight")

else:
    ## Compare GRAS electron, proton and total dose with SHIELDOSE-2Q
    if len(Labels) != 2:
        sys.exit("ERROR !!! --shieldose needs exactly two labels: <electron label> <proton label>")
    Scale = 1 / Months if Months else 1
    PerUnit = " per month" if Months else ""

    Electrons = SC[Labels[0]]
    Protons = SC[Labels[1]]
    Total = mergeTotalDose([Electrons, Protons])

    plt.figure(2)
    # SHIELDOSE-2Q: electron dose includes bremsstrahlung. The total is the trapped
    # part only (electrons + bremsstrahlung + trapped protons), as simulated in GRAS.
    # Both GRAS tiles and SHIELDOSE-2Q thicknesses are in mm of aluminium.
    SD = readSDQ2(ShieldoseFile)
    SDx = SD['Thickness']
    SDElectrons = SD['Electrons'] + SD['Bremsstrahlung']
    plt.plot(SDx, SDElectrons * Scale, 'C0--', label=ShieldoseName + " Electron Dose")
    plt.plot(SDx, SD['Trapped Protons'] * Scale, 'C1-', label=ShieldoseName + " Proton Dose")
    plt.plot(SDx, (SDElectrons + SD['Trapped Protons']) * Scale, 'C2:', label=ShieldoseName + " Total Dose")

    plt.plot(x, [1] * NumTiles, color='k', linewidth=2, label='1 krad' + PerUnit)

    plt.errorbar(x, Electrons['dose'] * Scale, Electrons['error'] * Scale, fmt='C0 ', capsize=4,
                 label="GRAS Electron Dose", elinewidth=1, capthick=1)
    plt.errorbar(x, Protons['dose'] * Scale, Protons['error'] * Scale, fmt='C1 ', capsize=4,
                 label="GRAS Proton Dose", elinewidth=1, capthick=1)
    plt.errorbar(x, Total['dose'] * Scale, Total['error'] * Scale, fmt='C2 ', capsize=4,
                 label="GRAS Total Ionizing Dose", elinewidth=1, capthick=1)

    plt.xlim(XLim if XLim else (-0.2, x[-1] + 0.2))
    plt.yscale("log")
    if YLim:
        plt.ylim(YLim)

    ax = plt.gca()
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(LogLocator(subs='all', numticks=10))

    plt.grid(which='both', linestyle='-', linewidth=0.5)
    plt.title(Title or "Ionising dose from trapped particles (AE9/AP9) on GTO", fontsize=14)
    plt.xlabel("Aluminium Thickness [mm]", fontsize=12)
    plt.ylabel("Ionizing dose in silicon [krad" + PerUnit + "]", fontsize=12)
    plt.legend(fontsize=12)

    plt.savefig(path.join(Path, Output or "ShieldoseComparison_" + Labels[0] + "-" + Labels[1] + ".pdf"),
                format='pdf', bbox_inches="tight")
