"""Per-tally diagnostic plots of a pooled GRAS dose result."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plotDose(dose, error, title, path):
    """Dose per tally with error bars, logarithmic axis, line at 1 kRad."""
    plt.figure()
    plt.errorbar(np.arange(len(dose)), dose, yerr=error, fmt=' ', capsize=5, elinewidth=1, capthick=1, label='Dose')
    plt.axhline(y=1, color='r', linestyle='--', label='1 kRad')
    plt.title(title, fontsize=10)
    plt.xlabel('Tile number')
    plt.ylabel('Dose [kRad]')
    plt.yscale('log')
    plt.grid(which='both')
    plt.legend()
    plt.savefig(path, format='pdf', bbox_inches="tight")
    plt.close()


def plotRelativeError(dose, error, path):
    """Relative error per tally in percent, on a linear axis starting at zero."""
    relative = np.divide(100 * error, dose, out=np.zeros(len(dose)), where=dose != 0)
    plt.figure()
    plt.plot(relative, '.', label='Relative Error')
    plt.axhline(y=1, color='r', linestyle='--', label='1% error')
    plt.title('Relative Error in %')
    plt.xlabel('Tile number')
    plt.ylabel('Relative Error [%]')
    plt.ylim(bottom=0)
    plt.grid(which='both')
    plt.legend()
    plt.savefig(path, format='pdf', bbox_inches="tight")
    plt.close()


def plotNonZeros(nonZeros, path):
    """Number of non-zero entries per tally, logarithmic axis."""
    plt.figure()
    plt.plot(nonZeros, '.', label='Non-zero entries')
    plt.axhline(y=1, color='r', linestyle='--', label='1 entry')
    plt.title('Number of non-zero entries')
    plt.xlabel('Tile number')
    plt.ylabel('Number of non-zero entries')
    plt.yscale('log')
    plt.grid(axis='x', which='both')
    plt.grid(axis='y', which='major')
    plt.legend()
    plt.savefig(path, format='pdf', bbox_inches="tight")
    plt.close()
