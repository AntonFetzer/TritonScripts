"""Fold mono-energetic response kernels with beam spectra and fit calibration scales.

A response kernel is the dose of each tally, simulated at a grid of
mono-energetic beam energies: arrays (energy) and (energy, tally) for the dose
and its Monte Carlo error. These functions are independent of the detector and
of the campaign; what the tallies are, and which measurements they are compared
with, is up to the caller.
"""

import csv
import os
import re

import numpy as np


def publishedSkandionSigma(energy):
    """Skandion proton energy spread [MeV] at a centre energy [MeV].

    Ardenfors et al., J. Phys.: Conf. Ser. 860 (2017) 012025: sigma_E/E falls
    linearly from 0.8 % at 60 MeV to 0.3 % at 226 MeV.
    """
    return energy * (0.8 - 0.5 * (energy - 60.0) / 166.0) / 100.0


def interpolateKernel(energies, kernel, evaluationEnergies):
    """Linearly interpolate every tally column of a kernel to new energies."""
    return np.column_stack([np.interp(evaluationEnergies, energies, kernel[:, i])
                            for i in range(kernel.shape[1])])


def foldKernel(energies, kernel, error, means, sigmas, step=0.025):
    """Fold a kernel with Gaussian spectra, truncated to the kernel's energy range.

    Args:
        energies (array): kernel energies [MeV], ascending.
        kernel, error (2D arrays): dose and Monte Carlo error, (energy, tally).
        means, sigmas (arrays): centre energy and width of each spectrum [MeV];
            a scalar sigma applies to every centre energy.
        step (float): spacing of the fine grid the kernel is interpolated to [MeV].

    Returns:
        tuple: folded dose and error, each (spectrum, tally). The interpolated
            errors are correlated within each kernel interval, so their
            quadrature sum is scaled by the number of fine points per interval;
            otherwise folding would shrink a single kernel point's error.
    """
    energies = np.asarray(energies, dtype=float)
    means = np.atleast_1d(np.asarray(means, dtype=float))
    sigmas = np.broadcast_to(np.asarray(sigmas, dtype=float), means.shape)
    fine = np.arange(energies.min(), energies.max() + step / 2, step)
    k = interpolateKernel(energies, kernel, fine)
    e = interpolateKernel(energies, error, fine)
    pointsPerInterval = np.interp(fine, energies[1:], np.diff(energies)) / step

    folded = np.empty((len(means), kernel.shape[1]))
    foldedError = np.empty_like(folded)
    for i, (mean, sigma) in enumerate(zip(means, sigmas)):
        w = np.exp(-0.5 * ((fine - mean) / sigma) ** 2)
        w /= w.sum()
        folded[i] = w @ k
        foldedError[i] = np.sqrt((w ** 2 * pointsPerInterval) @ e ** 2)
    return folded, foldedError


def fitScales(measured, sigma, model, modelError, groups, iterations=30):
    """Least-squares calibration scale per group of tallies: measured ~ scale * model.

    Each tally is weighted by 1 / (sigma^2 + (scale * modelError)^2); because the
    weight depends on the scale, the fit is iterated to convergence.

    Args:
        measured, sigma (arrays): measurement and its 1-sigma uncertainty.
        model, modelError (arrays): model prediction and its uncertainty.
        groups (dict): group name -> indices of the tallies sharing one scale.

    Returns:
        tuple: (scale by group, chi2, pull per tally, scale per tally); the pull
            is (scale * model - measured) / sigma_total.
    """
    scale = np.ones(len(measured))
    pulls = np.zeros(len(measured))
    scales = {}
    chi2 = 0.0
    for name, idx in groups.items():
        k, kerr, d, sd = model[idx], modelError[idx], measured[idx], sigma[idx]
        s = float(np.sum(d) / np.sum(k)) if np.sum(k) > 0 else 0.0
        for _ in range(iterations):
            w = 1.0 / (sd ** 2 + (s * kerr) ** 2)
            new = float(np.sum(w * d * k) / np.sum(w * k * k))
            converged = abs(new - s) <= 1e-12 * max(abs(s), 1e-300)
            s = new
            if converged:
                break
        w = 1.0 / (sd ** 2 + (s * kerr) ** 2)
        residual = s * k - d
        chi2 += float(np.sum(w * residual ** 2))
        pulls[idx] = residual * np.sqrt(w)
        scale[idx] = s
        scales[name] = s
    return scales, chi2, pulls, scale



def readKernelRuns(folder, pattern=r"^(\d+)p(\d+)MeV"):
    """Read a mono-energetic kernel from the TotalDose tables of a sweep folder.

    Every run folder below `folder` whose name matches `pattern` is one kernel
    energy, the match's two groups being the integer and fractional MeV (the
    RadEx convention 5p00MeVElectron). Its TotalDose_<folder>.csv, written by
    RadEx/RadExTID.py, holds one row per tally in the detector macro's order.

    Returns:
        tuple: energies [MeV] ascending, dose and error (energy, tally) in the
            tables' units (kRad per unit of the runs' normalisation fluence),
            and non-zero entries (energy, tally).
    """
    rows = []
    for name in sorted(os.listdir(folder)):
        match = re.match(pattern, name)
        table = os.path.join(folder, name, "TotalDose_" + name + ".csv")
        if match is None or not os.path.isfile(table):
            continue
        energy = float(match.group(1)) + float(match.group(2)) / 10 ** len(match.group(2))
        with open(table, newline="") as handle:
            reader = csv.reader(handle)
            next(reader)
            data = np.array([[float(v) for v in line[:4]] for line in reader if line])
        if np.any(data[:, 0] != np.arange(len(data))):
            raise ValueError(f"{table}: tiles are not 0..n-1 in order")
        rows.append((energy, data[:, 1], data[:, 2], data[:, 3]))
    if not rows:
        raise FileNotFoundError(f"{folder}: no run folder matching {pattern!r} with a TotalDose table")
    rows.sort(key=lambda r: r[0])
    energies = np.array([r[0] for r in rows])
    if np.any(np.diff(energies) <= 0):
        raise ValueError(f"{folder}: duplicate kernel energies")
    return (energies, np.vstack([r[1] for r in rows]), np.vstack([r[2] for r in rows]),
            np.vstack([r[3] for r in rows]))


def foldComponents(energies, kernel, error, componentEnergies, fractions, scales=1.0):
    """Fold a kernel with a discrete spectrum whose component energies are scaled.

    The kernel is interpolated linearly in energy, with zero dose at zero energy
    below its first point, and evaluated at scale * componentEnergies; every
    scale must keep the highest component inside the kernel's last energy.

    Args:
        energies (array): kernel energies [MeV], ascending.
        kernel, error (2D arrays): dose and Monte Carlo error, (energy, tally).
        componentEnergies, fractions (arrays): the spectrum's component energies
            [MeV] and their fractions (summing to one).
        scales (float or array): factor on every component energy; 1 is the
            spectrum as written.

    Returns:
        tuple: folded dose and error, each (scale, tally); the error is the
            quadrature sum of the components' interpolated errors.
    """
    energies = np.asarray(energies, dtype=float)
    scales = np.atleast_1d(np.asarray(scales, dtype=float))
    componentEnergies = np.asarray(componentEnergies, dtype=float)
    fractions = np.asarray(fractions, dtype=float)
    if scales.max() * componentEnergies.max() > energies.max() + 1e-9:
        raise ValueError("a scale takes the highest component beyond the kernel's last energy")
    e0 = np.concatenate([[0.0], energies])
    k0 = np.vstack([np.zeros(kernel.shape[1]), kernel])
    s0 = np.vstack([np.zeros(error.shape[1]), error])
    folded = np.empty((len(scales), kernel.shape[1]))
    foldedError = np.empty_like(folded)
    for i, scale in enumerate(scales):
        k = interpolateKernel(e0, k0, scale * componentEnergies)
        s = interpolateKernel(e0, s0, scale * componentEnergies)
        folded[i] = fractions @ k
        foldedError[i] = np.sqrt((fractions ** 2) @ s ** 2)
    return folded, foldedError


def scanScales(measured, sigma, models, modelErrors, groups):
    """fitScales at every row of a model grid.

    Args:
        measured, sigma (arrays): measurement and its 1-sigma uncertainty.
        models, modelErrors (2D arrays): predictions and their errors, (grid, tally).
        groups (dict): group name -> indices of the tallies sharing one scale.

    Returns:
        dict: 'chi2' (grid,), 'best' (index of the minimum), and for the best
            row 'scales' (by group), 'pulls', 'predicted' (scale * model) and
            'dof' (tallies minus groups). A row with a non-positive model sum
            in some group gets chi2 infinity.
    """
    chi2 = np.full(len(models), np.inf)
    for i in range(len(models)):
        if all(np.sum(models[i][idx]) > 0 for idx in groups.values()):
            chi2[i] = fitScales(measured, sigma, models[i], modelErrors[i], groups)[1]
    best = int(np.argmin(chi2))
    scales, _, pulls, scale = fitScales(measured, sigma, models[best], modelErrors[best], groups)
    return {"chi2": chi2, "best": best, "scales": scales, "pulls": pulls,
            "predicted": scale * models[best], "dof": len(measured) - len(groups)}


def chi2Interval(grid, chi2, level=1.0):
    """Best grid value and the interval where chi2 stays within `level` of its minimum.

    Returns:
        tuple: (best, low, high, atEdge); atEdge says the interval touches an
            end of the grid, so it is a bound, not an interval.
    """
    grid = np.asarray(grid, dtype=float)
    chi2 = np.asarray(chi2, dtype=float)
    best = int(np.argmin(chi2))
    inside = np.flatnonzero(chi2 <= chi2[best] + level)
    low, high = int(inside.min()), int(inside.max())
    return float(grid[best]), float(grid[low]), float(grid[high]), bool(low == 0 or high == len(grid) - 1)
