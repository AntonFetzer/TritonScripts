"""Fold mono-energetic response kernels with beam spectra and fit calibration scales.

A response kernel is the dose of each tally, simulated at a grid of
mono-energetic beam energies: arrays (energy) and (energy, tally) for the dose
and its Monte Carlo error. These functions are independent of the detector and
of the campaign; what the tallies are, and which measurements they are compared
with, is up to the caller.
"""

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
