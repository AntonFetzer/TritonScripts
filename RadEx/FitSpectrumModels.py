"""Spectral models with more freedom than one energy scale, on a campaign's measured doses.

    python3 -m RadEx.FitSpectrumModels --kernel-dir <sweep> --campaign <json> --output-dir <dir>
        [--layer oxide|die] [--phase late]

Models: the production spectrum scaled (alpha); scaled and stretched about its mean
(alpha, width); scaled with a free fluence fraction below --low-mev (alpha, fraction);
Gaussian (mu, sigma); a free non-negative spectrum on the kernel energies (NNLS,
absolute, no separate scale); the same with a second-difference smoothing penalty.
Same doses, kernel and error model as RadEx.FitScaledSpectrum, on the campaign's
measurement phase by default. Writes spectrum-models-<layer>.json and .pdf (the fitted
spectra as fluence densities and the pulls) into --output-dir.
"""
import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from Dependencies.EnergyKernel import (chi2Interval, fitSpectrumWeights, foldComponents, foldKernel,  # noqa: E402
                                       readKernelRuns, scaleSpectrum, scanScales)
from RadEx import KernelFit as kf  # noqa: E402
from Read.ReadGPSSources import readGPSSources  # noqa: E402

STEP = 0.025
WIDTHS, LOW_FRACTIONS = np.arange(0.4, 1.61, 0.05), np.arange(0.0, 0.501, 0.025)
SIGMAS = np.arange(0.05, 4.0001, 0.05)

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--kernel-dir", type=Path, required=True)
parser.add_argument("--campaign", type=Path, required=True)
parser.add_argument("--output-dir", type=Path, required=True)
parser.add_argument("--layer", help="detector layer of the campaign (default: its default_layer)")
parser.add_argument("--phase", help="measured phase (default: the campaign's measurement_phase)")
parser.add_argument("--low-mev", type=float, default=3.0, help="edge of the free low-energy fraction [MeV]")
parser.add_argument("--mu-min", type=float, default=3.0, help="lowest Gaussian centre [MeV]")
parser.add_argument("--alpha-min", type=float, default=0.80, help="lowest spectrum scale")
args = parser.parse_args()
campaign = kf.load_campaign(args.campaign)
layer = args.layer or campaign["default_layer"]
phase = args.phase or campaign["measurement_phase"]
LOW_MEV, SYSTEMATIC = args.low_mev, campaign["systematic"]
outDir = args.output_dir
outDir.mkdir(parents=True, exist_ok=True)
energies, dose, error, _ = readKernelRuns(args.kernel_dir)
kernel, kerr = dose[:, kf.layer_slice(campaign, layer)], error[:, kf.layer_slice(campaign, layer)]
gps = readGPSSources(str(campaign["spectrum"]))
E, F = gps["energy"], gps["fraction"]
meanEnergy = float(F @ E)
ALPHAS = np.arange(args.alpha_min, energies.max() / E.max() + 1e-9, 0.005)
MUS = np.arange(args.mu_min, energies.max() + STEP / 2, STEP)
meas = kf.measurement(campaign, phase)
measured, sigma, fluence = meas["value"], meas["sigma"], campaign["phases"][phase]["log_fluence"]
nChannels = len(campaign["channels"])
groupsOf = {model: kf.scale_groups(campaign, campaign["channels"], model) for model in ("shared", "per-type")}
spacing = float(np.median(np.diff(energies)))

def density(componentEnergies, fractions, edges):
    """Fluence fraction per MeV of a discrete spectrum, histogrammed on `edges`."""
    h, _ = np.histogram(componentEnergies, bins=edges, weights=fractions)
    return h / np.diff(edges)

edges = np.arange(0.0, energies.max() + 0.251, 0.25)
results, spectra = {}, {}
for model, groups in groupsOf.items():
    out = {}
    # 1. PDD scaled
    folded = foldComponents(energies, kernel, kerr, E, F, ALPHAS)
    scan = scanScales(measured, sigma, *folded, groups)
    a, lo, hi, edge = chi2Interval(ALPHAS, scan["chi2"])
    out["pdd_scaled"] = {"params": {"alpha": a}, "interval": {"alpha": [lo, hi]}, "chi2": float(scan["chi2"][scan["best"]]),
                         "dof": scan["dof"] - 1, "mean_energy_MeV": a * meanEnergy, "pulls": scan["pulls"].tolist()}
    spectra[(model, "pdd_scaled")] = (scaleSpectrum(E, a, 1.0, meanEnergy), F)
    # 2. PDD scaled and stretched
    best = None
    for w in WIDTHS:
        Es = scaleSpectrum(E, 1.0, w, meanEnergy)
        ok = ALPHAS[ALPHAS * Es.max() <= energies.max() + 1e-9]
        s = scanScales(measured, sigma, *foldComponents(energies, kernel, kerr, Es, F, ok), groups)
        if best is None or s["chi2"][s["best"]] < best[0]:
            best = (float(s["chi2"][s["best"]]), float(ok[s["best"]]), float(w), s)
    Es = scaleSpectrum(E, best[1], best[2], meanEnergy)
    out["pdd_scaled_stretched"] = {"params": {"alpha": best[1], "width": best[2]}, "chi2": best[0], "dof": best[3]["dof"] - 2,
                                   "mean_energy_MeV": float(F @ Es), "pulls": best[3]["pulls"].tolist()}
    spectra[(model, "pdd_scaled_stretched")] = (Es, F)
    # 3. PDD scaled with a free fraction below LOW_MEV
    low = E < LOW_MEV
    best = None
    for f in LOW_FRACTIONS:
        Fl = F.copy()
        Fl[low] = F[low] / F[low].sum() * f if F[low].sum() > 0 else 0.0
        Fl[~low] = F[~low] / F[~low].sum() * (1 - f)
        s = scanScales(measured, sigma, *foldComponents(energies, kernel, kerr, E, Fl, ALPHAS), groups)
        if best is None or s["chi2"][s["best"]] < best[0]:
            best = (float(s["chi2"][s["best"]]), float(ALPHAS[s["best"]]), float(f), s, Fl)
    out["pdd_scaled_lowfraction"] = {"params": {"alpha": best[1], f"fraction_below_{LOW_MEV:g}_MeV": best[2]}, "chi2": best[0],
                                     "dof": best[3]["dof"] - 2, "mean_energy_MeV": float(best[4] @ (best[1] * E)),
                                     "pdd_fraction_below": float(F[low].sum()), "pulls": best[3]["pulls"].tolist()}
    spectra[(model, "pdd_scaled_lowfraction")] = (best[1] * E, best[4])
    # 4. Gaussian
    grid = np.array([scanScales(measured, sigma, *foldKernel(energies, kernel, kerr, MUS, s, STEP), groups)["chi2"] for s in SIGMAS])
    j, i = np.unravel_index(int(np.nanargmin(grid)), grid.shape)
    s = scanScales(measured, sigma, *foldKernel(energies, kernel, kerr, np.array([MUS[i]]), SIGMAS[j], STEP), groups)
    out["gaussian"] = {"params": {"mu": float(MUS[i]), "sigma": float(SIGMAS[j])}, "chi2": float(grid[j, i]), "dof": s["dof"] - 2,
                       "mean_energy_MeV": float(MUS[i]), "pulls": s["pulls"].tolist()}
    fine = np.arange(energies.min(), energies.max() + STEP / 2, STEP)
    g = np.exp(-0.5 * ((fine - MUS[i]) / SIGMAS[j]) ** 2); g /= g.sum()
    spectra[(model, "gaussian")] = (fine, g)
    results[model] = out

# 5 and 6. free non-negative spectrum on the kernel energies, absolute (shared scale implicit)
for name, smoothing in (("free_nnls", 0.0), ("free_nnls_smoothed", None)):
    if smoothing is None:   # the smallest penalty that costs one unit of chi2 against the free fit
        base = results["shared"]["free_nnls"]["chi2"]
        smoothing = next((lam for lam in np.logspace(-6, 6, 97)
                          if fitSpectrumWeights(measured, sigma, kernel, kerr, lam)[3] >= base + 1.0), 1e6)
    w, pred, pulls, chi2 = fitSpectrumWeights(measured, sigma, kernel, kerr, smoothing)
    nonzero = int(np.sum(w > 0))
    results["shared"][name] = {"params": {"weights_electrons_cm2": w.tolist(), "smoothing": float(smoothing)}, "chi2": chi2,
                               "dof": nChannels - nonzero, "free_weights": nonzero, "fluence_electrons_cm2": float(w.sum()),
                               "calibration_ratio": float(w.sum() / fluence), "mean_energy_MeV": float(energies @ w / w.sum()),
                               "pulls": pulls.tolist()}
    spectra[("shared", name)] = (energies, w / w.sum())

with open(outDir / f"spectrum-models-{layer}.json", "w") as handle:
    json.dump({"layer": layer, "phase": phase, "systematic": SYSTEMATIC, "production_mean_energy_MeV": meanEnergy,
               "kernel_energies_MeV": energies.tolist(), **results}, handle, indent=2)

print(f"{'scale':9s} {'model':24s} {'chi2':>6s} {'dof':>4s}  mean E   parameters")
for model, out in results.items():
    for name, r in out.items():
        print(f"{model:9s} {name:24s} {r['chi2']:6.1f} {r['dof']:4d}  {r['mean_energy_MeV']:5.2f}   "
              + ", ".join(f"{k} {v:.3g}" if isinstance(v, float) else f"{k} {v}" for k, v in r["params"].items() if k != "weights_electrons_cm2"))
        if "weights_electrons_cm2" in r["params"]:
            w = np.array(r["params"]["weights_electrons_cm2"]); nz = w > 0
            print("            fractions:", ", ".join(f"{e:.2f} MeV {f:.3f}" for e, f in zip(energies[nz], w[nz] / w.sum())))

fig, axes = plt.subplots(2, 2, figsize=(17, 11))
styles = {"pdd_scaled": ("C0", "-"), "pdd_scaled_stretched": ("C1", "--"), "pdd_scaled_lowfraction": ("C2", "-."),
          "gaussian": ("C4", ":"), "free_nnls": ("C3", "-"), "free_nnls_smoothed": ("C5", "--")}
for col, model in enumerate(groupsOf):
    ax = axes[0, col]
    ax.stairs(density(E, F, edges), edges, color="black", lw=2, label=f"PDD spectrum as produced: {meanEnergy:.2f} MeV")
    for name, (cE, cF) in spectra.items():
        if name[0] != model: continue
        r = results[model][name[1]]; colour, ls = styles[name[1]]
        label = f"{name[1]}: {r['mean_energy_MeV']:.2f} MeV, chi2 {r['chi2']:.1f}/{r['dof']}"
        if name[1] == "gaussian":
            ax.plot(cE, cF / STEP, ls, color=colour, lw=1.5, label=label)
        elif name[1].startswith("free"):
            ax.bar(cE, cF / spacing, width=0.9 * spacing, color=colour, alpha=0.35, label=label + f" ({r['free_weights']} weights)")
        else:
            ax.stairs(density(cE, cF, edges), edges, color=colour, ls=ls, lw=1.5, label=label)
    ax.set_xlim(0, edges[-1]); ax.set_xlabel("electron energy (MeV)"); ax.set_ylabel("fluence fraction per MeV")
    ax.set_title(f"{phase} exposure, {layer} layer, {model} calibration scale: fitted spectra"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    ax = axes[1, col]
    names = [n for (m, n) in spectra if m == model]
    for k, name in enumerate(names):
        pulls = np.array(results[model][name]["pulls"])
        ax.bar(np.arange(nChannels) + (k - 0.5 * (len(names) - 1)) * 0.8 / len(names), pulls, 0.8 / len(names), color=styles[name][0], label=name)
    for level in (-2, 0, 2): ax.axhline(level, color="grey" if level else "black", ls="--" if level else "-", lw=0.8)
    ax.set_xticks(range(nChannels)); ax.set_xlabel("channel"); ax.set_ylabel("pull (model - measured) / sigma"); ax.grid(alpha=0.25, axis="y"); ax.legend(fontsize=7, ncol=3)
fig.tight_layout()
fig.savefig(outDir / f"spectrum-models-{layer}.pdf", bbox_inches="tight")
