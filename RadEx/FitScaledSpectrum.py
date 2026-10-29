"""Beam energy from a campaign's measured RadFET doses and a mono-energetic dose kernel.

    python3 -m RadEx.FitScaledSpectrum --kernel-dir <sweep> --campaign <json> --output-dir <dir>
        [--layer oxide|die] [--exclude 9,10] [--systematic 0.10] [--tag name]

Measured doses and errors, the device types and the calibration-scale groups come from
the campaign file through RadEx.KernelFit, the core shared with the proton fits; every
phase in it is fitted. The kernel is read from the sweep's TotalDose tables with
Dependencies.EnergyKernel.readKernelRuns, in dose per unit of the runs' normalisation
fluence; --layer picks its tally block from the campaign. Beam models: mono-energetic,
Gaussian (mu, sigma), and the campaign's production spectrum with every component
energy scaled by alpha. Calibration scales shared or per device type; intervals are
delta-chi2 <= 1, profiled. --exclude leaves channels out of the chi2 (they are still
predicted and plotted); --tag suffixes the output names so variants sit beside the
default fit. Writes spectrum-fit-<layer>[-tag]-{summary.json,channel-residuals.csv,
<phase>.pdf,<phase>-chi2.pdf} into --output-dir.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from Dependencies.EnergyKernel import (chi2Interval, foldComponents, foldKernel, interpolateKernel,
                                       readKernelRuns, scanScales)
from Plotting.KernelFitPlots import plotChannelFit, plotChi2Scans
from RadEx import KernelFit as kf
from Read.ReadGPSSources import readGPSSources

STEP, ALPHA_STEP = 0.025, 0.0025
SIGMAS = np.arange(0.05, 4.0001, 0.05)

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--kernel-dir", type=Path, required=True)
parser.add_argument("--campaign", type=Path, required=True)
parser.add_argument("--output-dir", type=Path, required=True)
parser.add_argument("--layer", help="detector layer of the campaign (default: its default_layer)")
parser.add_argument("--exclude", default="", help="channels left out of the chi2, e.g. 9,10")
parser.add_argument("--systematic", type=float, help="model fraction in sigma (default: the campaign's)")
parser.add_argument("--tag", default="", help="suffix for the output names")
parser.add_argument("--mu-min", type=float, default=3.0, help="lowest mono and Gaussian centre [MeV]")
parser.add_argument("--alpha-min", type=float, default=0.80, help="lowest spectrum scale")
args = parser.parse_args()
default = kf.load_campaign(args.campaign)["systematic"]
campaign = kf.load_campaign(args.campaign, systematic=args.systematic)
layer = args.layer or campaign["default_layer"]
systematic = campaign["systematic"]
channels = campaign["channels"]
TYPE = campaign["device"]
EXCLUDED = sorted(int(c) for c in args.exclude.split(",") if c.strip())
TAG = ("-" + args.tag) if args.tag else ""
if (EXCLUDED or systematic != default) and not TAG:
    parser.error("--exclude or a non-default --systematic needs --tag")
out = args.output_dir
out.mkdir(parents=True, exist_ok=True)

energies, dose, error, _ = readKernelRuns(args.kernel_dir)
kernel, kerr = dose[:, kf.layer_slice(campaign, layer)], error[:, kf.layer_slice(campaign, layer)]
gps = readGPSSources(str(campaign["spectrum"]))
meanEnergy = float(gps["fraction"] @ gps["energy"])
mus = np.arange(args.mu_min, energies.max() + STEP / 2, STEP)
alphas = np.arange(args.alpha_min, energies.max() / gps["energy"].max() + 1e-9, ALPHA_STEP)
models = {"mono": (interpolateKernel(energies, kernel, mus), interpolateKernel(energies, kerr, mus)),
          "scaled": foldComponents(energies, kernel, kerr, gps["energy"], gps["fraction"], alphas)}
gaussian = [foldKernel(energies, kernel, kerr, mus, s, STEP) for s in SIGMAS]
groupsOf = {model: kf.scale_groups(campaign, channels, model, EXCLUDED) for model in ("shared", "per-type")}


def complete(fit, groups):
    """Predict and pull the excluded channels with their group's scale (not in chi2)."""
    if not EXCLUDED:
        return fit
    scale = np.ones(len(channels))
    for name in groups:
        members = [c for c in channels if (name == "all" or TYPE[c] == name)]
        scale[members] = fit["scales"][name]
    predicted = fit["predicted"].copy()
    pulls = fit["pulls"].copy()
    predicted[EXCLUDED] = scale[EXCLUDED] * fit["model"][EXCLUDED]
    pulls[EXCLUDED] = (predicted[EXCLUDED] - fit["measured"][EXCLUDED]) / fit["sigma"][EXCLUDED]
    return dict(fit, predicted=predicted, pulls=pulls)


def scan(measured, sigma, models, errors, groups):
    fit = scanScales(measured, sigma, models, errors, groups)
    fit.update(model=models[fit["best"]], measured=measured, sigma=sigma)
    return complete(fit, groups)


summary = {"layer": layer, "systematic": systematic, "excluded_channels": EXCLUDED,
           "kernel_energies_MeV": energies.tolist(), "campaign": str(campaign["path"])}
residuals = []
measurement = campaign["measurement_phase"]
for phase, cfg in sorted(campaign["phases"].items(), key=lambda p: p[0] != measurement):
    meas = kf.measurement(campaign, phase)
    measured, sigma, fluence, table = meas["value"], meas["sigma"], cfg["log_fluence"], cfg["table"]
    summary[phase] = {(cfg["log_fluence_column"] or "log_fluence").removeprefix("phase_"): fluence,
                      "note": str(table.loc[channels[0]].get("start_time", "")) + " to "
                      + str(table.loc[channels[0]].get("end_time", ""))}
    if phase != measurement:
        summary[phase]["note"] += "; " + cfg["note"]
    scans, predictions = [], []
    for model, groups in groupsOf.items():
        result = {}
        mono = scan(measured, sigma, *models["mono"], groups)
        best, lo, hi, edge = chi2Interval(mus, mono["chi2"])
        result["mono"] = {"energy_MeV": best, "interval_MeV": [lo, hi], "at_edge": edge,
                          "chi2": float(mono["chi2"][mono["best"]]), "dof": mono["dof"] - 1}
        grid = np.array([scanScales(measured, sigma, g[0], g[1], groups)["chi2"] for g in gaussian])  # (sigma, mu)
        j, i = np.unravel_index(int(np.nanargmin(grid)), grid.shape)
        best, lo, hi, edge = chi2Interval(mus, grid.min(axis=0))
        result["gaussian"] = {"mean_energy_MeV": float(mus[i]), "sigma_energy_MeV": float(SIGMAS[j]),
                              "mean_interval_MeV": [lo, hi], "at_edge": edge, "chi2": float(grid[j, i]),
                              "dof": mono["dof"] - 2}
        scaled = scan(measured, sigma, *models["scaled"], groups)
        best, lo, hi, edge = chi2Interval(alphas, scaled["chi2"])
        k = int(np.argmin(np.abs(alphas - 1.0)))
        produced = scan(measured, sigma, models["scaled"][0][[k]], models["scaled"][1][[k]], groups)
        result["scaled_pdd"] = {"alpha": best, "alpha_interval": [lo, hi], "at_edge": edge,
                                "mean_energy_MeV": best * meanEnergy,
                                "mean_energy_interval_MeV": [lo * meanEnergy, hi * meanEnergy],
                                "chi2": float(scaled["chi2"][scaled["best"]]), "dof": scaled["dof"] - 1,
                                "calibration_ratio": {g: s / fluence for g, s in scaled["scales"].items()},
                                "as_produced": {"chi2": float(produced["chi2"][0]), "dof": produced["dof"],
                                                "calibration_ratio": {g: s / fluence for g, s in produced["scales"].items()}}}
        summary[phase][model] = result
        for c in channels:
            residuals.append({"phase": phase, "scale_model": model, "channel": c, "type": TYPE[c],
                              "in_fit": c not in EXCLUDED, "measured_krad": measured[c], "sigma_krad": sigma[c],
                              "scaled_predicted_krad": scaled["predicted"][c], "scaled_pull": scaled["pulls"][c],
                              "produced_predicted_krad": produced["predicted"][c], "produced_pull": produced["pulls"][c],
                              "mono_predicted_krad": mono["predicted"][c], "mono_pull": mono["pulls"][c]})
        s = result["scaled_pdd"]
        predictions.append((scaled["predicted"], f"{model}: spectrum x {s['alpha']:.4f} ({s['mean_energy_MeV']:.2f} MeV), "
                            f"chi2 {s['chi2']:.1f}/{s['dof']}", scaled["pulls"]))
        predictions.append((produced["predicted"], f"{model}: as produced ({meanEnergy:.2f} MeV), "
                            f"chi2 {s['as_produced']['chi2']:.1f}/{s['as_produced']['dof']}", None))
        scans.append((alphas * meanEnergy, [(scaled["chi2"], f"{model}: spectrum scaled")],
                      f"mean energy of the scaled spectrum (MeV), {phase}", [(meanEnergy, "as produced")]))
        print(f"{phase:5s} {layer} {model:8s} mono {result['mono']['energy_MeV']:.2f} MeV chi2 {result['mono']['chi2']:.1f}/{result['mono']['dof']} | "
              f"Gaussian {result['gaussian']['mean_energy_MeV']:.2f} ({result['gaussian']['mean_interval_MeV'][0]:.2f}-"
              f"{result['gaussian']['mean_interval_MeV'][1]:.2f}) sigma {result['gaussian']['sigma_energy_MeV']:.2f} "
              f"chi2 {result['gaussian']['chi2']:.1f}/{result['gaussian']['dof']} | "
              f"spectrum x{s['alpha']:.4f} ({s['alpha_interval'][0]:.4f}-{s['alpha_interval'][1]:.4f}) = {s['mean_energy_MeV']:.2f} MeV "
              f"chi2 {s['chi2']:.1f}/{s['dof']}{' EDGE' if s['at_edge'] else ''} | as produced chi2 {s['as_produced']['chi2']:.1f} | "
              "calibration " + ", ".join(f"{g} {v:.3f}" for g, v in s["calibration_ratio"].items()))
    without = f", without Ch{', '.join(map(str, EXCLUDED))}" if EXCLUDED else ""
    plotChannelFit(channels, measured, sigma, predictions, str(out / f"spectrum-fit-{layer}{TAG}-{phase}.pdf"),
                   [TYPE[c] for c in channels],
                   f"{phase} exposure, {layer} layer{without}, model term {systematic:.0%}: measured doses and the scaled spectrum")
    plotChi2Scans(scans, str(out / f"spectrum-fit-{layer}{TAG}-{phase}-chi2.pdf"), title=f"{phase} exposure, {layer} layer")

with open(out / f"spectrum-fit-{layer}{TAG}-summary.json", "w") as handle:
    json.dump(summary, handle, indent=2)
with open(out / f"spectrum-fit-{layer}{TAG}-channel-residuals.csv", "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(residuals[0]))
    writer.writeheader()
    writer.writerows(residuals)
