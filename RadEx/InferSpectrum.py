"""Infer effective proton energies from RadFET phase-end measurements and a response kernel.

The campaign file, the measurement models and the calibration-scale groups are
RadEx.KernelFit, shared with every other kernel fit; this module holds the proton
beam models (published and fixed Gaussian widths) and the kernel folding that
FitAllChannels and AluminiumRowEnergyScan share.

The campaign is described by a JSON file kept next to the measured doses (for
Uppsala, GRAS/3D/RadEx/RadEx-Uppsala/MeasuredDoses/campaign.json), so the same
code fits any campaign of these devices. Per phase it names the measured
phase-end table, the dose column, the phase it increments on (if any), the
accelerator-log per-channel fluence, the log fluence, the nominal energy, the
searched centre-energy range and the calibration-scale models (each of which may
be given per fit space, as {"dose": ..., "shift": ...}); for the campaign
it gives the fit space, the systematic device scatter, the response exponents
and the channels that faded between phases. Relative paths resolve against the
campaign file. The channels and their device types come from the first phase's
table (columns `channel` and `group`).

Fit spaces (`fit_space`, default "dose"; --space overrides it):

    dose    Every channel enters in linear dose space: the measured dose is the
            phase-end shift inverted through the vendor calibration, and the
            model is scale * K, weighted by

                sigma^2 = (systematic_fraction * D)^2 + sigma_noise^2 + (scale * sigma_K)^2

            where sigma_noise is the endpoint-fit RMS voltage `fit_rms_mv`
            turned into dose with the local slope of the datasheet calibration,
            D_tot / (b * V_tot) at the phase end (for a phase that increments on
            a previous one, the noise of both endpoints adds in quadrature). A
            channel whose endpoint shift is below its noise enters with zero
            dose and the noise converted through the lowest-range fit.

    shift   Every channel enters as its phase-end threshold-voltage SHIFT, and
            the model is

                dV = A_g * ((D + D0)^n - D0^n),    D = D_prev + K * F_log

            with one amplitude A per scale group, the response exponent n and
            the knee dose D0 per device type from the campaign's
            `shift_exponent` and `shift_offset_krad`, and D_prev the previous
            phase's unscaled model dose at its fitted energy (so a phase on
            pre-dosed devices is fitted as a total shift, with no dose-space
            subtraction and no inverted calibration). The response is linear
            below D0 and the power law above it; D0 = 0 is the pure power law.
            Defaults: VT01 n = 0.51, the 64 MeV phase's own fitted exponent on
            the sensor board, and D0 = 0.36 krad, which puts the vendor's 1 and
            10 krad calibration points in their published ratio at that
            exponent; VT05 n = 0.64 (the datasheet exponent, all the single
            VT05 channel with 64 MeV dose allows) and D0 = 0. The weights are

                sigma^2 = noise^2 + (d(dV)/dD * D)^2 * (systematic_fraction^2 + (sigma_K / D)^2)

            with noise = `fit_rms_mv` of the phase's own endpoint. Nothing of
            the vendor calibration enters the fit. Amplitudes are always split
            by device type in this space, because the two types' responses
            have different exponents and amplitudes: `shared` becomes
            `per-type` and `fading` becomes `fading-per-type`.

            `shift_response_table` (optional, type -> CSV with dose_krad and
            shift_v, relative to the campaign file) replaces the power law of
            that device type by the tabulated curve, log-log interpolated
            (Dependencies.EnergyKernel.tabulatedResponse). With the vendor
            calibration as the table the amplitude is the ratio of the
            devices to the calibration, and the response carries the curve's
            own falling slope instead of one exponent.

Calibration scales are solved by iteratively reweighted least squares. Each
phase's own per-channel delivered fluence multiplies the kernel.

Calibration-scale models (each phase's `scale_model` is its default):

    shared           one scale for all channels.
    fading           two scales: the campaign's `faded_channels` and all others.
                     At Uppsala these are the pre-dosed VT01, which lost 31-48 mV
                     in the 15 min beam-off pause before the 85 MeV phase,
                     against no drift for VT05 Ch7 (pre-dosed, 208 mV) and the
                     fresh channels. The grouping comes from the pause, not from
                     the fit.
    per-type         one scale per device type (VT01: Co-60, rad(Si); VT05:
                     X-ray, rad(air)).
    fading-per-type  faded and unfaded channels split again by device type: at
                     Uppsala the pre-dosed VT01 (Ch0-6), the fresh VT01 (Ch10)
                     and the VT05 (Ch7-9, 11). Ch10 reproduces the 64 MeV
                     fresh-VT01 scale (0.29 against 0.31), so it belongs with
                     neither the faded VT01 nor the VT05.

Three beam models are compared:

1. a monoenergetic beam;
2. a Gaussian energy distribution over (mu, sigma), truncated to the kernel range;
3. a Gaussian with the published width at each centre energy (as in FitAllChannels).

Optionally, beside them (off by default, so the published width stays the
reference):

4. --fixed-width-mev W: a Gaussian of fixed width W at every centre energy (at
   Skandion the facility's quoted energy-spread ceiling, 0.3 MeV), written under
   `fixed_width` in the summary and as `fixed_width_*` columns;
5. --fixed-width-cut-sigma C (with 4): the same Gaussian cut off above
   mu + C sigma, the facility's description of a spectrum targeted at the
   Bragg-peak position, under `fixed_width_truncated`. Its mu is the Gaussian's
   centre; `spectrum_mean_energy_MeV` is the mean of the cut spectrum.

Like the published width, each fixed-width model is a mu scan of the phase's
default scale model; in shift space a phase on pre-dosed devices starts from the
previous phase's model dose of the SAME fixed-width model, so each model is a
self-consistent beam hypothesis over the run.

The summary reports the default model's fits and, under "diagnostics", those of
the phase's diagnostic models; the figure shows the default fits and marks the
diagnostic Gaussians.

Intervals come from delta-chi2 <= 1, profiled over the other parameter. The
searched centre energies default to each phase's `mu_range_MeV` (at Uppsala the
nominal setting +-3 MeV, 61-67 and 82-88 MeV), which also sets the energy axis
of the misfit maps, and can be changed with --mu-range PHASE LO HI.

The folding and the scale and amplitude solvers are Dependencies.EnergyKernel; the campaign
file, the measurement models and the scale groups are RadEx.KernelFit, shared with the
electron fits; this module holds the proton beam models and the fits.

Run from the repository root:

    python3 -m RadEx.InferSpectrum --kernel <energy-kernel-dose-summary-oxide.csv> \\
        --campaign <campaign.json> --output <stem> [--space dose|shift] \\
        [--fixed-width-mev 0.3 [--fixed-width-cut-sigma 1.0]]

writes <stem>-summary.json, <stem>-channel-residuals.csv and <stem>-diagnostics.pdf.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PYTHON_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_ROOT))

from Dependencies import EnergyKernel as ek  # noqa: E402
from RadEx.KernelFit import (  # noqa: E402
    FINE_STEP, MARKER, SIGMAS, add_mu_range_argument, add_space_argument, channel_weights, gaussian_grid,
    load_campaign, load_kernel, measurement, mu_range_overrides, phase_models, scale_groups, scale_report, solve)


DIAGNOSTIC_COLOURS = ("grey", "magenta", "green")
#: Colours of the optional fixed-width models (full Gaussian, cut Gaussian).
FIXED_WIDTH_COLOURS = ("darkorange", "saddlebrown")


def published_sigma(energy_mev):
    """Skandion energy spread [MeV] at a centre energy [MeV] (Dependencies.EnergyKernel)."""
    return ek.publishedSkandionSigma(energy_mev)


def fixed_width_specs(width_mev: float | None, cut_sigma: float | None = None) -> list[dict]:
    """The optional fixed-width beam models: none, the Gaussian, and the Gaussian cut above mu + cut sigma.

    Each spec has 'key' (summary key and column prefix), 'sigma' [MeV], 'cut'
    (None or sigmas above the centre), 'label' and 'colour'.
    """
    if width_mev is None:
        if cut_sigma is not None:
            raise SystemExit("--fixed-width-cut-sigma needs --fixed-width-mev")
        return []
    if width_mev <= 0 or (cut_sigma is not None and cut_sigma <= 0):
        raise SystemExit("--fixed-width-mev and --fixed-width-cut-sigma must be positive")
    specs = [{"key": "fixed_width", "sigma": float(width_mev), "cut": None,
              "label": f"{width_mev:g} MeV width", "colour": FIXED_WIDTH_COLOURS[0]}]
    if cut_sigma is not None:
        specs.append({"key": "fixed_width_truncated", "sigma": float(width_mev), "cut": float(cut_sigma),
                      "label": f"{width_mev:g} MeV width cut at +{cut_sigma:g} sigma",
                      "colour": FIXED_WIDTH_COLOURS[1]})
    return specs


def spectrum_mean(mu: float, sigma: float, cut: float | None) -> float:
    """Mean energy [MeV] of a Gaussian (mu, sigma), cut off above mu + cut * sigma if cut is given."""
    if cut is None:
        return float(mu)
    phi = np.exp(-0.5 * cut ** 2) / np.sqrt(2.0 * np.pi)
    big_phi = 0.5 * (1.0 + math.erf(cut / np.sqrt(2.0)))
    return float(mu - sigma * phi / big_phi)


def add_width_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--fixed-width-mev", type=float, default=None,
                        help="also fit a Gaussian of this fixed energy spread [MeV] beside the published "
                             "width (e.g. 0.3, the facility's ceiling); off by default")
    parser.add_argument("--fixed-width-cut-sigma", type=float, default=None,
                        help="with --fixed-width-mev, also fit that Gaussian cut off above "
                             "mu + this many sigma (e.g. 1.0)")


def fold(energies, kernel, error, mus, sigma=None, cut=None):
    """Fold weighted kernel columns with the published width at each centre energy.

    A given `sigma` [MeV] replaces the published width by a fixed one; `cut`
    cuts each Gaussian off above mu + cut * sigma (Dependencies.EnergyKernel.foldKernel).
    """
    mus = np.asarray(mus, dtype=float)
    sigmas = ek.publishedSkandionSigma(mus) if sigma is None else float(sigma)
    return ek.foldKernel(energies, kernel, error, mus, sigmas, step=FINE_STEP, upperCut=cut)


def weighted_kernel(kernel_path: Path, campaign: dict, phase: str, channels: list[int]):
    """Kernel energies, doses and errors for `channels`, times each channel's fluence weight."""
    energies, kernel, error = load_kernel(kernel_path)
    w = channel_weights(campaign, phase, channels)
    return energies, kernel[:, channels] * w, error[:, channels] * w


def fit_model(meas, energies, k, e, mus, mu_range, groups, log_fluence):
    """Mono and Gaussian fits for one calibration-scale model; returns (summary, arrays)."""
    lo, hi = mu_range
    inside = (energies >= lo) & (energies <= hi)
    if not inside.any():
        raise ValueError(f"no kernel energy inside {lo}-{hi} MeV")
    mono_chi2 = np.array([solve(meas, k[i], e[i], groups)[1] if inside[i] else np.inf
                          for i in range(len(energies))])
    m = int(np.argmin(mono_chi2))
    mono_scales, _, mono_pulls, mono_predicted = solve(meas, k[m], e[m], groups)

    means, errs = gaussian_grid(energies, k, e, mus, SIGMAS)
    chi2 = np.empty((len(mus), len(SIGMAS)))
    for i in range(len(mus)):
        for j in range(len(SIGMAS)):
            chi2[i, j] = solve(meas, means[i, j], errs[i, j], groups)[1]
    i, j = np.unravel_index(int(np.argmin(chi2)), chi2.shape)
    scales, best_chi2, pulls, predicted = solve(meas, means[i, j], errs[i, j], groups)
    mu_ok = mus[chi2.min(axis=1) <= best_chi2 + 1.0]
    sigma_ok = SIGMAS[chi2.min(axis=0) <= best_chi2 + 1.0]

    pw_means, pw_errs = fold(energies, k, e, mus)
    pw_chi2 = np.array([solve(meas, pw_means[p], pw_errs[p], groups)[1] for p in range(len(mus))])
    p = int(np.argmin(pw_chi2))
    pw_scales, _, pw_pulls, pw_predicted = solve(meas, pw_means[p], pw_errs[p], groups)
    pw_ok = mus[pw_chi2 <= pw_chi2[p] + 1.0]

    n = len(meas["value"])
    summary = {
        "mono": {
            "energy_MeV": float(energies[m]), "chi2": float(mono_chi2[m]), "dof": n - 1 - len(groups),
            **scale_report(meas, mono_scales, log_fluence),
        },
        "gaussian": {
            "mean_energy_MeV": float(mus[i]), "sigma_energy_MeV": float(SIGMAS[j]),
            "chi2": float(best_chi2), "dof": n - 2 - len(groups),
            "delta_chi2_vs_mono": float(best_chi2 - mono_chi2[m]),
            "mean_energy_interval_MeV": [float(mu_ok.min()), float(mu_ok.max())],
            "sigma_energy_interval_MeV": [float(sigma_ok.min()), float(sigma_ok.max())],
            **scale_report(meas, scales, log_fluence),
        },
        "published_width": {
            "mean_energy_MeV": float(mus[p]), "sigma_energy_MeV": published_sigma(float(mus[p])),
            "chi2": float(pw_chi2[p]), "dof": n - 1 - len(groups),
            "mean_energy_interval_MeV": [float(pw_ok.min()), float(pw_ok.max())],
            **scale_report(meas, pw_scales, log_fluence),
        },
    }
    arrays = {
        "objective": chi2,
        "mono_predicted": mono_predicted, "mono_pull": mono_pulls,
        "gaussian_predicted": predicted, "gaussian_pull": pulls,
        "published_predicted": pw_predicted, "published_pull": pw_pulls,
        # The unscaled model dose [rad] at the published-width best energy: what
        # a later phase on the same devices starts from in shift space.
        "published_model_rad": pw_means[p] * log_fluence,
    }
    return summary, arrays


def fit_fixed_width(meas, energies, k, e, mus, groups, log_fluence, spec):
    """Centre-energy scan of one fixed-width model (fixed_width_specs); returns (summary, arrays)."""
    means, errs = fold(energies, k, e, mus, spec["sigma"], spec["cut"])
    chi2 = np.array([solve(meas, means[p], errs[p], groups)[1] for p in range(len(mus))])
    p = int(np.argmin(chi2))
    scales, _, pulls, predicted = solve(meas, means[p], errs[p], groups)
    ok = mus[chi2 <= chi2[p] + 1.0]
    summary = {
        "label": spec["label"],
        "mean_energy_MeV": float(mus[p]), "sigma_energy_MeV": spec["sigma"],
        "upper_cut_sigma": spec["cut"],
        "spectrum_mean_energy_MeV": spectrum_mean(float(mus[p]), spec["sigma"], spec["cut"]),
        "chi2": float(chi2[p]), "dof": len(meas["value"]) - 1 - len(groups),
        "mean_energy_interval_MeV": [float(ok.min()), float(ok.max())],
        "at_edge": bool(p in (0, len(mus) - 1)),
        **scale_report(meas, scales, log_fluence),
    }
    arrays = {"chi2": chi2, "predicted": predicted, "pull": pulls, "model_rad": means[p] * log_fluence}
    return summary, arrays


def fit_phase(
    campaign: dict,
    label: str,
    kernel_path: Path,
    mu_range: tuple[float, float],
    previous_model=None,
    fixed_widths=(),
    previous_fixed=None,
) -> tuple[dict, pd.DataFrame, dict]:
    """Fits of one phase to all channels: default scale model plus diagnostics.

    `previous_model` is the previous phase's unscaled model dose per channel
    [rad] at its fitted energy; it enters the shift-space model of a phase on
    pre-dosed devices and is ignored in dose space. `fixed_widths` are the
    optional fixed-width models (fixed_width_specs), fitted with the default
    scale model; `previous_fixed` gives, by spec key, the previous phase's
    model dose of that same model (`previous_model` when missing).
    """
    cfg = campaign["phases"][label]
    channels = campaign["channels"]
    meas = measurement(campaign, label, previous_model)
    energies, k, e = weighted_kernel(kernel_path, campaign, label, channels)
    weights = channel_weights(campaign, label, channels)
    lo, hi = mu_range
    mus = np.arange(lo, hi + FINE_STEP / 2, FINE_STEP)
    default, *diagnostic_models = phase_models(campaign, label)
    fits = {model: fit_model(meas, energies, k, e, mus, mu_range,
                             scale_groups(campaign, channels, model), cfg["log_fluence"])
            for model in (default, *diagnostic_models)}
    main, main_arrays = fits[default]
    fixed = {}
    for spec in fixed_widths:
        previous = (previous_fixed or {}).get(spec["key"], previous_model)
        fixed_meas = measurement(campaign, label, previous)
        fixed[spec["key"]] = (spec, fixed_meas, *fit_fixed_width(
            fixed_meas, energies, k, e, mus, scale_groups(campaign, channels, default), cfg["log_fluence"], spec))

    unit = meas["unit"]
    scale = 1000.0 if unit == "rad" else 1.0
    shown = "krad" if unit == "rad" else "V"
    sigma = meas["sigma"] if meas["space"] == "dose" else meas["noise"]
    table = cfg["table"]
    rows = []
    for c in channels:
        row = {
            "phase": label, "channel": c, "group": table.loc[c, "group"], "sample": table.loc[c, "sample"],
            "fit_space": meas["space"],
            "delivered_fluence_protons_cm2": weights[c] * cfg["log_fluence"],
            f"measured_{shown}": meas["value"][c] / scale, f"sigma_{shown}": sigma[c] / scale,
            f"mono_predicted_{shown}": main_arrays["mono_predicted"][c] / scale,
            "mono_pull": main_arrays["mono_pull"][c],
            f"gaussian_predicted_{shown}": main_arrays["gaussian_predicted"][c] / scale,
            "gaussian_pull": main_arrays["gaussian_pull"][c],
            f"published_width_predicted_{shown}": main_arrays["published_predicted"][c] / scale,
            "published_width_pull": main_arrays["published_pull"][c],
        }
        if meas["space"] == "shift":
            row["previous_model_dose_krad"] = meas["previous"][c] / 1000.0
            row["exponent"] = meas["exponent"][c]
            row["offset_krad"] = meas["offset"][c] / 1000.0
        for model in diagnostic_models:
            key = model.replace("-", "_")
            row[f"{key}_gaussian_predicted_{shown}"] = fits[model][1]["gaussian_predicted"][c] / scale
            row[f"{key}_gaussian_pull"] = fits[model][1]["gaussian_pull"][c]
        for key, (spec, fixed_meas, _, fixed_arrays) in fixed.items():
            row[f"{key}_predicted_{shown}"] = fixed_arrays["predicted"][c] / scale
            row[f"{key}_pull"] = fixed_arrays["pull"][c]
            row[f"{key}_model_dose_krad"] = fixed_arrays["model_rad"][c] / 1000.0
            if meas["space"] == "shift":
                row[f"{key}_previous_model_dose_krad"] = fixed_meas["previous"][c] / 1000.0
        if fixed:
            row["published_width_model_dose_krad"] = main_arrays["published_model_rad"][c] / 1000.0
        rows.append(row)
    summary = {
        "phase": label,
        "fit_space": meas["space"],
        "systematic_fraction": campaign["systematic"],
        "noise_floor_mv": 1e3 * campaign["noise_floor"],
        "channels": channels,
        "nominal_energy_MeV": cfg["nominal_energy"],
        "mu_range_MeV": [lo, hi],
        "scale_model": default,
        "scale_groups": {g: [channels[i] for i in idx]
                         for g, idx in scale_groups(campaign, channels, default).items()},
        **main,
        **{key: fit[2] for key, fit in fixed.items()},
        "diagnostics": {model: fits[model][0] for model in diagnostic_models},
    }
    if meas["space"] == "shift":
        types = sorted(set(campaign["device"].values()))
        summary["exponent"] = {t: campaign["exponent"][t] for t in types}
        summary["offset_krad"] = {t: campaign["offset_krad"][t] for t in types}
    arrays = {
        "mus": mus, "objective": main_arrays["objective"], "unit": shown,
        "channels": np.asarray(channels), "device": [campaign["device"][c] for c in channels],
        "measured": meas["value"] / scale, "sigma": sigma / scale,
        "mono_predicted": main_arrays["mono_predicted"] / scale,
        "gaussian_predicted": main_arrays["gaussian_predicted"] / scale,
        "published_predicted": main_arrays["published_predicted"] / scale,
        "published_model_rad": main_arrays["published_model_rad"],
        "diagnostic_gaussian_predicted": {
            model: fits[model][1]["gaussian_predicted"] / scale for model in diagnostic_models},
        # Optional fixed-width models: spec, prediction (shown unit) and the unscaled
        # model dose [rad] a later phase of the same model starts from.
        "fixed": {key: {"spec": spec, "predicted": fixed_arrays["predicted"] / scale,
                        "model_rad": fixed_arrays["model_rad"]}
                  for key, (spec, _, _, fixed_arrays) in fixed.items()},
    }
    return summary, pd.DataFrame(rows), arrays


def plot_results(results: dict[str, tuple[dict, pd.DataFrame, dict]], output: Path) -> None:
    fig, axes = plt.subplots(2, len(results), figsize=(6.0 * len(results), 8.2),
                             constrained_layout=True, squeeze=False)
    for column, (label, (summary, rows, arrays)) in enumerate(results.items()):
        ax = axes[0, column]
        x = arrays["channels"]
        unit = arrays["unit"]
        device = np.asarray(arrays["device"])
        for dev in sorted(set(device)):
            sel = device == dev
            y, s = arrays["measured"][sel], arrays["sigma"][sel]
            ax.errorbar(x[sel], y, yerr=[np.minimum(s, 0.95 * np.abs(y)), s], fmt=MARKER[dev], color="C0",
                        capsize=3, label=f"RadFET measured, {dev}")
        mono, gauss, pw = summary["mono"], summary["gaussian"], summary["published_width"]
        model = summary["scale_model"]
        curves = [
            (arrays["mono_predicted"], "s", "C1", "--",
             f"best monoenergy {mono['energy_MeV']:.2f} MeV, chi2 {mono['chi2']:.1f}/{mono['dof']}"),
            (arrays["gaussian_predicted"], "^", "C2", "--",
             f"best Gaussian {gauss['mean_energy_MeV']:.2f} MeV, sigma {gauss['sigma_energy_MeV']:.2f}, "
             f"chi2 {gauss['chi2']:.1f}/{gauss['dof']}"),
            (arrays["published_predicted"], "D", "C3", "-",
             f"best Gaussian, published width: {pw['mean_energy_MeV']:.2f} MeV, "
             f"sigma {pw['sigma_energy_MeV']:.3f}, chi2 {pw['chi2']:.1f}/{pw['dof']}"),
        ]
        for (diag_model, diag), colour in zip(summary["diagnostics"].items(), DIAGNOSTIC_COLOURS):
            g = diag["gaussian"]
            curves.append((arrays["diagnostic_gaussian_predicted"][diag_model], "x", colour, ":",
                           f"diagnostic, {diag_model} scale: {g['mean_energy_MeV']:.2f} MeV, "
                           f"chi2 {g['chi2']:.1f}/{g['dof']}"))
        for key, fixed in arrays["fixed"].items():
            f = summary[key]
            curves.append((fixed["predicted"], "v", fixed["spec"]["colour"], "-.",
                           f"best Gaussian, {fixed['spec']['label']}: {f['mean_energy_MeV']:.2f} MeV, "
                           f"chi2 {f['chi2']:.1f}/{f['dof']}"))
        for y, marker, colour, style, name in curves:
            ax.plot(x, np.where(y > 0, y, np.nan), marker + style, color=colour, mfc="none", label=name)
        ax.set_yscale("log")
        ax.set_xticks(x, [str(v) for v in x])
        ax.set_xlabel("RadFET channel")
        ax.set_ylabel("Phase dose (krad)" if unit == "krad" else "Phase-end threshold shift (V)")
        ax.set_title(f"{label.replace('mev', ' MeV phase')}: all {len(x)} channels, {model} scale, "
                     f"{summary['fit_space']} space")
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=8)

        ax = axes[1, column]
        delta = arrays["objective"] - arrays["objective"].min()
        image = ax.pcolormesh(arrays["mus"], SIGMAS, delta.T, shading="auto",
                              cmap="viridis_r", vmin=0, vmax=25, rasterized=True)
        ax.plot(arrays["mus"], [published_sigma(mu) for mu in arrays["mus"]], "w--", lw=1,
                label="published width")
        ax.plot(gauss["mean_energy_MeV"], gauss["sigma_energy_MeV"], "r*", markersize=12,
                label="best Gaussian")
        ax.plot(pw["mean_energy_MeV"], pw["sigma_energy_MeV"], "D", color="C3", mfc="none", markersize=8,
                mew=2, label="best Gaussian, published width")
        for (diag_model, diag), colour in zip(summary["diagnostics"].items(), DIAGNOSTIC_COLOURS):
            g = diag["gaussian"]
            ax.plot(g["mean_energy_MeV"], g["sigma_energy_MeV"], "x", color=colour, markersize=9, mew=2,
                    label=f"{diag_model} diagnostic")
        lines = set()
        for key, fixed in arrays["fixed"].items():
            spec, f = fixed["spec"], summary[key]
            if spec["sigma"] not in lines:
                lines.add(spec["sigma"])
                ax.axhline(spec["sigma"], color=spec["colour"], ls=":", lw=1.5,
                           label=f"sigma = {spec['sigma']:g} MeV (fixed width)")
            ax.plot(f["mean_energy_MeV"], spec["sigma"], "v", color=spec["colour"], mfc="none", markersize=9,
                    mew=2, label=f"best Gaussian, {spec['label']}"
                    + (" (map: full Gaussians)" if spec["cut"] is not None else ""))
        ax.set_xlim(*summary["mu_range_MeV"])
        ax.set_xlabel("Gaussian mean energy (MeV)")
        ax.set_ylabel("Gaussian sigma (MeV)")
        ax.set_title(f"Weighted 12-channel misfit, {model} scale")
        ax.legend(fontsize=8, loc="upper right")
        fig.colorbar(image, ax=ax, label="Delta chi2 (clipped at 25)")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--kernel", type=Path, required=True,
                        help="unweighted kernel summary from BuildKernelSummary")
    parser.add_argument("--campaign", type=Path, required=True,
                        help="campaign JSON next to the measured doses")
    parser.add_argument("--output", type=Path, required=True,
                        help="output stem, e.g. Results/spectrum-fit-oxide")
    add_mu_range_argument(parser)
    add_space_argument(parser)
    add_width_argument(parser)
    args = parser.parse_args()
    specs = fixed_width_specs(args.fixed_width_mev, args.fixed_width_cut_sigma)
    campaign = load_campaign(args.campaign, args.space, args.systematic, args.noise_floor_mv)
    mu_ranges = mu_range_overrides(campaign, args.mu_range)
    stem = args.output
    stem.parent.mkdir(parents=True, exist_ok=True)

    # Phases in order: a phase on pre-dosed devices starts, in shift space,
    # from the previous phase's unscaled model dose at its fitted energy
    # (each fixed-width model from its own previous model dose).
    results = {}
    model_by_phase = {}
    fixed_by_phase = {}
    for label, cfg in campaign["phases"].items():
        previous = model_by_phase.get(cfg["previous_phase"]) if cfg["previous_phase"] else None
        previous_fixed = fixed_by_phase.get(cfg["previous_phase"]) if cfg["previous_phase"] else None
        results[label] = fit_phase(campaign, label, args.kernel, mu_ranges[label], previous,
                                   specs, previous_fixed)
        model_by_phase[label] = results[label][2]["published_model_rad"] + (
            previous if previous is not None else 0.0)
        fixed_by_phase[label] = {
            key: fixed["model_rad"] + (previous_fixed[key] if previous_fixed else 0.0)
            for key, fixed in results[label][2]["fixed"].items()}

    summary = {label: result[0] for label, result in results.items()}
    stem.with_name(f"{stem.name}-summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    pd.concat([result[1] for result in results.values()], ignore_index=True).to_csv(
        stem.with_name(f"{stem.name}-channel-residuals.csv"), index=False
    )
    plot_results(results, stem.with_name(f"{stem.name}-diagnostics.pdf"))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
