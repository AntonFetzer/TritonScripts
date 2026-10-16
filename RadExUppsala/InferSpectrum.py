"""Infer effective Uppsala proton energies from all twelve RadFET phase-end doses.

This module also holds the dose model that FitAllChannels and
AluminiumRowEnergyScan share, so that every fit uses the same constants,
uncertainties and kernel folding.

Every channel enters the fit in linear dose space, weighted by

    sigma^2 = (SYSTEMATIC * D)^2 + sigma_noise^2 + (scale * sigma_K)^2

where SYSTEMATIC is the ~10 % device-to-device scatter of the VT01 channels,
sigma_noise is the endpoint-fit RMS voltage `fit_rms_mv` turned into dose with
the local slope of the datasheet calibration, D_tot / (b * V_tot) at the phase
end (for the 85 MeV increment the noise of both endpoints adds in quadrature),
and sigma_K is the kernel's Monte Carlo error. One calibration scale, shared by
all twelve channels, is solved by iteratively reweighted least squares; both
device types are converted to dose with their datasheet calibrations, so a
single constant relates measured and simulated dose. As a diagnostic the same
fits are repeated with one scale per device type (VT01: Co-60, rad(Si); VT05:
X-ray, rad(air)), which absorbs any VT05/VT01 calibration mismatch. Each phase's
own per-channel delivered fluence (accelerator-log radfet_location_fluence.csv)
multiplies the kernel.

Two beam models are compared:

1. a monoenergetic beam;
2. a Gaussian energy distribution over (mu, sigma), truncated to the kernel range.

The summary reports the shared-scale fits and, under "per_type_diagnostic", the
per-type ones; the figure shows the shared fits and marks the per-type Gaussian.

Intervals come from delta-chi2 <= 1, profiled over the other parameter. The
searched centre energies default to PHASES' mu_range and can be narrowed with
--mu-range-64 and --mu-range-85.

Run from the repository root:

    python3 -m RadExUppsala.InferSpectrum --kernel <energy-kernel-dose-summary-oxide.csv> \\
        --dose-64 <phase-end-dose-64mev.csv> --dose-85 <phase-end-dose-85mev.csv> \\
        --fluence-64 <RadEx_64MeV/radfet_location_fluence.csv> \\
        --fluence-85 <Radex_85MeV/radfet_location_fluence.csv> --output-dir <dir>

writes spectrum-fit-summary.json, spectrum-fit-channel-residuals.csv and
spectrum-fit-diagnostics.png into <dir>.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PHASES = {
    "64mev": {
        "dose_column": "dose_krad",
        "mu_range": (60.0, 65.0),
        "nominal_energy": 64.0,
        "log_fluence": 6.942e10,
    },
    "85mev": {
        "dose_column": "dose_krad_85mev_increment",
        "mu_range": (80.0, 90.0),
        "nominal_energy": 85.0,
        "log_fluence": 9.152e10,
    },
}
CHANNELS = list(range(12))
DEVICE = {c: "VT01" for c in (0, 1, 2, 3, 4, 5, 6, 10)}
DEVICE.update({c: "VT05" for c in (7, 8, 9, 11)})
MARKER = {"VT01": "o", "VT05": "s"}
SYSTEMATIC = 0.10
#: Exponents B of the Varadis datasheet fits dV = A * dose_rad**B (QF 16 test
#: records, 13/05/2025, as in fs2_radex ProtonAnalysis.DATASHEET_FITS), by the
#: upper end in krad of the dose range each fit covers.
VT01_EXPONENT_BY_RANGE = ((1.0, 0.91193), (5.0, 0.70621), (10.0, 0.61638),
                          (50.0, 0.47515), (100.0, 0.45509))
VT05_EXPONENT = 0.64146
#: Step of the fine energy grid for interpolation, folding and the mu search [MeV].
FINE_STEP = 0.025
SIGMAS = np.arange(0.05, 4.0001, 0.05)


def published_sigma(energy_mev: float) -> float:
    """Skandion energy spread [MeV] at a centre energy [MeV].

    Ardenfors et al., J. Phys.: Conf. Ser. 860 (2017) 012025, as interpolated by
    the production spectrum macros: sigma_E/E falls linearly from 0.8 % at 60 MeV
    to 0.3 % at 226 MeV, which reproduces the 0.50429 and 0.61600 MeV in the
    64 and 85 MeV macros.
    """
    return energy_mev * (0.8 - 0.5 * (energy_mev - 60.0) / 166.0) / 100.0


def load_kernel(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    frame = pd.read_csv(path)
    required = {"energy_MeV", "channel", "mean_dose_rad", "combined_error_rad"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Kernel is missing columns: {sorted(missing)}")
    dose = frame.pivot(index="energy_MeV", columns="channel", values="mean_dose_rad")
    err = frame.pivot(index="energy_MeV", columns="channel", values="combined_error_rad")
    dose = dose.sort_index()
    err = err.reindex(index=dose.index, columns=dose.columns)
    if dose.isna().any().any():
        raise ValueError("Kernel is not a complete energy-by-channel matrix")
    return dose.index.to_numpy(float), dose.to_numpy(float), err.to_numpy(float)


def interpolate_kernel(
    energies: np.ndarray, kernel: np.ndarray, evaluation_energy: np.ndarray
) -> np.ndarray:
    return np.column_stack(
        [np.interp(evaluation_energy, energies, kernel[:, i]) for i in range(kernel.shape[1])]
    )


def channel_weights(fluence_path: Path, phase: str, channels: list[int]) -> np.ndarray:
    """Delivered fluence per channel relative to the phase's log fluence."""
    frame = pd.read_csv(fluence_path).set_index("channel")
    fluence = frame.loc[channels, "delivered_fluence_protons_cm2"].to_numpy(float)
    if np.any(fluence <= 0):
        raise ValueError(f"{fluence_path}: non-positive delivered fluence")
    return fluence / PHASES[phase]["log_fluence"]


def local_exponent(device: str, dose_krad: float) -> float:
    if device == "VT05":
        return VT05_EXPONENT
    for top, exponent in VT01_EXPONENT_BY_RANGE:
        if dose_krad <= top:
            return exponent
    return VT01_EXPONENT_BY_RANGE[-1][1]


def measurements(dose_64: Path, dose_85: Path, phase: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Measured dose [rad], its 1-sigma uncertainty [rad] and the noise part, per channel."""
    t64 = pd.read_csv(dose_64).set_index("channel").loc[CHANNELS]
    t85 = pd.read_csv(dose_85).set_index("channel").loc[CHANNELS]
    if phase == "64mev":
        dose = t64["dose_krad"].to_numpy(float) * 1000.0
        total = dose
        v_total = t64["shift_end_v"].to_numpy(float)
        noise_v = t64["fit_rms_mv"].to_numpy(float) / 1000.0
    else:
        dose = t85["dose_krad_85mev_increment"].to_numpy(float) * 1000.0
        total = t85["dose_krad"].to_numpy(float) * 1000.0
        v_total = t85["shift_end_v"].to_numpy(float)
        noise_v = np.hypot(t85["fit_rms_mv"].to_numpy(float), t64["fit_rms_mv"].to_numpy(float)) / 1000.0
    if np.any(v_total <= 0) or np.any(total <= 0):
        raise RuntimeError(f"{phase}: non-positive endpoint shift or dose; the noise model needs both positive")
    exponent = np.array([local_exponent(DEVICE[c], total[i] / 1000.0) for i, c in enumerate(CHANNELS)])
    noise = total / (exponent * v_total) * noise_v
    sigma = np.hypot(SYSTEMATIC * dose, noise)
    return dose, sigma, noise


#: Calibration-scale models: the shared scale is the default, per-type a diagnostic.
SCALE_MODELS = ("shared", "per-type")


def scale_groups(channels: list[int], model: str = "shared") -> dict[str, np.ndarray]:
    """Indices into `channels` that share one calibration scale, by model."""
    if model == "shared":
        return {"all": np.arange(len(channels))}
    if model != "per-type":
        raise ValueError(f"unknown scale model {model!r}; use one of {SCALE_MODELS}")
    groups = {dev: np.array([i for i, c in enumerate(channels) if DEVICE[c] == dev])
              for dev in ("VT01", "VT05")}
    return {k: v for k, v in groups.items() if len(v)}


def solve_scales(dose, sigma, model, model_err, groups):
    """Iteratively reweighted least squares for one scale per group.

    Returns (scales by group, chi2, pulls, per-channel scale).
    """
    scale = np.ones(len(dose))
    chi2 = 0.0
    pulls = np.zeros(len(dose))
    scales = {}
    for name, idx in groups.items():
        k, kerr, d, sd = model[idx], model_err[idx], dose[idx], sigma[idx]
        s = float(np.sum(d) / np.sum(k)) if np.sum(k) > 0 else 0.0
        for _ in range(30):
            w = 1.0 / (sd**2 + (s * kerr) ** 2)
            new = float(np.sum(w * d * k) / np.sum(w * k * k))
            if abs(new - s) <= 1e-12 * max(abs(s), 1e-300):
                s = new
                break
            s = new
        w = 1.0 / (sd**2 + (s * kerr) ** 2)
        resid = s * k - d
        chi2 += float(np.sum(w * resid**2))
        pulls[idx] = resid * np.sqrt(w)
        scale[idx] = s
        scales[name] = s
    return scales, chi2, pulls, scale


def _fine_grid(energies, kernel, error):
    fine = np.arange(float(energies.min()), float(energies.max()) + FINE_STEP / 2, FINE_STEP)
    k = interpolate_kernel(energies, kernel, fine)
    e = interpolate_kernel(energies, error, fine)
    # Interpolated errors are correlated over each kernel interval; scale the
    # quadrature sum by the number of fine points per kernel interval so that
    # folding does not shrink a single kernel point's error.
    points_per_interval = np.interp(fine, energies[1:], np.diff(energies)) / FINE_STEP
    return fine, k, e, points_per_interval


def fold(energies, kernel, error, mus):
    """Fold weighted kernel columns with the published width at each centre energy."""
    fine, k, e, ppi = _fine_grid(energies, kernel, error)
    means = np.empty((len(mus), kernel.shape[1]))
    errs = np.empty_like(means)
    for i, mu in enumerate(mus):
        w = np.exp(-0.5 * ((fine - mu) / published_sigma(mu)) ** 2)
        w /= w.sum()
        means[i] = w @ k
        errs[i] = np.sqrt((w**2 * ppi) @ e**2)
    return means, errs


def gaussian_grid(energies, kernel, error, mus, sigmas):
    """Fold weighted kernel columns with Gaussians on a (mu, sigma) grid."""
    fine, k, e, ppi = _fine_grid(energies, kernel, error)
    means = np.empty((len(mus), len(sigmas), kernel.shape[1]))
    errs = np.empty_like(means)
    for i, mu in enumerate(mus):
        w = np.exp(-0.5 * ((fine[None, :] - mu) / sigmas[:, None]) ** 2)
        w /= w.sum(axis=1, keepdims=True)
        means[i] = w @ k
        errs[i] = np.sqrt((w**2 * ppi) @ e**2)
    return means, errs


def weighted_kernel(kernel_path: Path, fluence_path: Path, phase: str, channels: list[int]):
    """Kernel energies, doses and errors for `channels`, times each channel's fluence weight."""
    energies, kernel, error = load_kernel(kernel_path)
    w = channel_weights(fluence_path, phase, channels)
    return energies, kernel[:, channels] * w, error[:, channels] * w


def fit_model(dose, sigma, energies, k, e, mus, mu_range, groups, log_fluence):
    """Mono and Gaussian fits for one calibration-scale model; returns (summary, arrays)."""
    lo, hi = mu_range
    inside = (energies >= lo) & (energies <= hi)
    if not inside.any():
        raise ValueError(f"no kernel energy inside {lo}-{hi} MeV")
    mono_chi2 = np.array([solve_scales(dose, sigma, k[i], e[i], groups)[1] if inside[i] else np.inf
                          for i in range(len(energies))])
    m = int(np.argmin(mono_chi2))
    mono_scales, _, mono_pulls, mono_scale = solve_scales(dose, sigma, k[m], e[m], groups)

    means, errs = gaussian_grid(energies, k, e, mus, SIGMAS)
    chi2 = np.empty((len(mus), len(SIGMAS)))
    for i in range(len(mus)):
        for j in range(len(SIGMAS)):
            chi2[i, j] = solve_scales(dose, sigma, means[i, j], errs[i, j], groups)[1]
    i, j = np.unravel_index(int(np.argmin(chi2)), chi2.shape)
    scales, best_chi2, pulls, scale = solve_scales(dose, sigma, means[i, j], errs[i, j], groups)
    mu_ok = mus[chi2.min(axis=1) <= best_chi2 + 1.0]
    sigma_ok = SIGMAS[chi2.min(axis=0) <= best_chi2 + 1.0]

    n = len(dose)
    summary = {
        "mono": {
            "energy_MeV": float(energies[m]), "chi2": float(mono_chi2[m]), "dof": n - 1 - len(groups),
            "calibration_ratio": {g: s / log_fluence for g, s in mono_scales.items()},
        },
        "gaussian": {
            "mean_energy_MeV": float(mus[i]), "sigma_energy_MeV": float(SIGMAS[j]),
            "chi2": float(best_chi2), "dof": n - 2 - len(groups),
            "delta_chi2_vs_mono": float(best_chi2 - mono_chi2[m]),
            "mean_energy_interval_MeV": [float(mu_ok.min()), float(mu_ok.max())],
            "sigma_energy_interval_MeV": [float(sigma_ok.min()), float(sigma_ok.max())],
            "calibration_ratio": {g: s / log_fluence for g, s in scales.items()},
        },
    }
    arrays = {
        "objective": chi2,
        "mono_predicted_rad": mono_scale * k[m], "mono_pull": mono_pulls,
        "gaussian_predicted_rad": scale * means[i, j], "gaussian_pull": pulls,
    }
    return summary, arrays


def fit_phase(
    label: str,
    dose_64: Path,
    dose_85: Path,
    fluence_path: Path,
    kernel_path: Path,
    mu_range: tuple[float, float],
) -> tuple[dict, pd.DataFrame, dict]:
    """Shared-scale fits of one phase to all twelve channels, plus the per-type diagnostic."""
    cfg = PHASES[label]
    channels = CHANNELS
    dose, sigma, _noise = measurements(dose_64, dose_85, label)
    energies, k, e = weighted_kernel(kernel_path, fluence_path, label, channels)
    weights = channel_weights(fluence_path, label, channels)
    lo, hi = mu_range
    mus = np.arange(lo, hi + FINE_STEP / 2, FINE_STEP)
    fits = {model: fit_model(dose, sigma, energies, k, e, mus, mu_range,
                             scale_groups(channels, model), cfg["log_fluence"])
            for model in SCALE_MODELS}
    shared, shared_arrays = fits["shared"]
    per_type, per_type_arrays = fits["per-type"]

    table = pd.read_csv(dose_64 if label == "64mev" else dose_85).set_index("channel")
    rows = [{
        "phase": label, "channel": c, "group": table.loc[c, "group"], "sample": table.loc[c, "sample"],
        "delivered_fluence_protons_cm2": weights[c] * cfg["log_fluence"],
        "measured_dose_krad": dose[c] / 1000.0, "sigma_krad": sigma[c] / 1000.0,
        "mono_predicted_krad": shared_arrays["mono_predicted_rad"][c] / 1000.0,
        "mono_pull": shared_arrays["mono_pull"][c],
        "gaussian_predicted_krad": shared_arrays["gaussian_predicted_rad"][c] / 1000.0,
        "gaussian_pull": shared_arrays["gaussian_pull"][c],
        "per_type_gaussian_predicted_krad": per_type_arrays["gaussian_predicted_rad"][c] / 1000.0,
        "per_type_gaussian_pull": per_type_arrays["gaussian_pull"][c],
    } for c in channels]
    summary = {
        "phase": label,
        "channels": channels,
        "nominal_energy_MeV": cfg["nominal_energy"],
        "mu_range_MeV": [lo, hi],
        "scale_model": "shared",
        **shared,
        "per_type_diagnostic": per_type,
    }
    arrays = {
        "mus": mus, "objective": shared_arrays["objective"],
        "channels": np.asarray(channels), "device": [DEVICE[c] for c in channels],
        "measured_krad": dose / 1000.0, "sigma_krad": sigma / 1000.0,
        "mono_predicted_krad": shared_arrays["mono_predicted_rad"] / 1000.0,
        "gaussian_predicted_krad": shared_arrays["gaussian_predicted_rad"] / 1000.0,
        "per_type_gaussian_predicted_krad": per_type_arrays["gaussian_predicted_rad"] / 1000.0,
    }
    return summary, pd.DataFrame(rows), arrays


def plot_results(results: dict[str, tuple[dict, pd.DataFrame, dict]], output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.2), constrained_layout=True)
    for column, (label, (summary, rows, arrays)) in enumerate(results.items()):
        ax = axes[0, column]
        x = arrays["channels"]
        device = np.asarray(arrays["device"])
        for dev in ("VT01", "VT05"):
            sel = device == dev
            y, s = arrays["measured_krad"][sel], arrays["sigma_krad"][sel]
            ax.errorbar(x[sel], y, yerr=[np.minimum(s, 0.95 * y), s], fmt=MARKER[dev], color="C0",
                        capsize=3, label=f"RadFET inferred, {dev}")
        mono, gauss = summary["mono"], summary["gaussian"]
        diag = summary["per_type_diagnostic"]["gaussian"]
        for key, marker, colour, style, name in (
                ("mono_predicted_krad", "s", "C1", "--",
                 f"best monoenergy {mono['energy_MeV']:.2f} MeV, chi2 {mono['chi2']:.1f}/{mono['dof']}"),
                ("gaussian_predicted_krad", "^", "C2", "--",
                 f"best Gaussian {gauss['mean_energy_MeV']:.2f} MeV, sigma {gauss['sigma_energy_MeV']:.2f}, "
                 f"chi2 {gauss['chi2']:.1f}/{gauss['dof']}"),
                ("per_type_gaussian_predicted_krad", "x", "grey", ":",
                 f"diagnostic, per-type scales: {diag['mean_energy_MeV']:.2f} MeV, "
                 f"chi2 {diag['chi2']:.1f}/{diag['dof']}")):
            y = arrays[key]
            ax.plot(x, np.where(y > 0, y, np.nan), marker + style, color=colour, mfc="none", label=name)
        ax.set_yscale("log")
        ax.set_xticks(x, [str(v) for v in x])
        ax.set_xlabel("RadFET channel")
        ax.set_ylabel("Phase dose (krad)")
        ax.set_title(f"{label.replace('mev', ' MeV phase')}: all 12 channels, shared scale")
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=8)

        ax = axes[1, column]
        delta = arrays["objective"] - arrays["objective"].min()
        image = ax.pcolormesh(arrays["mus"], SIGMAS, delta.T, shading="auto",
                              cmap="viridis_r", vmin=0, vmax=25)
        ax.plot(arrays["mus"], [published_sigma(mu) for mu in arrays["mus"]], "w--", lw=1,
                label="published width")
        ax.plot(gauss["mean_energy_MeV"], gauss["sigma_energy_MeV"], "r*", markersize=12,
                label="best Gaussian")
        ax.plot(diag["mean_energy_MeV"], diag["sigma_energy_MeV"], "wx", markersize=9, mew=2,
                label="per-type diagnostic")
        ax.set_xlabel("Gaussian mean energy (MeV)")
        ax.set_ylabel("Gaussian sigma (MeV)")
        ax.set_title("Weighted 12-channel misfit, shared scale")
        ax.legend(fontsize=8, loc="upper right")
        fig.colorbar(image, ax=ax, label="Delta chi2 (clipped at 25)")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--kernel", type=Path, required=True,
                        help="unweighted kernel summary from BuildKernelSummary")
    parser.add_argument("--dose-64", type=Path, required=True)
    parser.add_argument("--dose-85", type=Path, required=True)
    parser.add_argument("--fluence-64", type=Path, required=True)
    parser.add_argument("--fluence-85", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mu-range-64", type=float, nargs=2, default=PHASES["64mev"]["mu_range"])
    parser.add_argument("--mu-range-85", type=float, nargs=2, default=PHASES["85mev"]["mu_range"])
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    results = {}
    for label, fluence, mu_range in (("64mev", args.fluence_64, args.mu_range_64),
                                     ("85mev", args.fluence_85, args.mu_range_85)):
        results[label] = fit_phase(label, args.dose_64, args.dose_85, fluence, args.kernel,
                                   tuple(mu_range))

    summary = {label: result[0] for label, result in results.items()}
    (args.output_dir / "spectrum-fit-summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    pd.concat([result[1] for result in results.values()], ignore_index=True).to_csv(
        args.output_dir / "spectrum-fit-channel-residuals.csv", index=False
    )
    plot_results(results, args.output_dir / "spectrum-fit-diagnostics.png")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
