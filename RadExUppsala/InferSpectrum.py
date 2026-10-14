"""Infer effective Uppsala proton energies from RadFET phase-end doses.

The absolute fluence is a nuisance parameter. Fits therefore use the pattern of
log dose across trusted VT01 channels. Two models are compared:

1. a monoenergetic beam;
2. a Gaussian energy distribution, truncated to the simulated kernel range.

This is deliberately a low-dimensional, auditable inverse model. A free spectrum
over the 60 simulated energies is not identifiable from seven or eight channels.
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
        "channels": [0, 1, 2, 3, 4, 5, 6],
        "mu_range": (58.0, 70.0),
        "nominal_energy": 64.0,
        "log_fluence": 6.942e10,
    },
    "85mev": {
        "dose_column": "dose_krad_85mev_increment",
        "channels": [0, 1, 2, 3, 4, 5, 6, 10],
        "mu_range": (68.0, 92.0),
        "nominal_energy": 85.0,
        "log_fluence": 9.152e10,
    },
}


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


def fit_scale_log(measured_rad: np.ndarray, response: np.ndarray) -> tuple[float, np.ndarray, float]:
    if np.any(measured_rad <= 0) or np.any(response <= 0):
        raise ValueError("Log-pattern fit requires strictly positive doses")
    log_scale = float(np.mean(np.log(measured_rad) - np.log(response)))
    scale = float(np.exp(log_scale))
    predicted = scale * response
    residual = np.log(predicted / measured_rad)
    return scale, predicted, float(np.sum(residual**2))


def gaussian_model_grid(
    energies: np.ndarray,
    kernel: np.ndarray,
    channels: list[int],
    mu_range: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    fine_energy = np.arange(float(energies.min()), float(energies.max()) + 0.0125, 0.025)
    fine_kernel = interpolate_kernel(energies, kernel[:, channels], fine_energy)
    mus = np.arange(mu_range[0], mu_range[1] + 0.0125, 0.025)
    sigmas = np.arange(0.05, 4.0001, 0.05)
    models = np.empty((len(mus), len(sigmas), len(channels)), dtype=float)
    for i, mu in enumerate(mus):
        z = (fine_energy[None, :] - mu) / sigmas[:, None]
        weight = np.exp(-0.5 * z**2)
        weight /= weight.sum(axis=1, keepdims=True)
        models[i] = weight @ fine_kernel
    return mus, sigmas, models


def centered_log_objective(measured_rad: np.ndarray, models: np.ndarray) -> np.ndarray:
    observed = np.log(measured_rad)
    observed -= observed.mean()
    model_log = np.log(models)
    model_log -= model_log.mean(axis=-1, keepdims=True)
    return np.sum((model_log - observed) ** 2, axis=-1)


def fit_phase(
    label: str,
    dose_path: Path,
    energies: np.ndarray,
    kernel: np.ndarray,
    bootstrap_draws: int,
) -> tuple[dict, pd.DataFrame, dict[str, np.ndarray]]:
    cfg = PHASES[label]
    table = pd.read_csv(dose_path).set_index("channel")
    channels = cfg["channels"]
    measured_krad = table.loc[channels, cfg["dose_column"]].to_numpy(float)
    measured_rad = measured_krad * 1000.0

    mono_response = kernel[:, channels]
    mono_obj = centered_log_objective(measured_rad, mono_response)
    mono_index = int(np.argmin(mono_obj))
    mono_energy = float(energies[mono_index])
    mono_scale, mono_predicted, mono_rss = fit_scale_log(measured_rad, mono_response[mono_index])

    mus, sigmas, gaussian_models = gaussian_model_grid(
        energies, kernel, channels, cfg["mu_range"]
    )
    gaussian_obj = centered_log_objective(measured_rad, gaussian_models)
    best_index = np.unravel_index(int(np.argmin(gaussian_obj)), gaussian_obj.shape)
    best_mu = float(mus[best_index[0]])
    best_sigma = float(sigmas[best_index[1]])
    best_response = gaussian_models[best_index]
    gaussian_scale, gaussian_predicted, gaussian_rss = fit_scale_log(measured_rad, best_response)

    n = len(channels)
    mono_bic = n * np.log(max(mono_rss, 1e-15) / n) + 2 * np.log(n)
    gaussian_bic = n * np.log(max(gaussian_rss, 1e-15) / n) + 3 * np.log(n)

    loo = []
    for omitted in range(n):
        keep = np.arange(n) != omitted
        obj = centered_log_objective(measured_rad[keep], gaussian_models[..., keep])
        idx = np.unravel_index(int(np.argmin(obj)), obj.shape)
        loo.append((float(mus[idx[0]]), float(sigmas[idx[1]]), channels[omitted]))

    rng = np.random.default_rng(20260818 + int(label[:2]))
    boot_mu = np.empty(bootstrap_draws)
    boot_sigma = np.empty(bootstrap_draws)
    # The README reports about 11% channel-to-channel scatter for the trusted
    # VT01 shielding comparison. Treat it as independent log-scale uncertainty
    # for a sensitivity interval; it is not a formal calibration uncertainty.
    for draw in range(bootstrap_draws):
        perturbed = measured_rad * np.exp(rng.normal(0.0, 0.11, size=n))
        obj = centered_log_objective(perturbed, gaussian_models)
        idx = np.unravel_index(int(np.argmin(obj)), obj.shape)
        boot_mu[draw] = mus[idx[0]]
        boot_sigma[draw] = sigmas[idx[1]]

    rows = []
    for i, channel in enumerate(channels):
        rows.append(
            {
                "phase": label,
                "channel": channel,
                "group": table.loc[channel, "group"],
                "sample": table.loc[channel, "sample"],
                "measured_dose_krad": measured_krad[i],
                "mono_predicted_krad": mono_predicted[i] / 1000.0,
                "mono_residual_percent": 100.0 * (mono_predicted[i] / measured_rad[i] - 1.0),
                "gaussian_predicted_krad": gaussian_predicted[i] / 1000.0,
                "gaussian_residual_percent": 100.0
                * (gaussian_predicted[i] / measured_rad[i] - 1.0),
            }
        )

    summary = {
        "phase": label,
        "trusted_channels": channels,
        "nominal_energy_MeV": cfg["nominal_energy"],
        "log_file_fluence_protons_cm2": cfg["log_fluence"],
        "mono": {
            "energy_MeV": mono_energy,
            "fluence_protons_cm2": mono_scale,
            "fluence_ratio_to_log": mono_scale / cfg["log_fluence"],
            "geometric_rms_fractional_residual": float(
                np.sqrt(np.mean(np.log(mono_predicted / measured_rad) ** 2))
            ),
            "rss_log": mono_rss,
            "bic": float(mono_bic),
        },
        "gaussian": {
            "mean_energy_MeV": best_mu,
            "sigma_energy_MeV": best_sigma,
            "fwhm_energy_MeV": 2.354820045 * best_sigma,
            "fluence_protons_cm2": gaussian_scale,
            "fluence_ratio_to_log": gaussian_scale / cfg["log_fluence"],
            "geometric_rms_fractional_residual": float(
                np.sqrt(np.mean(np.log(gaussian_predicted / measured_rad) ** 2))
            ),
            "rss_log": gaussian_rss,
            "bic": float(gaussian_bic),
            "delta_bic_vs_mono": float(gaussian_bic - mono_bic),
            "sensitivity_interval_11pct_scatter": {
                "mean_energy_MeV_16_50_84": np.quantile(boot_mu, [0.16, 0.5, 0.84]).tolist(),
                "sigma_energy_MeV_16_50_84": np.quantile(boot_sigma, [0.16, 0.5, 0.84]).tolist(),
            },
            "leave_one_channel_out": [
                {"omitted_channel": ch, "mean_energy_MeV": mu, "sigma_energy_MeV": sigma}
                for mu, sigma, ch in loo
            ],
        },
    }
    arrays = {
        "mus": mus,
        "sigmas": sigmas,
        "objective": gaussian_obj,
        "boot_mu": boot_mu,
        "boot_sigma": boot_sigma,
        "channels": np.asarray(channels),
        "measured_krad": measured_krad,
        "mono_predicted_krad": mono_predicted / 1000.0,
        "gaussian_predicted_krad": gaussian_predicted / 1000.0,
    }
    return summary, pd.DataFrame(rows), arrays


def plot_results(results: dict[str, tuple[dict, pd.DataFrame, dict]], output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.2), constrained_layout=True)
    for column, (label, (summary, rows, arrays)) in enumerate(results.items()):
        ax = axes[0, column]
        x = np.arange(len(rows))
        ax.plot(x, arrays["measured_krad"], "o", label="RadFET inferred")
        ax.plot(x, arrays["mono_predicted_krad"], "s--", label="best monoenergy")
        ax.plot(x, arrays["gaussian_predicted_krad"], "^--", label="best Gaussian")
        ax.set_xticks(x, [str(v) for v in arrays["channels"]])
        ax.set_xlabel("RadFET channel")
        ax.set_ylabel("Phase dose (krad)")
        ax.set_title(f"{label.replace('mev', ' MeV phase')}: trusted VT01 channels")
        ax.grid(alpha=0.25)
        ax.legend()

        ax = axes[1, column]
        objective = arrays["objective"]
        delta = objective - objective.min()
        image = ax.pcolormesh(
            arrays["mus"], arrays["sigmas"], delta.T, shading="auto", cmap="viridis_r"
        )
        best = summary["gaussian"]
        ax.plot(best["mean_energy_MeV"], best["sigma_energy_MeV"], "r*", markersize=12)
        ax.set_xlabel("Gaussian mean energy (MeV)")
        ax.set_ylabel("Gaussian sigma (MeV)")
        ax.set_title("Relative-dose pattern misfit")
        fig.colorbar(image, ax=ax, label="Delta log-RSS")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kernel", type=Path, required=True)
    parser.add_argument("--dose-64", type=Path, required=True)
    parser.add_argument("--dose-85", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=500)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    energies, kernel, _kernel_error = load_kernel(args.kernel)
    results = {}
    all_rows = []
    for label, dose_path in (("64mev", args.dose_64), ("85mev", args.dose_85)):
        result = fit_phase(label, dose_path, energies, kernel, args.bootstrap)
        results[label] = result
        all_rows.append(result[1])

    summary = {label: result[0] for label, result in results.items()}
    (args.output_dir / "spectrum-fit-summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    pd.concat(all_rows, ignore_index=True).to_csv(
        args.output_dir / "spectrum-fit-channel-residuals.csv", index=False
    )
    plot_results(results, args.output_dir / "spectrum-fit-diagnostics.png")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()


