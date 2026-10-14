"""Does the refitted centre energy improve agreement with the measured doses?

The -NewGaussian- runs repeat the -Gaussian- NewRadFET-Cuts productions with the
centre energy refitted from the EnergyKernel-NewRadFET response kernel:
62.07 MeV instead of nominal 64, and 83.47 MeV instead of nominal 85. Geometry,
detector model, production cut, beam macro, normalisation and history count are
identical, so a per-channel ratio between the two isolates the spectrum change.

The test that matters is against the measurement. For each simulated spectrum,
the per-channel prediction is

    D_ch = C * Phi_ch * K_ch

with K_ch the simulated dose coefficient, Phi_ch the delivered fluence at that
channel's location from the accelerator logs, and C one calibration constant
shared by every channel. If the refitted energy is a better description of the
beam, the RMS fractional residual must fall.

Aggregation is Dependencies.TotalDose.totalDose via AnalyseNewRadFET.aggregate;
nothing here re-parses a GRAS file.

ENVIRONMENT NOTE. This runs on the system python3, which carries numpy,
matplotlib and uncertainties but not pandas. scicomp-python-env is the other way
round, so InferSpectrum cannot be imported alongside Dependencies.TotalDose.
fit_scale_log below is therefore a transcription of
RadExUppsala/InferSpectrum.fit_scale_log rather
than an import of it; keep the two in step if either changes. It is scale
invariant, so working in kRad here and in rad there makes no difference.
"""

import argparse
import csv
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PYTHON_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_ROOT))

from RadExUppsala import AnalyseNewRadFET as anr  # noqa: E402

BASE_PATH = anr.BASE_PATH

N_CHANNELS = anr.N_CHANNELS
INPUTS = BASE_PATH / "EnergyKernel-FinalGeometry" / "SpectrumInference" / "Inputs"
LOGS = Path("/scratch/work/fetzera1/Radiation Testing/2025-06 Proton/"
            "Uppsala Accelerator Logs/results/final/aggregates")

CONFIG = {
    64: {
        "nominal": "64MeVProton-Gaussian-MeasuredGeometry-NewRadFET-Cuts",
        "refit": "64MeVProton-NewGaussian-MeasuredGeometry-NewRadFET-Cuts",
        "nominal_job": 20313480, "refit_job": 20328954,
        "nominal_spectrum": (64.00, 0.50429), "refit_spectrum": (62.07, 0.49269),
        "dose_csv": INPUTS / "phase-end-dose-64mev.csv",
        "dose_column": "dose_krad",
        "fluence_csv": LOGS / "RadEx_64MeV" / "radfet_location_fluence.csv",
        "channel_sets": {"ch0-6": [0, 1, 2, 3, 4, 5, 6]},
    },
    85: {
        "nominal": "85MeVProton-Gaussian-MeasuredGeometry-NewRadFET-Cuts",
        "refit": "85MeVProton-NewGaussian-MeasuredGeometry-NewRadFET-Cuts",
        "nominal_job": 20309279, "refit_job": 20328955,
        "nominal_spectrum": (85.00, 0.61600), "refit_spectrum": (83.47, 0.60875),
        "dose_csv": INPUTS / "phase-end-dose-85mev.csv",
        "dose_column": "dose_krad_85mev_increment",
        "fluence_csv": LOGS / "Radex_85MeV" / "radfet_location_fluence.csv",
        # The refitted 85 MeV centre came from the no-Ch10 fit, so both sets are
        # reported: Ch10 is the back-mounted channel and the known outlier.
        "channel_sets": {"ch0-6 + ch10": [0, 1, 2, 3, 4, 5, 6, 10],
                         "ch0-6 (no Ch10)": [0, 1, 2, 3, 4, 5, 6]},
    },
}


def read_column(path, key_column, value_columns):
    """Minimal CSV reader; pandas is not available in this interpreter."""
    table = {}
    with Path(path).open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            table[int(row[key_column])] = {c: row[c] for c in value_columns}
    if not table:
        raise RuntimeError(f"{path}: no rows")
    return table


def fit_scale_log(measured, response):
    """Transcription of InferSpectrum.fit_scale_log. See the module docstring."""
    if np.any(measured <= 0) or np.any(response <= 0):
        raise ValueError("Log-pattern fit requires strictly positive doses")
    scale = float(np.exp(np.mean(np.log(measured) - np.log(response))))
    predicted = scale * response
    residual = np.log(predicted / measured)
    return scale, predicted, float(np.sum(residual**2))


def fit(measured, coefficients, fluence, channels):
    observed = np.array([measured[c] for c in channels], dtype=float)
    model = np.array([coefficients[c] * fluence[c] for c in channels], dtype=float)
    scale, predicted, rss = fit_scale_log(observed, model)
    return {"scale": scale,
            "rms_percent": 100.0 * float(np.sqrt(rss / len(channels))),
            "residual_percent": 100.0 * (predicted / observed - 1.0),
            "channels": channels}


def analyse(energy, rows):
    config = CONFIG[energy]
    dose_table = read_column(config["dose_csv"], "channel",
                             [config["dose_column"], "sample"])
    measured = {c: float(v[config["dose_column"]]) for c, v in dose_table.items()}
    samples = {c: v["sample"] for c, v in dose_table.items()}
    fluence_table = read_column(config["fluence_csv"], "channel",
                                ["delivered_fluence_protons_cm2"])
    fluence = {c: float(v["delivered_fluence_protons_cm2"])
               for c, v in fluence_table.items()}
    if any(value <= 0 for value in fluence.values()):
        raise RuntimeError(f"{config['fluence_csv']}: non-positive fluence")

    print(f"\n{'=' * 74}\n{energy} MeV phase\n{'=' * 74}")
    print(f"  nominal  {config['nominal_spectrum'][0]:.2f} MeV "
          f"sigma {config['nominal_spectrum'][1]:.5f}  job {config['nominal_job']}")
    print(f"  refitted {config['refit_spectrum'][0]:.2f} MeV "
          f"sigma {config['refit_spectrum'][1]:.5f}  job {config['refit_job']}")

    runs = {}
    for key in ("nominal", "refit"):
        results, dose, error = anr.aggregate(config[key], 2 * N_CHANNELS)
        runs[key] = {"oxide": dose[:N_CHANNELS], "oxide_error": error[:N_CHANNELS],
                     "die": dose[N_CHANNELS:], "die_error": error[N_CHANNELS:],
                     "results": results}

    full_dose = np.concatenate([runs["refit"]["oxide"], runs["refit"]["die"]])
    full_error = np.concatenate([runs["refit"]["oxide_error"],
                                 runs["refit"]["die_error"]])
    anr.write_total_dose_csv(config["refit"], full_dose, full_error,
                             runs["refit"]["results"]["non-zeros"])
    anr.write_diagnostic_plots(config["refit"], full_dose, full_error,
                               runs["refit"]["results"]["non-zeros"])

    ratios = {}
    for layer in ("oxide", "die"):
        a, ae = runs["nominal"][layer], runs["nominal"][f"{layer}_error"]
        b, be = runs["refit"][layer], runs["refit"][f"{layer}_error"]
        ratio = b / a
        ratios[layer] = (ratio, ratio * np.hypot(ae / a, be / b))

    print("\nSpectrum effect, refitted / nominal [%]")
    print(f"{'ch':>3s} {'sample':<24s}{'oxide':>9s}{'+-':>7s}{'die':>9s}{'+-':>7s}")
    for ch in range(N_CHANNELS):
        mount = " (back)" if ch in anr.BACK_MOUNTED else ""
        print(f"{ch:>3d} {samples.get(ch, '?') + mount:<24s}"
              f"{100 * (ratios['oxide'][0][ch] - 1):>+8.2f}"
              f"{100 * ratios['oxide'][1][ch]:>7.2f}"
              f"{100 * (ratios['die'][0][ch] - 1):>+8.2f}"
              f"{100 * ratios['die'][1][ch]:>7.2f}")

    for set_name, channels in config["channel_sets"].items():
        print(f"\nAgreement with measured doses, {set_name}")
        print(f"{'':22s}{'RMS residual':>14s}{'calibration C':>16s}")
        fits = {}
        for layer in ("oxide", "die"):
            for key in ("nominal", "refit"):
                fits[(layer, key)] = fit(measured, runs[key][layer], fluence,
                                         channels)
                print(f"  {layer + ', ' + key:<20s}"
                      f"{fits[(layer, key)]['rms_percent']:>12.2f} %"
                      f"{fits[(layer, key)]['scale']:>16.4f}")
            change = (fits[(layer, "refit")]["rms_percent"]
                      - fits[(layer, "nominal")]["rms_percent"])
            verdict = "better" if change < 0 else "worse"
            print(f"  {'-> change':<20s}{change:>+12.2f} pp   ({verdict})")

        print(f"\n  Per-channel residuals [%]")
        for layer in ("oxide", "die"):
            for key in ("nominal", "refit"):
                result = fits[(layer, key)]
                pairs = "  ".join(
                    f"ch{c}:{v:+.1f}"
                    for c, v in zip(result["channels"], result["residual_percent"]))
                print(f"    {layer:<6s}{key:<9s}{pairs}")

        for index, channel in enumerate(channels):
            rows.append({
                "energy_MeV": energy, "channel_set": set_name, "channel": channel,
                "sample": samples.get(channel, ""),
                "measured_dose_krad": f"{measured[channel]:.9e}",
                "delivered_fluence_protons_cm2": f"{fluence[channel]:.9e}",
                **{f"{layer}_{key}_coefficient": f"{runs[key][layer][channel]:.9e}"
                   for layer in ("oxide", "die") for key in ("nominal", "refit")},
                **{f"{layer}_{key}_residual_percent":
                   f"{fits[(layer, key)]['residual_percent'][index]:+.6f}"
                   for layer in ("oxide", "die") for key in ("nominal", "refit")},
                **{f"{layer}_{key}_rms_percent":
                   f"{fits[(layer, key)]['rms_percent']:.6f}"
                   for layer in ("oxide", "die") for key in ("nominal", "refit")},
            })
    return runs, config


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--energies", type=int, nargs="+", default=[64, 85])
    parser.add_argument("--output", type=Path,
                        default=BASE_PATH / "RadEx-Uppsala-NewGaussian-comparison.csv")
    parser.add_argument("--plot", type=Path,
                        default=BASE_PATH / "NewGaussian_vs_Gaussian.pdf")
    args = parser.parse_args()

    rows = []
    panels = {}
    for energy in args.energies:
        panels[energy] = analyse(energy, rows)

    with args.output.open("w", newline="", encoding="utf-8") as handle:
        lead = ["energy_MeV", "channel_set", "channel", "sample",
                "measured_dose_krad", "delivered_fluence_protons_cm2"]
        rest = sorted({k for row in rows for k in row} - set(lead))
        writer = csv.DictWriter(handle, fieldnames=lead + rest)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {args.output}")

    figure, axes = plt.subplots(len(panels), 1, figsize=(8, 4 * len(panels)),
                                squeeze=False)
    for axis, (energy, (runs, config)) in zip(axes[:, 0], sorted(panels.items())):
        channels = np.arange(N_CHANNELS)
        for layer, offset, marker in (("oxide", -0.08, "o"), ("die", 0.08, "s")):
            a, ae = runs["nominal"][layer], runs["nominal"][f"{layer}_error"]
            b, be = runs["refit"][layer], runs["refit"][f"{layer}_error"]
            ratio = b / a
            axis.errorbar(channels + offset, 100.0 * (ratio - 1.0),
                          yerr=100.0 * ratio * np.hypot(ae / a, be / b),
                          fmt=marker, capsize=4, elinewidth=1, capthick=1,
                          label=layer)
        axis.axhline(0.0, color="r", linestyle="--", label="no change")
        axis.axvline(9.5, color="grey", linestyle=":")
        axis.set_title(f"{energy} MeV: {config['refit_spectrum'][0]:.2f} MeV refitted "
                       f"vs {config['nominal_spectrum'][0]:.0f} MeV nominal",
                       fontsize=10)
        axis.set_xlabel("Channel")
        axis.set_ylabel("Change in dose [%]")
        axis.set_xticks(channels)
        axis.grid(which="both", alpha=0.3)
        axis.legend(fontsize=8)
    figure.savefig(args.plot, format="pdf", bbox_inches="tight")
    print(f"Wrote {args.plot}")


if __name__ == "__main__":
    main()
