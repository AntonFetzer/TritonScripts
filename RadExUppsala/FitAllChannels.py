"""Published-width energy fits to all twelve RadFETs, and without Ch6.

Uses the dose model of InferSpectrum (uncertainty weights, each phase's
per-channel fluence) with the beam width fixed to the published Skandion
relation, InferSpectrum.published_sigma. For each centre energy the calibration
scale is refitted and the energy minimises chi2 on a 0.025 MeV grid; the
interval is delta-chi2 <= 1, and also <= chi2/dof when the fit is worse than the
assumed uncertainties.

The default calibration model is one scale shared by all channels; the per-type
model (separate VT01 and VT05 scales) is reported alongside as a diagnostic.
Two channel sets are fitted per phase: all twelve, and all but Ch6, whose 64 MeV
dose is about 30 % below the model under its 2 mm absorber. --probe PHASE:MeV
prints every channel's prediction at a given centre energy for both models.

Run from the repository root:

    python3 -m RadExUppsala.FitAllChannels --kernel <energy-kernel-dose-summary-oxide.csv> \\
        --dose-64 <phase-end-dose-64mev.csv> --dose-85 <phase-end-dose-85mev.csv> \\
        --fluence-64 <RadEx_64MeV/radfet_location_fluence.csv> \\
        --fluence-85 <Radex_85MeV/radfet_location_fluence.csv> --output <fit.json>

writes <fit.json> and the figure <fit.pdf> (doses, pulls and chi2 curves).
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

PYTHON_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_ROOT))

from RadExUppsala import InferSpectrum as inf  # noqa: E402

CHANNEL_SETS = (("all 12", inf.CHANNELS), ("without Ch6", [c for c in inf.CHANNELS if c != 6]))
MODEL_COLOUR = {"shared": "C0", "per-type": "grey"}
SET_STYLE = {"all 12": "-", "without Ch6": ":"}


def fit(kernel_path, dose, sigma, fluence_path, phase, channels, model, mu_range):
    lo, hi = mu_range
    mus = np.arange(lo, hi + inf.FINE_STEP / 2, inf.FINE_STEP)
    energies, kernel, error = inf.weighted_kernel(kernel_path, fluence_path, phase, channels)
    means, errs = inf.fold(energies, kernel, error, mus)
    sel = np.array(channels)
    d, sd = dose[sel], sigma[sel]
    groups = inf.scale_groups(channels, model)
    chi2 = np.array([inf.solve_scales(d, sd, means[i], errs[i], groups)[1] for i in range(len(mus))])
    best = int(np.argmin(chi2))
    scales, chi2_min, pulls, scale = inf.solve_scales(d, sd, means[best], errs[best], groups)
    dof = len(channels) - 1 - len(groups)
    per_dof = chi2_min / dof if dof > 0 else float("nan")
    inside = mus[chi2 <= chi2_min + 1.0]
    inside_scaled = mus[chi2 <= chi2_min + max(1.0, per_dof)]
    predicted = scale * means[best]
    log_fluence = inf.PHASES[phase]["log_fluence"]
    return {
        "phase": phase, "model": model, "channels": channels, "mu_range_MeV": [lo, hi],
        "mu_MeV": float(mus[best]), "mu_lo": float(inside.min()), "mu_hi": float(inside.max()),
        "mu_lo_scaled": float(inside_scaled.min()), "mu_hi_scaled": float(inside_scaled.max()),
        "at_edge": bool(best in (0, len(mus) - 1)),
        "chi2": chi2_min, "dof": dof, "chi2_per_dof": per_dof,
        "calibration_ratio": {k: v / log_fluence for k, v in scales.items()},
        "pull": dict(zip(map(int, channels), map(float, pulls))),
        "residual_percent": dict(zip(map(int, channels), map(float, 100.0 * (predicted / d - 1.0)))),
        "measured_rad": dict(zip(map(int, channels), map(float, d))),
        "sigma_rad": dict(zip(map(int, channels), map(float, sd))),
        # Arrays for the figure only; not written to the JSON.
        "_mus": mus, "_chi2": chi2, "_predicted_rad": predicted,
    }


def probe(kernel_path, dose, sigma, fluence_path, phase, energy):
    """Every channel's prediction at one centre energy, published width, both scale models."""
    channels = inf.CHANNELS
    energies, kernel, error = inf.weighted_kernel(kernel_path, fluence_path, phase, channels)
    means, errs = inf.fold(energies, kernel, error, np.array([energy]))
    out = {}
    for model in inf.SCALE_MODELS:
        _, chi2, pulls, scale = inf.solve_scales(dose[channels], sigma[channels], means[0], errs[0],
                                                 inf.scale_groups(channels, model))
        out[model] = (chi2, pulls, scale * means[0])
    print(f"   probe {phase} at {energy:.2f} MeV: chi2 shared {out['shared'][0]:.1f}, "
          f"per-type {out['per-type'][0]:.1f}")
    for c in channels:
        print(f"      ch{c:<2d} {inf.DEVICE[c]}  measured {dose[c] / 1000:10.4g} krad  "
              f"shared {out['shared'][2][c] / 1000:10.4g} (pull {out['shared'][1][c]:+6.1f})  "
              f"per-type {out['per-type'][2][c] / 1000:10.4g} (pull {out['per-type'][1][c]:+6.1f})")


def plot(results, output):
    """Doses, pulls and chi2 curves, one column per phase; shared scale first, per-type grey."""
    phases = list(dict.fromkeys(r["phase"] for r in results))
    fig, axes = plt.subplots(3, len(phases), figsize=(6.2 * len(phases), 11.5), squeeze=False)
    for column, phase in enumerate(phases):
        rows = {(r["set"], r["model"]): r for r in results if r["phase"] == phase}
        full = rows[("all 12", "shared")]
        x = np.array(inf.CHANNELS)
        measured = np.array([full["measured_rad"][c] for c in inf.CHANNELS]) / 1000.0
        sigma = np.array([full["sigma_rad"][c] for c in inf.CHANNELS]) / 1000.0

        ax = axes[0, column]
        for device in ("VT01", "VT05"):
            sel = np.array([inf.DEVICE[c] == device for c in inf.CHANNELS])
            lower = np.minimum(sigma[sel], 0.95 * measured[sel])
            ax.errorbar(x[sel], measured[sel], yerr=[lower, sigma[sel]], fmt=inf.MARKER[device],
                        color="black", capsize=3, label=f"measured {device}")
        for model in inf.SCALE_MODELS:
            r = rows[("all 12", model)]
            pred = r["_predicted_rad"] / 1000.0
            name = "shared scale" if model == "shared" else "diagnostic, per-type scales"
            ax.plot(x, np.where(pred > 0, pred, np.nan), "_", color=MODEL_COLOUR[model], ms=16, mew=2,
                    label=f"{name}, {r['mu_MeV']:.2f} MeV")
        ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xlabel("RadFET channel")
        ax.set_ylabel("Phase dose (krad)")
        ax.set_title(f"{phase.replace('mev', ' MeV phase')}: all 12 channels, published width")
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=8)

        ax = axes[1, column]
        width = 0.38
        for k, model in enumerate(inf.SCALE_MODELS):
            r = rows[("all 12", model)]
            ax.bar(x + (k - 0.5) * width, [r["pull"][c] for c in inf.CHANNELS], width,
                   color=MODEL_COLOUR[model], label=f"{model}: chi2 {r['chi2']:.1f}/{r['dof']}")
        for level in (-2, 2):
            ax.axhline(level, color="grey", ls="--", lw=0.8)
        ax.axhline(0, color="black", lw=0.8)
        ax.set_xticks(x)
        ax.set_xlabel("RadFET channel")
        ax.set_ylabel("Pull (predicted - measured) / sigma")
        ax.grid(alpha=0.25, axis="y")
        ax.legend(fontsize=8)

        ax = axes[2, column]
        for (label, model), r in rows.items():
            ax.plot(r["_mus"], r["_chi2"] - r["chi2"], SET_STYLE[label], color=MODEL_COLOUR[model],
                    label=f"{label}, {model}: {r['mu_MeV']:.2f} MeV")
        ax.axhline(1.0, color="grey", lw=0.8)
        ax.set_ylim(0, 25)
        ax.set_xlabel("Centre energy (MeV)")
        ax.set_ylabel("Delta chi2 (each curve to its own minimum)")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--kernel", type=Path, required=True)
    parser.add_argument("--dose-64", type=Path, required=True)
    parser.add_argument("--dose-85", type=Path, required=True)
    parser.add_argument("--fluence-64", type=Path, required=True)
    parser.add_argument("--fluence-85", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mu-range-64", type=float, nargs=2, default=inf.PHASES["64mev"]["mu_range"])
    parser.add_argument("--mu-range-85", type=float, nargs=2, default=inf.PHASES["85mev"]["mu_range"])
    parser.add_argument("--probe", action="append", default=[],
                        help="PHASE:MeV, e.g. 85mev:84.0; repeatable")
    args = parser.parse_args()
    probes = [(s.split(":")[0], float(s.split(":")[1])) for s in args.probe]

    results = []
    for phase, fluence_path in (("64mev", args.fluence_64), ("85mev", args.fluence_85)):
        dose, sigma, noise = inf.measurements(args.dose_64, args.dose_85, phase)
        mu_range = tuple(args.mu_range_64 if phase == "64mev" else args.mu_range_85)
        print(f"== {phase}: measured dose [rad] +- total (noise part); "
              f"centre energy searched over {mu_range[0]:g}-{mu_range[1]:g} MeV")
        print("   " + "  ".join(f"ch{c}:{dose[c]:.4g}+-{sigma[c]:.2g}({noise[c]:.1g})" for c in inf.CHANNELS))
        for label, channels in CHANNEL_SETS:
            for model in inf.SCALE_MODELS:
                r = fit(args.kernel, dose, sigma, fluence_path, phase, channels, model, mu_range)
                r["set"] = label
                results.append(r)
                cal = ", ".join(f"{k} {v:.3f}" for k, v in r["calibration_ratio"].items())
                edge = "  AT GRID EDGE" if r["at_edge"] else ""
                print(f"   {label:<12s}{model:<9s} mu {r['mu_MeV']:6.2f} MeV  "
                      f"[{r['mu_lo']:.2f}-{r['mu_hi']:.2f}; scaled {r['mu_lo_scaled']:.2f}-{r['mu_hi_scaled']:.2f}]  "
                      f"chi2/dof {r['chi2']:.1f}/{r['dof']}  C: {cal}{edge}")
                print("      pulls " + " ".join(f"ch{c}:{p:+.1f}" for c, p in r["pull"].items()))
        for probe_phase, energy in probes:
            if probe_phase == phase:
                probe(args.kernel, dose, sigma, fluence_path, phase, energy)
    public = [{k: v for k, v in r.items() if not k.startswith("_")} for r in results]
    args.output.write_text(json.dumps(public, indent=2), encoding="utf-8")
    figure = args.output.with_suffix(".pdf")
    plot(results, figure)
    print(f"Wrote {args.output} and {figure}")


if __name__ == "__main__":
    main()
