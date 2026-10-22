"""Published-width energy fits to all RadFETs of a campaign.

Uses the measurement model of InferSpectrum (fit space, uncertainty weights,
each phase's per-channel fluence) with the beam width fixed to the published
Skandion relation, InferSpectrum.published_sigma. For each centre energy the
calibration scales (dose space) or response amplitudes (shift space) are
refitted and the energy minimises chi2 on a 0.025 MeV grid; the interval is
delta-chi2 <= 1, and also <= chi2/dof when the fit is worse than the assumed
uncertainties.

Phases are fitted in order. In shift space a phase on pre-dosed devices starts
from the previous phase's unscaled model dose at the energy its default model
fitted, so its measured quantity is the total shift at its end.

Each phase uses the default calibration-scale model its campaign file names and
reports the phase's diagnostic models alongside. --probe PHASE:MeV prints every
channel's prediction at a given centre energy for every model.
--scan-exponent TYPE profiles chi2 of each phase's default model over the
shift-space response exponent of one device type (0.30-1.00 in steps of 0.01),
with the centre energy refitted at each exponent, and prints the best exponent
and its delta-chi2 <= 1 interval; it only makes sense in shift space.

Run from the repository root:

    python3 -m RadEx.FitAllChannels --kernel <energy-kernel-dose-summary-oxide.csv> \\
        --campaign <campaign.json> --output <fit.json> [--space dose|shift]

writes <fit.json> and the figure <fit.pdf> (measurements, pulls and chi2 curves).
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

from RadEx import InferSpectrum as inf  # noqa: E402

#: Colours by role: the phase's default model, then its diagnostics in order.
ROLE_COLOURS = ("C0", "grey", "magenta", "green")
#: Exponent grid of --scan-exponent.
EXPONENTS = np.round(np.arange(0.30, 1.0001, 0.01), 2)


def folded(kernel_path, campaign, phase, mus):
    channels = campaign["channels"]
    energies, kernel, error = inf.weighted_kernel(kernel_path, campaign, phase, channels)
    return inf.fold(energies, kernel, error, mus)


def fit(meas, means, errs, mus, campaign, phase, model):
    channels = campaign["channels"]
    groups = inf.scale_groups(campaign, channels, model)
    chi2 = np.array([inf.solve(meas, means[i], errs[i], groups)[1] for i in range(len(mus))])
    best = int(np.argmin(chi2))
    scales, chi2_min, pulls, predicted = inf.solve(meas, means[best], errs[best], groups)
    dof = len(channels) - 1 - len(groups)
    per_dof = chi2_min / dof if dof > 0 else float("nan")
    inside = mus[chi2 <= chi2_min + 1.0]
    inside_scaled = mus[chi2 <= chi2_min + max(1.0, per_dof)]
    log_fluence = campaign["phases"][phase]["log_fluence"]
    value = meas["value"]
    sigma = meas["sigma"] if meas["space"] == "dose" else meas["noise"]
    with np.errstate(divide="ignore", invalid="ignore"):
        residual = np.where(value != 0, 100.0 * (predicted / value - 1.0), np.nan)
    out = {
        "phase": phase, "model": model, "fit_space": meas["space"], "unit": meas["unit"],
        "systematic_fraction": campaign["systematic"], "noise_floor_mv": 1e3 * campaign["noise_floor"],
        "channels": channels, "mu_range_MeV": [float(mus[0]), float(mus[-1])],
        "mu_MeV": float(mus[best]), "mu_lo": float(inside.min()), "mu_hi": float(inside.max()),
        "mu_lo_scaled": float(inside_scaled.min()), "mu_hi_scaled": float(inside_scaled.max()),
        "at_edge": bool(best in (0, len(mus) - 1)),
        "chi2": chi2_min, "dof": dof, "chi2_per_dof": per_dof,
        **inf.scale_report(meas, scales, log_fluence),
        "pull": dict(zip(map(int, channels), map(float, pulls))),
        "residual_percent": dict(zip(map(int, channels), map(float, residual))),
        "measured": dict(zip(map(int, channels), map(float, value))),
        "sigma": dict(zip(map(int, channels), map(float, sigma))),
        "predicted": dict(zip(map(int, channels), map(float, predicted))),
        # Arrays for the figure and the next phase only; not written to the JSON.
        "_mus": mus, "_chi2": chi2, "_predicted": predicted, "_model_rad": means[best] * log_fluence,
    }
    if meas["space"] == "shift":
        out["exponent"] = dict(zip(map(int, channels), map(float, meas["exponent"])))
        out["offset_krad"] = dict(zip(map(int, channels), map(float, meas["offset"] / 1000.0)))
        out["previous_model_rad"] = dict(zip(map(int, channels), map(float, meas["previous"])))
    return out


def scales_text(r):
    key = "calibration_ratio" if "calibration_ratio" in r else "amplitude_v_per_rad_pow_n"
    return ", ".join(f"{k} {v:.4g}" for k, v in r[key].items())


def probe(meas, kernel_path, campaign, phase, energy):
    """Every channel's prediction at one centre energy, published width, every model of the phase."""
    channels = campaign["channels"]
    means, errs = folded(kernel_path, campaign, phase, np.array([energy]))
    out = {}
    for model in inf.phase_models(campaign, phase):
        _, chi2, pulls, predicted = inf.solve(meas, means[0], errs[0],
                                              inf.scale_groups(campaign, channels, model))
        out[model] = (chi2, pulls, predicted)
    scale = 1000.0 if meas["unit"] == "rad" else 1.0
    unit = "krad" if meas["unit"] == "rad" else "V"
    print(f"   probe {phase} at {energy:.2f} MeV: chi2 " + ", ".join(f"{m} {o[0]:.1f}" for m, o in out.items()))
    for i, c in enumerate(channels):
        print(f"      ch{c:<2d} {campaign['device'][c]}  measured {meas['value'][i] / scale:10.4g} {unit}  " +
              "  ".join(f"{m} {o[2][i] / scale:10.4g} (pull {o[1][i]:+6.1f})" for m, o in out.items()))


def scan_exponent(kernel_path, campaign, phase, device, previous_model, means, errs, mus):
    """chi2 of the phase's default model against the exponent of one device type."""
    model = inf.phase_models(campaign, phase)[0]
    saved = campaign["exponent"][device]
    chi2 = np.empty(len(EXPONENTS))
    energy = np.empty(len(EXPONENTS))
    for i, n in enumerate(EXPONENTS):
        campaign["exponent"][device] = float(n)
        meas = inf.measurement(campaign, phase, previous_model)
        r = fit(meas, means, errs, mus, campaign, phase, model)
        chi2[i], energy[i] = r["chi2"], r["mu_MeV"]
    campaign["exponent"][device] = saved
    best = int(np.argmin(chi2))
    inside = EXPONENTS[chi2 <= chi2[best] + 1.0]
    print(f"   exponent scan {device}, {model} model: best n {EXPONENTS[best]:.2f} "
          f"[{inside.min():.2f}-{inside.max():.2f}], chi2 {chi2[best]:.1f} at {energy[best]:.2f} MeV; "
          f"at n = {saved:.2f}: chi2 {chi2[np.argmin(np.abs(EXPONENTS - saved))]:.1f}")
    return {"phase": phase, "device": device, "model": model, "exponent": EXPONENTS.tolist(),
            "chi2": chi2.tolist(), "mu_MeV": energy.tolist(), "best_exponent": float(EXPONENTS[best]),
            "exponent_interval": [float(inside.min()), float(inside.max())]}


def plot(results, campaign, output):
    """Measurements, pulls and chi2 curves, one column per phase; default model first, diagnostics after."""
    channels = campaign["channels"]
    phases = list(dict.fromkeys(r["phase"] for r in results))
    fig, axes = plt.subplots(3, len(phases), figsize=(6.2 * len(phases), 11.5), squeeze=False)
    for column, phase in enumerate(phases):
        rows = {r["model"]: r for r in results if r["phase"] == phase}
        models = inf.phase_models(campaign, phase)
        colour = dict(zip(models, ROLE_COLOURS))
        full = rows[models[0]]
        scale = 1000.0 if full["unit"] == "rad" else 1.0
        unit = "krad" if full["unit"] == "rad" else "V"
        x = np.array(channels)
        measured = np.array([full["measured"][c] for c in channels]) / scale
        sigma = np.array([full["sigma"][c] for c in channels]) / scale

        ax = axes[0, column]
        for device in sorted(set(campaign["device"].values())):
            sel = np.array([campaign["device"][c] == device for c in channels])
            lower = np.minimum(sigma[sel], 0.95 * np.abs(measured[sel]))
            ax.errorbar(x[sel], measured[sel], yerr=[lower, sigma[sel]], fmt=inf.MARKER[device],
                        color="black", capsize=3, label=f"measured {device}")
        for model in models:
            r = rows[model]
            pred = r["_predicted"] / scale
            name = f"{model} scale" if model == models[0] else f"diagnostic, {model} scale"
            ax.plot(x, np.where(pred > 0, pred, np.nan), "_", color=colour[model], ms=16, mew=2,
                    label=f"{name}, {r['mu_MeV']:.2f} MeV")
        ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xlabel("RadFET channel")
        ax.set_ylabel("Phase dose (krad)" if unit == "krad" else "Phase-end threshold shift (V)")
        ax.set_title(f"{phase.replace('mev', ' MeV phase')}: all {len(channels)} channels, "
                     f"published width, {full['fit_space']} space")
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=8)

        ax = axes[1, column]
        width = 0.8 / len(models)
        for k, model in enumerate(models):
            r = rows[model]
            ax.bar(x + (k - (len(models) - 1) / 2) * width, [r["pull"][c] for c in channels], width,
                   color=colour[model], label=f"{model}: chi2 {r['chi2']:.1f}/{r['dof']}")
        for level in (-2, 2):
            ax.axhline(level, color="grey", ls="--", lw=0.8)
        ax.axhline(0, color="black", lw=0.8)
        ax.set_xticks(x)
        ax.set_xlabel("RadFET channel")
        ax.set_ylabel("Pull (predicted - measured) / sigma")
        ax.grid(alpha=0.25, axis="y")
        ax.legend(fontsize=8)

        ax = axes[2, column]
        for model in models:
            r = rows[model]
            ax.plot(r["_mus"], r["_chi2"] - r["chi2"], "-", color=colour[model],
                    label=f"{model}: {r['mu_MeV']:.2f} MeV")
        ax.axhline(1.0, color="grey", lw=0.8)
        ax.set_ylim(0, 25)
        ax.set_xlim(*full["mu_range_MeV"])
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
    parser.add_argument("--campaign", type=Path, required=True,
                        help="campaign JSON next to the measured doses")
    parser.add_argument("--output", type=Path, required=True)
    inf.add_mu_range_argument(parser)
    inf.add_space_argument(parser)
    parser.add_argument("--probe", action="append", default=[],
                        help="PHASE:MeV, e.g. 85mev:84.0; repeatable")
    parser.add_argument("--scan-exponent", action="append", default=[], metavar="TYPE",
                        help="profile chi2 over the shift-space exponent of a device type; repeatable")
    args = parser.parse_args()
    campaign = inf.load_campaign(args.campaign, args.space, args.systematic, args.noise_floor_mv)
    mu_ranges = inf.mu_range_overrides(campaign, args.mu_range)
    probes = [(s.split(":")[0], float(s.split(":")[1])) for s in args.probe]
    unknown = sorted({p for p, _ in probes} - set(campaign["phases"]))
    if unknown:
        parser.error(f"--probe: unknown phases {unknown}; the campaign has {list(campaign['phases'])}")
    unknown = sorted(set(args.scan_exponent) - set(campaign["exponent"]))
    if unknown:
        parser.error(f"--scan-exponent: unknown device types {unknown}; the campaign has "
                     f"{sorted(campaign['exponent'])}")
    if args.scan_exponent and campaign["space"] != "shift":
        parser.error("--scan-exponent needs the shift space")

    print(f"error model: device fraction {campaign['systematic']:.2f}, noise floor {1e3 * campaign['noise_floor']:.1f} mV")
    print(f"fit space: {campaign['space']}" + (
        "; response dV = A ((D + D0)^n - D0^n) with " + ", ".join(
            f"{t}: n {campaign['exponent'][t]:.3f}, D0 {campaign['offset_krad'][t]:.2f} krad"
            for t in sorted(campaign["exponent"]))
        if campaign["space"] == "shift" else ""))
    results = []
    scans = []
    model_by_phase = {}
    for phase, cfg in campaign["phases"].items():
        previous = model_by_phase.get(cfg["previous_phase"]) if cfg["previous_phase"] else None
        meas = inf.measurement(campaign, phase, previous)
        mu_range = mu_ranges[phase]
        mus = np.arange(mu_range[0], mu_range[1] + inf.FINE_STEP / 2, inf.FINE_STEP)
        means, errs = folded(args.kernel, campaign, phase, mus)
        scale = 1000.0 if meas["unit"] == "rad" else 1.0
        unit = "krad" if meas["unit"] == "rad" else "V"
        sigma = meas["sigma"] if meas["space"] == "dose" else meas["noise"]
        print(f"== {phase}: measured {unit} +- sigma" +
              (" (noise part)" if meas["space"] == "dose" else "") +
              f"; centre energy searched over {mu_range[0]:g}-{mu_range[1]:g} MeV")
        print("   " + "  ".join(
            f"ch{c}:{meas['value'][i] / scale:.4g}+-{sigma[i] / scale:.2g}" +
            (f"({meas['noise'][i] / scale:.1g})" if meas["space"] == "dose" else "")
            for i, c in enumerate(campaign["channels"])))
        if meas["space"] == "shift" and previous is not None:
            print("   starts from the previous phase's model dose [krad]: " +
                  "  ".join(f"ch{c}:{meas['previous'][i] / 1000:.3g}"
                            for i, c in enumerate(campaign["channels"])))
        for model in inf.phase_models(campaign, phase):
            r = fit(meas, means, errs, mus, campaign, phase, model)
            results.append(r)
            edge = "  AT GRID EDGE" if r["at_edge"] else ""
            print(f"   {model:<16s} mu {r['mu_MeV']:6.2f} MeV  "
                  f"[{r['mu_lo']:.2f}-{r['mu_hi']:.2f}; scaled {r['mu_lo_scaled']:.2f}-{r['mu_hi_scaled']:.2f}]  "
                  f"chi2/dof {r['chi2']:.1f}/{r['dof']}  scales: {scales_text(r)}{edge}")
            print("      pulls " + " ".join(f"ch{c}:{p:+.1f}" for c, p in r["pull"].items()))
        default = next(r for r in results if r["phase"] == phase)
        model_by_phase[phase] = default["_model_rad"] + (previous if previous is not None else 0.0)
        for probe_phase, energy in probes:
            if probe_phase == phase:
                probe(meas, args.kernel, campaign, phase, energy)
        for device in args.scan_exponent:
            scans.append(scan_exponent(args.kernel, campaign, phase, device, previous, means, errs, mus))
    public = [{k: v for k, v in r.items() if not k.startswith("_")} for r in results]
    if scans:
        public.append({"exponent_scans": scans})
    args.output.write_text(json.dumps(public, indent=2), encoding="utf-8")
    figure = args.output.with_suffix(".pdf")
    plot(results, campaign, figure)
    print(f"Wrote {args.output} and {figure}")


if __name__ == "__main__":
    main()
