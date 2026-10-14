"""Fit the DUT pose as a nuisance parameter alongside the beam energy.

radfet_location_fluence.csv samples the isocentre heatmap at the NOMINAL RadFET
coordinates: "GDML coordinates are assumed to be nominally centred and aligned
with the machine axes; no hypothetical DUT translation or rotation is applied."
The residual that survives per-channel fluence weighting is largest on the
lateral extremes and changes sign between them, which is what a small offset or
rotation of the instrument would produce.

This script samples the same heatmap at transformed positions

    x' = x cos(theta) - y sin(theta) + dx
    y' = x sin(theta) + y cos(theta) + dy

folds the resulting per-channel fluence into the kernel exactly as
ApplyChannelFluence.py does, and minimises the same centred-log objective that
InferSpectrum.py uses, over (dx, dy, theta) as well as energy.

Everything numerical is imported from InferSpectrum so the fits stay directly
comparable: the objective, the scale fit, the phase configuration and the BIC
convention are its own.

Two pose models are reported:
  per-phase  -- independent pose for the 64 and 85 MeV irradiations
  joint      -- one shared pose, separate energy and scale per phase

The joint model is the physically motivated one if the instrument was not
re-seated between the two phases, and it is far better conditioned: 15 measured
channels against 7 parameters, versus 8 against 5 for a single phase. Both are
reported because whether the DUT was re-seated is not recorded.
"""

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize

PYTHON_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_ROOT))

from RadExUppsala import InferSpectrum as inf  # noqa: E402

AGG = Path("/scratch/work/fetzera1/Radiation Testing/2025-06 Proton/"
           "Uppsala Accelerator Logs/results/final/aggregates")
PHASE_DIR = {"64mev": AGG / "RadEx_64MeV", "85mev": AGG / "Radex_85MeV"}
PHASE_REFERENCE = {"64mev": 6.942e10, "85mev": 9.152e10}

# Nominal RadFET positions in GDML/machine mm, from radfet_location_fluence.csv.
CHANNEL_XY = {}


class Heatmap:
    """Bilinear sampler on the isocentre fluence grid."""

    def __init__(self, path):
        data = np.load(path)
        self.f = data["fluence_protons_cm2"].astype(float)   # indexed [y, x]
        self.x = 0.5 * (data["x_edges_mm"][:-1] + data["x_edges_mm"][1:])
        self.y = 0.5 * (data["y_edges_mm"][:-1] + data["y_edges_mm"][1:])

    def __call__(self, x, y):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        y = np.atleast_1d(np.asarray(y, dtype=float))
        ix = np.clip(np.searchsorted(self.x, x) - 1, 0, len(self.x) - 2)
        iy = np.clip(np.searchsorted(self.y, y) - 1, 0, len(self.y) - 2)
        tx = (x - self.x[ix]) / (self.x[ix + 1] - self.x[ix])
        ty = (y - self.y[iy]) / (self.y[iy + 1] - self.y[iy])
        return ((1 - tx) * (1 - ty) * self.f[iy, ix]
                + tx * (1 - ty) * self.f[iy, ix + 1]
                + (1 - tx) * ty * self.f[iy + 1, ix]
                + tx * ty * self.f[iy + 1, ix + 1])


def load_channel_xy(phase):
    rows = list(csv.DictReader((PHASE_DIR[phase] / "radfet_location_fluence.csv").open()))
    xy = {int(r["channel"]): (float(r["gdml_x_mm"]), float(r["gdml_y_mm"])) for r in rows}
    published = {int(r["channel"]): float(r["delivered_fluence_protons_cm2"]) for r in rows}
    return xy, published


def transform(xy, dx, dy, theta_deg):
    t = np.radians(theta_deg)
    x = np.array([p[0] for p in xy])
    y = np.array([p[1] for p in xy])
    return x * np.cos(t) - y * np.sin(t) + dx, x * np.sin(t) + y * np.cos(t) + dy


class Phase:
    """Everything needed to score one irradiation phase at a given pose."""

    def __init__(self, label, kernel_path, dose_path):
        self.label = label
        self.cfg = inf.PHASES[label]
        self.channels = self.cfg["channels"]
        self.reference = PHASE_REFERENCE[label]
        self.energies, self.kernel, _ = inf.load_kernel(kernel_path)

        import pandas as pd
        table = pd.read_csv(dose_path).set_index("channel")
        self.measured_rad = table.loc[self.channels, self.cfg["dose_column"]].to_numpy(float) * 1000.0

        xy, self.published = load_channel_xy(label)
        self.xy = [xy[c] for c in self.channels]
        self.heatmap = Heatmap(PHASE_DIR[label] / "heatmap_isocentre.npz")

        # The sampler must reproduce the published nominal fluences exactly.
        check = self.heatmap(*transform(self.xy, 0.0, 0.0, 0.0))
        want = np.array([self.published[c] for c in self.channels])
        worst = float(np.max(np.abs(check / want - 1.0)))
        if worst > 1e-6:
            raise RuntimeError(
                f"{label}: sampler disagrees with published fluences by {worst:.2e}")
        self.sampler_check = worst

    def score(self, dx, dy, theta):
        """Return (rss, energy, scale, fluence) for the best energy at this pose."""
        fluence = self.heatmap(*transform(self.xy, dx, dy, theta))
        if np.any(fluence <= 0):
            return np.inf, np.nan, np.nan, fluence
        response = self.kernel[:, self.channels] * (fluence / self.reference)
        obj = inf.centered_log_objective(self.measured_rad, response)
        index = int(np.argmin(obj))
        scale, _, rss = inf.fit_scale_log(self.measured_rad, response[index])
        return rss, float(self.energies[index]), scale, fluence

    def residuals(self, dx, dy, theta):
        fluence = self.heatmap(*transform(self.xy, dx, dy, theta))
        response = self.kernel[:, self.channels] * (fluence / self.reference)
        obj = inf.centered_log_objective(self.measured_rad, response)
        index = int(np.argmin(obj))
        scale, predicted, _ = inf.fit_scale_log(self.measured_rad, response[index])
        return 100.0 * (predicted / self.measured_rad - 1.0), fluence


def bic(rss, n, k):
    return n * np.log(max(rss, 1e-15) / n) + k * np.log(n)


def optimise(objective, bounds, grid_steps=(13, 13, 9)):
    """Coarse grid then Nelder-Mead, because the energy grid makes the surface kinked."""
    axes = [np.linspace(lo, hi, n) for (lo, hi), n in zip(bounds, grid_steps)]
    best, best_val = None, np.inf
    for dx in axes[0]:
        for dy in axes[1]:
            for th in axes[2]:
                v = objective((dx, dy, th))
                if v < best_val:
                    best_val, best = v, (dx, dy, th)
    # Nelder-Mead is unconstrained by default and will happily walk out of the
    # plausibility box, so the bounds are passed explicitly.
    result = minimize(objective, np.array(best), method="Nelder-Mead",
                      bounds=bounds,
                      options={"xatol": 1e-3, "fatol": 1e-10, "maxiter": 4000})
    if result.fun < best_val:
        return np.asarray(result.x, dtype=float), float(result.fun)
    return np.asarray(best, dtype=float), float(best_val)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kernel-64", type=Path, required=True)
    parser.add_argument("--kernel-85", type=Path, required=True)
    parser.add_argument("--dose-64", type=Path, required=True)
    parser.add_argument("--dose-85", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-shift-mm", type=float, default=15.0)
    parser.add_argument("--max-rotation-deg", type=float, default=10.0)
    parser.add_argument("--dx-range", type=float, nargs=2, default=None,
                        metavar=("LO", "HI"))
    parser.add_argument("--dy-range", type=float, nargs=2, default=None,
                        metavar=("LO", "HI"))
    parser.add_argument("--theta-range", type=float, nargs=2, default=None,
                        metavar=("LO", "HI"))
    parser.add_argument("--fix-pose", type=float, nargs=3, default=None,
                        metavar=("DX", "DY", "THETA"),
                        help="evaluate this pose with NO free pose parameters, "
                             "so the BIC is directly comparable to the nominal fit")
    args = parser.parse_args()

    phases = {
        "64mev": Phase("64mev", args.kernel_64, args.dose_64),
        "85mev": Phase("85mev", args.kernel_85, args.dose_85),
    }
    for label, phase in phases.items():
        print(f"{label}: sampler reproduces published fluences to "
              f"{phase.sampler_check:.2e}")

    s = args.max_shift_mm
    r = args.max_rotation_deg
    bounds = [tuple(args.dx_range) if args.dx_range else (-s, s),
              tuple(args.dy_range) if args.dy_range else (-s, s),
              tuple(args.theta_range) if args.theta_range else (-r, r)]
    report = {"bounds": {"dx_mm": list(bounds[0]), "dy_mm": list(bounds[1]),
                         "theta_deg": list(bounds[2])}, "phases": {}, "joint": {}}

    print("\n=== No pose (nominal), per-channel fluence ===")
    baseline = {}
    for label, phase in phases.items():
        rss, energy, scale, _ = phase.score(0.0, 0.0, 0.0)
        n = len(phase.channels)
        baseline[label] = (rss, n)
        print(f"  {label}: E {energy:.2f} MeV  RMS {100 * np.sqrt(rss / n):.2f} %  "
              f"C {scale / phase.reference:.3f}  BIC {bic(rss, n, 2):.2f}")
        report["phases"][label] = {"nominal": {
            "energy_MeV": energy, "rms_percent": 100 * float(np.sqrt(rss / n)),
            "calibration_ratio": scale / phase.reference, "bic": bic(rss, n, 2)}}

    if args.fix_pose is not None:
        dx, dy, theta = args.fix_pose
        print(f"\n=== Fixed pose dx {dx:+.2f} mm, dy {dy:+.2f} mm, "
              f"theta {theta:+.2f} deg -- NO free pose parameters ===")
        report["fixed_pose"] = {"dx_mm": dx, "dy_mm": dy, "theta_deg": theta,
                                "phases": {}}
        for label, phase in phases.items():
            rss, energy, scale, fluence = phase.score(dx, dy, theta)
            n = len(phase.channels)
            rss_nom = baseline[label][0]
            res_pose, _ = phase.residuals(dx, dy, theta)
            res_nom, flu_nom = phase.residuals(0.0, 0.0, 0.0)
            print(f"\n  {label}: E {energy:.2f} MeV  "
                  f"RMS {100 * np.sqrt(rss / n):.2f} % "
                  f"(nominal {100 * np.sqrt(rss_nom / n):.2f} %)  "
                  f"C {scale / phase.reference:.3f}")
            print(f"        BIC {bic(rss, n, 2):.2f} "
                  f"(nominal {bic(rss_nom, n, 2):.2f}, same parameter count)")
            print(f"  {'ch':>4s}{'resid nominal':>15s}{'resid fixed':>13s}"
                  f"{'flux nominal':>14s}{'flux fixed':>12s}")
            for i, channel in enumerate(phase.channels):
                print(f"  {channel:>4d}{res_nom[i]:>14.2f}%{res_pose[i]:>12.2f}%"
                      f"{flu_nom[i]:>14.4e}{fluence[i]:>12.4e}")
            report["fixed_pose"]["phases"][label] = {
                "energy_MeV": energy,
                "rms_percent": 100 * float(np.sqrt(rss / n)),
                "rms_percent_nominal": 100 * float(np.sqrt(rss_nom / n)),
                "calibration_ratio": scale / phase.reference,
                "bic": bic(rss, n, 2), "bic_nominal": bic(rss_nom, n, 2),
                "channels": phase.channels,
                "residual_percent_nominal": [float(v) for v in res_nom],
                "residual_percent_fixed": [float(v) for v in res_pose],
                "fluence_nominal": [float(v) for v in flu_nom],
                "fluence_fixed": [float(v) for v in fluence]}
        args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(f"\nWrote {args.output}")
        return

    print("\n=== Independent pose per phase ===")
    for label, phase in phases.items():
        pose, rss = optimise(lambda p: phase.score(*p)[0], bounds)
        _, energy, scale, fluence = phase.score(*pose)
        n = len(phase.channels)
        atbound = any(abs(abs(v) - b[1]) < 1e-3 for v, b in zip(pose, bounds))
        print(f"  {label}: dx {pose[0]:+.2f} mm  dy {pose[1]:+.2f} mm  "
              f"theta {pose[2]:+.2f} deg{'  [AT BOUND]' if atbound else ''}")
        print(f"          E {energy:.2f} MeV  RMS {100 * np.sqrt(rss / n):.2f} %  "
              f"C {scale / phase.reference:.3f}  BIC {bic(rss, n, 5):.2f}")
        report["phases"][label]["independent_pose"] = {
            "dx_mm": pose[0], "dy_mm": pose[1], "theta_deg": pose[2],
            "at_bound": atbound, "energy_MeV": energy,
            "rms_percent": 100 * float(np.sqrt(rss / n)),
            "calibration_ratio": scale / phase.reference, "bic": bic(rss, n, 5)}

    print("\n=== Joint pose, shared between phases ===")
    def joint_objective(p):
        return sum(phase.score(*p)[0] for phase in phases.values())

    pose, rss_total = optimise(joint_objective, bounds)
    n_total = sum(len(p.channels) for p in phases.values())
    rss_nominal_total = sum(baseline[l][0] for l in phases)
    atbound = any(abs(abs(v) - b[1]) < 1e-3 for v, b in zip(pose, bounds))
    print(f"  dx {pose[0]:+.2f} mm  dy {pose[1]:+.2f} mm  theta {pose[2]:+.2f} deg"
          f"{'  [AT BOUND]' if atbound else ''}")
    print(f"  pooled RMS {100 * np.sqrt(rss_total / n_total):.2f} % "
          f"(nominal {100 * np.sqrt(rss_nominal_total / n_total):.2f} %)")
    print(f"  BIC {bic(rss_total, n_total, 7):.2f} "
          f"(nominal, no pose: {bic(rss_nominal_total, n_total, 4):.2f})")
    report["joint"] = {
        "dx_mm": pose[0], "dy_mm": pose[1], "theta_deg": pose[2], "at_bound": atbound,
        "pooled_rms_percent": 100 * float(np.sqrt(rss_total / n_total)),
        "pooled_rms_percent_nominal": 100 * float(np.sqrt(rss_nominal_total / n_total)),
        "bic": bic(rss_total, n_total, 7),
        "bic_nominal_no_pose": bic(rss_nominal_total, n_total, 4),
        "per_phase": {}}

    for label, phase in phases.items():
        rss, energy, scale, fluence = phase.score(*pose)
        n = len(phase.channels)
        res_pose, _ = phase.residuals(*pose)
        res_nom, flu_nom = phase.residuals(0.0, 0.0, 0.0)
        print(f"\n  {label}: E {energy:.2f} MeV  RMS {100 * np.sqrt(rss / n):.2f} %  "
              f"C {scale / phase.reference:.3f}")
        print(f"  {'ch':>4s}{'resid nominal':>15s}{'resid pose':>12s}"
              f"{'flux nominal':>14s}{'flux pose':>12s}")
        for i, channel in enumerate(phase.channels):
            print(f"  {channel:>4d}{res_nom[i]:>14.2f}%{res_pose[i]:>11.2f}%"
                  f"{flu_nom[i]:>14.4e}{fluence[i]:>12.4e}")
        report["joint"]["per_phase"][label] = {
            "energy_MeV": energy, "rms_percent": 100 * float(np.sqrt(rss / n)),
            "calibration_ratio": scale / phase.reference,
            "channels": phase.channels,
            "residual_percent_nominal": [float(v) for v in res_nom],
            "residual_percent_pose": [float(v) for v in res_pose],
            "fluence_nominal": [float(v) for v in flu_nom],
            "fluence_pose": [float(v) for v in fluence]}

    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
