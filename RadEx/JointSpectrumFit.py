"""Joint fit of one electron spectrum to a measured water PDD and a campaign's RadFET doses.

    python3 -m RadEx.JointSpectrumFit --kernel-dir <sweep> --campaign <json> --output-dir <dir>
        [--layer oxide] [--include-excluded] [--model 0.10] [--per-type]

Unknown: non-negative fluence weights x on the energy grid of the water response
matrix named in the campaign's joint_fit block (columns E<energy>_MeV, rows depths),
smoothed with lambda ||D2 x||^2 as in the PDD-only fit. The PDD block is A x ~ y (the
measured PDD in % of its maximum up to depth_max_mm, one free scale absorbed in x,
sigma sigma_pdd_pp). The RadFET block is r F K x ~ d: the layer's kernel K
(interpolated to the water grid), the measurement phase's fluence F, the measured
doses d with sigma^2 = noise^2 + conversion^2 + (model x d)^2, and one amplitude r
(shared calibration scale) per group; per-type scales are a variant. The campaign's
excluded_channels are left out unless --include-excluded.

Profiles: the fraction of fluence below --low-mev is forced to a grid of values through
a heavy constraint row, and the PDD and RadFET chi2 are recorded, to show which
fractions each data set allows. Writes joint-fit-<tag>.json/.csv/.pdf into --output-dir.
"""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from Dependencies.EnergyKernel import interpolateKernel, readKernelRuns  # noqa: E402
from RadEx.DoseCampaign import deviceTypes, layerSlice, loadCampaign, phaseDoses  # noqa: E402

try:
    from scipy.optimize import nnls as _nnls

    def nnls(A, b):
        return _nnls(A, b, maxiter=50 * A.shape[1])[0]
except ImportError:  # the repository's Lawson-Hanson implementation
    from Dependencies.EnergyKernel import nnls

FRACTIONS = np.array([0.02, 0.04, 0.06, 0.085, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50])

parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("--kernel-dir", type=Path, required=True)
parser.add_argument("--campaign", type=Path, required=True)
parser.add_argument("--output-dir", type=Path, required=True)
parser.add_argument("--layer", help="detector layer of the campaign (default: its default_layer)")
parser.add_argument("--include-excluded", action="store_true", help="fit the campaign's excluded_channels too")
parser.add_argument("--model", type=float, help="model fraction in the RadFET sigma (default: the campaign's)")
parser.add_argument("--per-type", action="store_true", help="one calibration scale per RadFET type")
parser.add_argument("--low-mev", type=float, default=3.0, help="edge of the profiled low-energy fraction [MeV]")
args = parser.parse_args()
campaign = loadCampaign(args.campaign)
joint_cfg = campaign["joint_fit"]
layer = args.layer or campaign["default_layer"]
model = campaign["systematic_fraction"] if args.model is None else args.model
SIGMA_PDD, DEPTH_MAX, LAMBDA = joint_cfg["sigma_pdd_pp"], joint_cfg["depth_max_mm"], joint_cfg["smoothing"]
excluded = [] if args.include_excluded else list(joint_cfg.get("excluded_channels", []))
low = f"{args.low_mev:g}MeV"
tag = (layer + ("-with-" + "-".join(map(str, joint_cfg.get("excluded_channels", []))) if args.include_excluded else "")
       + (f"-model{model:.2f}" if model != campaign["systematic_fraction"] else "") + ("-per-type" if args.per_type else ""))
out = args.output_dir
out.mkdir(parents=True, exist_ok=True)

# ---- water PDD block ---------------------------------------------------------------------
rows = list(csv.reader(open(joint_cfg["response_matrix"])))
grid = np.array([float(h[1:-4]) for h in rows[0][1:]])            # E0.25_MeV -> 0.25
depth = np.array([float(r[0]) for r in rows[1:]])
A = np.array([[float(v) for v in r[1:]] for r in rows[1:]])
pred = [r for r in csv.DictReader(open(joint_cfg["measured_pdd"]))
        if all(r[k] == v for k, v in joint_cfg["measured_pdd_selection"].items())]
y = np.array([float(r["measured_percent"]) for r in sorted(pred, key=lambda r: float(r["depth_mm"]))])
use = depth <= DEPTH_MAX
A_pdd, y_pdd = A[use], y[use]
n = len(grid)
D2 = np.zeros((n - 2, n))
for i in range(n - 2):
    D2[i, i:i + 3] = (1, -2, 1)
below3 = grid < args.low_mev

# ---- RadFET block -----------------------------------------------------------------------
energies, dose, error, _ = readKernelRuns(args.kernel_dir)
K = interpolateKernel(energies, dose[:, layerSlice(campaign, layer)], grid)   # (grid, channel), dose per e/cm2
doses = phaseDoses(campaign, campaign["measurement_phase"], model)
d, sig, F = doses["measured"], doses["sigma"], doses["fluence"]
nChannels = len(campaign["channels"])
channels = [c for c in campaign["channels"] if c not in excluded]
TYPE = campaign["devices"]
groups = {t: [c for c in channels if TYPE[c] == t] for t in deviceTypes(campaign)} if args.per_type else {"all": channels}


def solve(r, lam=LAMBDA, fraction=None, pdd_weight=1.0):
    """NNLS on the stacked blocks; r: amplitude per group (dict). Returns x and the chi2 parts."""
    blocks_A, blocks_b = [A_pdd * pdd_weight / SIGMA_PDD], [y_pdd * pdd_weight / SIGMA_PDD]
    for g, idx in groups.items():
        blocks_A.append((r[g] * F * K[:, idx]).T / sig[idx][:, None])
        blocks_b.append(d[idx] / sig[idx])
    blocks_A.append(np.sqrt(lam) * D2 * np.sqrt(len(y_pdd)) / SIGMA_PDD)   # same relative strength as the PDD-only fit
    blocks_b.append(np.zeros(n - 2))
    if fraction is not None:
        row = np.where(below3, 1.0 - fraction, -fraction) * 1e3 / SIGMA_PDD
        blocks_A.append(row[None, :]); blocks_b.append(np.zeros(1))
    x = nnls(np.vstack(blocks_A), np.concatenate(blocks_b))
    chi_pdd = float(np.sum(((A_pdd @ x - y_pdd) / SIGMA_PDD) ** 2))
    chi_rf = 0.0
    pulls = np.zeros(nChannels)
    for g, idx in groups.items():
        p = r[g] * F * (K[:, idx].T @ x)
        pulls[idx] = (p - d[idx]) / sig[idx]
        chi_rf += float(np.sum(pulls[idx] ** 2))
    return x, chi_pdd, chi_rf, pulls


def profile(fraction=None, lam=LAMBDA, pdd_weight=1.0, iterations=12):
    """Minimise the total chi2: NNLS for x at fixed amplitudes, then the weighted
    least-squares amplitude of each group at fixed x, alternated to convergence."""
    r = {g: 1.0 for g in groups}
    x = nnls(np.vstack([A_pdd / SIGMA_PDD, np.sqrt(lam) * D2 * np.sqrt(len(y_pdd)) / SIGMA_PDD]), np.concatenate([y_pdd / SIGMA_PDD, np.zeros(n - 2)]))
    for _ in range(iterations):
        for g, idx in groups.items():
            k = F * (K[:, idx].T @ x)
            r[g] = float(np.sum(d[idx] * k / sig[idx] ** 2) / np.sum(k * k / sig[idx] ** 2))
        x_new, cp, cr, pulls = solve(r, lam, fraction, pdd_weight)
        done = np.allclose(x_new, x, rtol=1e-6, atol=1e-12)
        x = x_new
        if done:
            break
    f = x / x.sum()
    return {"x": x, "fraction": f, "chi2_pdd": cp, "chi2_radfet": cr, "pulls": pulls, "amplitude": r,
            "calibration": {g: r[g] * x.sum() for g in groups},
            "mean_energy_MeV": float(f @ grid), "below_" + low: float(f[below3].sum()), "below_2MeV": float(f[grid < 2].sum()),
            "above_7MeV": float(f[grid > 7].sum()), "surface_percent": float((A[0] @ x)), "rmse_pdd_pp": float(np.sqrt(np.mean((A_pdd @ x - y_pdd) ** 2)))}


results = {"tag": tag, "layer": layer, "channels": channels, "model_fraction": model, "sigma_pdd_pp": SIGMA_PDD, "lambda": LAMBDA, "depth_max_mm": DEPTH_MAX}
# PDD only (RadFET block weight zero) and the joint fits
pdd_only = profile(pdd_weight=1.0)
results["pdd_only_reference"] = joint_cfg.get("pdd_only_reference")
joint = profile()
results["joint"] = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in joint.items()}
print(f"joint {tag}: mean {joint['mean_energy_MeV']:.3f} MeV, <{low} {joint['below_' + low]:.3f}, <2 MeV {joint['below_2MeV']:.3f}, "
      f"chi2 PDD {joint['chi2_pdd']:.1f}/{use.sum()} (rmse {joint['rmse_pdd_pp']:.3f} pp), RadFET {joint['chi2_radfet']:.1f}/{len(channels)}, "
      f"calibration {joint['calibration']}")
print("   pulls " + " ".join(f"Ch{c}:{joint['pulls'][c]:+.2f}" for c in channels))
# profile over the fraction below 3 MeV
prof = []
for phi in FRACTIONS:
    p = profile(fraction=phi)
    prof.append({"fraction_below_" + low: float(phi), "mean_energy_MeV": p["mean_energy_MeV"], "chi2_pdd": p["chi2_pdd"], "chi2_radfet": p["chi2_radfet"],
                 "surface_percent": p["surface_percent"], "rmse_pdd_pp": p["rmse_pdd_pp"], "calibration": p["calibration"],
                 **{f"pull_ch{c}": (float(p["pulls"][c]) if c in channels else None) for c in joint_cfg.get("excluded_channels", [])},
                 "spectrum": p["fraction"].tolist()})
    print(f"   forced <{low} {phi:.3f}: mean {p['mean_energy_MeV']:.2f} chi2 PDD {p['chi2_pdd']:.1f} RadFET {p['chi2_radfet']:.1f} surface {p['surface_percent']:.2f} %"
          + "".join(f" Ch{c} {p['pulls'][c]:+.1f}" for c in joint_cfg.get("excluded_channels", []) if c in channels))
results["fraction_profile"] = prof
with open(out / f"joint-fit-{tag}.json", "w") as handle:
    json.dump(results, handle, indent=2)
with open(out / f"joint-fit-{tag}.csv", "w", newline="") as handle:
    w = csv.writer(handle)
    w.writerow(["energy_MeV", "pdd_only_fraction", "joint_fraction"] + [f"forced_below{args.low_mev:g}_{phi:.3f}" for phi in FRACTIONS])
    for i, e in enumerate(grid):
        w.writerow([e, f"{pdd_only['fraction'][i]:.5f}", f"{joint['fraction'][i]:.5f}"] + [f"{p['spectrum'][i]:.5f}" for p in prof])

fig, ax = plt.subplots(2, 2, figsize=(11, 8))
ax[0, 0].step(grid, pdd_only["fraction"], where="mid", label=f"PDD only: {pdd_only['mean_energy_MeV']:.2f} MeV")
ax[0, 0].step(grid, joint["fraction"], where="mid", label=f"joint: {joint['mean_energy_MeV']:.2f} MeV, chi2 PDD {joint['chi2_pdd']:.0f}, RadFET {joint['chi2_radfet']:.1f}")
for phi in (0.20, 0.40):
    p = prof[[round(q["fraction_below_" + low], 3) for q in prof].index(phi)]
    ax[0, 0].step(grid, p["spectrum"], where="mid", ls="--", label=f"forced {phi:.0%} below {args.low_mev:g} MeV: chi2 PDD {p['chi2_pdd']:.0f}, RadFET {p['chi2_radfet']:.1f}")
ax[0, 0].set_xlabel("energy (MeV)"); ax[0, 0].set_ylabel(f"fluence fraction per {np.median(np.diff(grid)):g} MeV"); ax[0, 0].legend(fontsize=7); ax[0, 0].grid(alpha=0.3)
ax[0, 1].plot([p["fraction_below_" + low] for p in prof], [p["chi2_pdd"] for p in prof], "o-", label=f"PDD ({use.sum()} depths, sigma {SIGMA_PDD} pp)")
ax[0, 1].plot([p["fraction_below_" + low] for p in prof], [p["chi2_radfet"] for p in prof], "s-", label=f"RadFET ({len(channels)} channels)")
reference = (joint_cfg.get("pdd_only_reference") or {}).get("below_" + low)
if reference is not None:
    ax[0, 1].axvline(reference, color="k", lw=0.8, ls=":", label=f"PDD-only spectrum ({reference:.1%})")
ax[0, 1].set_xlabel(f"fluence fraction below {args.low_mev:g} MeV (forced)"); ax[0, 1].set_ylabel("chi2"); ax[0, 1].set_yscale("log"); ax[0, 1].legend(fontsize=7); ax[0, 1].grid(alpha=0.3)
x = joint["x"]
ax[1, 0].plot(depth, y, "k.", label="measured PDD")
ax[1, 0].plot(depth, A @ x, "-", label="joint spectrum")
ax[1, 0].plot(depth, A @ pdd_only["x"], "--", label="PDD-only spectrum")
ax[1, 0].set_xlabel("depth in water (mm)"); ax[1, 0].set_ylabel("PDD (%)"); ax[1, 0].legend(fontsize=7); ax[1, 0].grid(alpha=0.3)
ax[1, 1].bar(np.arange(nChannels), joint["pulls"], color=["C3" if c not in channels else "C0" for c in range(nChannels)])
ax[1, 1].set_xticks(np.arange(nChannels)); ax[1, 1].set_xlabel("channel (red: not in the fit)"); ax[1, 1].set_ylabel("RadFET pull, joint spectrum"); ax[1, 1].grid(alpha=0.3)
fig.suptitle(f"Joint PDD + RadFET spectrum fit, {tag}")
fig.tight_layout()
fig.savefig(out / f"joint-fit-{tag}.pdf")
print("wrote", out / f"joint-fit-{tag}.*")
