"""Compare a simulation variant with a reference run and with a campaign's measured doses.

    python3 -m RadEx.CompareVariant <variant run folder> --campaign <json> [--reference <run folder>]

Reads TotalDose_<run>.csv of both runs (RadEx/RadExTID.py output) and the campaign's
measured doses of its measurement phase (RadEx.KernelFit; sigma^2 = noise^2 +
conversion^2 + (systematic x dose)^2). The reference defaults to the campaign's
reference_run. Writes Comparison_<run>.csv and Comparison_<run>.pdf into the variant
folder: the per-channel variant/reference ratios of every layer, and for each run the
per-type calibration-scale fit of the default layer to the measured doses, with its
chi2 and pulls.
"""
import argparse
import csv
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from RadEx import KernelFit as kf  # noqa: E402


def total_dose(folder):
    rows = np.loadtxt(folder / f"TotalDose_{folder.name}.csv", delimiter=",", skiprows=1)
    return rows[:, 1], rows[:, 2]


def scale_fit(coef, measured, sigma, fluence, channels):
    D, S = measured[channels], sigma[channels]
    m = coef[channels] * fluence
    s = float(np.sum(D * m / S ** 2) / np.sum(m ** 2 / S ** 2))
    pulls = (D - s * m) / S
    return s, pulls, float(np.sum(pulls ** 2))


parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
parser.add_argument("variant", type=Path)
parser.add_argument("--campaign", type=Path, required=True)
parser.add_argument("--reference", type=Path, help="reference run folder (default: the campaign's reference_run)")
args = parser.parse_args()
campaign = kf.load_campaign(args.campaign)
variant = args.variant.resolve()
reference = (args.reference or Path(campaign["reference_run"])).resolve()
channels = campaign["channels"]
phase = campaign["measurement_phase"]
meas = kf.measurement(campaign, phase)
cfg = campaign["phases"][phase]
measured, sigma, fluence = meas["value"], meas["sigma"], cfg["log_fluence"]
sample = {c: str(cfg["table"].loc[c].get("sample", "")) for c in channels}
TYPE = campaign["device"]
types = {t: [c for c in channels if TYPE[c] == t] for t in kf.device_types(campaign)}
layers = sorted(campaign["layers"], key=campaign["layers"].get)
fitLayer = kf.layer_slice(campaign, campaign["default_layer"])

dv, ev = total_dose(variant)
dr, er = total_dose(reference)
vname, rname = variant.name, reference.name
ratio = dv / dr
ratio_err = ratio * np.sqrt((ev / dv) ** 2 + (er / dr) ** 2)

lines = [f"# {vname} against {rname}; measured {phase} doses from {cfg['dose'].parent.name}/{cfg['dose'].name}"]
fits = {}
for run, coef in ((rname, dr[fitLayer]), (vname, dv[fitLayer])):
    fits[run] = {}
    chi2 = 0.0
    for t, chans in types.items():
        s, pulls, c2 = scale_fit(coef, measured, sigma, fluence, chans)
        fits[run][t] = (s, dict(zip(chans, pulls)), c2)
        chi2 += c2
    fits[run]["chi2"] = chi2
    dof = len(channels) - len(types)
    lines.append(f"# {run}: scales " + ", ".join(f"{t} {fits[run][t][0]:.3f}" for t in types)
                 + f"; chi2 {chi2:.1f} for {dof} dof")

with open(variant / f"Comparison_{vname}.csv", "w", newline="") as f:
    for line in lines:
        f.write(line + "\n")
    w = csv.writer(f)
    header = ["channel", "type", "sample"]
    for layer in layers:
        header += [f"ref_{layer}", f"var_{layer}", f"{layer}_ratio", f"{layer}_ratio_err"]
    w.writerow(header + [f"measured_{phase}_krad", "sigma_krad", "ref_pull", "var_pull"])
    for c in channels:
        row = [c, TYPE[c], sample[c]]
        for layer in layers:
            k = campaign["layers"][layer] + c
            row += [f"{dr[k]:.4e}", f"{dv[k]:.4e}", f"{ratio[k]:.4f}", f"{ratio_err[k]:.4f}"]
        w.writerow(row + [f"{measured[c]:.2f}", f"{sigma[c]:.2f}",
                          f"{fits[rname][TYPE[c]][1][c]:+.2f}", f"{fits[vname][TYPE[c]][1][c]:+.2f}"])

for line in lines:
    print(line)
print("ch  type  sample           " + "".join(f"{layer + ' ratio':16s}" for layer in layers) + "ref pull  var pull")
for c in channels:
    print(f"{c:2d}  {TYPE[c]}  {sample[c]:<16} "
          + "".join(f"{ratio[campaign['layers'][layer] + c]:.3f}+-{ratio_err[campaign['layers'][layer] + c]:.3f}   " for layer in layers)
          + f"{fits[rname][TYPE[c]][1][c]:+.2f}     {fits[vname][TYPE[c]][1][c]:+.2f}")

fig, (a, b) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
x = np.arange(len(channels))
for k, (layer, marker) in enumerate(zip(layers, "osd^")):
    part = kf.layer_slice(campaign, layer)
    a.errorbar(x + 0.2 * (k - 0.5 * (len(layers) - 1)), ratio[part], ratio_err[part], fmt=marker, label=layer)
a.axhline(1, color="k", lw=0.8)
a.set_ylabel(f"{vname} / {rname}")
a.legend()
a.grid(alpha=0.3)
b.plot(x - 0.1, [fits[rname][TYPE[c]][1][c] for c in channels], "o", label=f"{rname}, chi2 {fits[rname]['chi2']:.1f}")
b.plot(x + 0.1, [fits[vname][TYPE[c]][1][c] for c in channels], "s", label=f"{vname}, chi2 {fits[vname]['chi2']:.1f}")
b.axhline(0, color="k", lw=0.8)
b.set_ylabel("pull (measured - scale x model) / sigma")
b.set_xlabel("channel")
b.set_xticks(x)
b.set_xticklabels([f"{c}\n{sample[c][:8]}" for c in channels], fontsize=7)
b.legend()
b.grid(alpha=0.3)
fig.suptitle(f"{phase.capitalize()}-phase measured doses: per-type calibration scale fits ({campaign['default_layer']} tiles)")
fig.tight_layout()
fig.savefig(variant / f"Comparison_{vname}.pdf")
print("wrote", variant / f"Comparison_{vname}.csv", "and .pdf")
