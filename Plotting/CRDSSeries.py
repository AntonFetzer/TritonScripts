"""Compare every run of one CRDS campaign, both instruments side by side.

The runs of a campaign share the detector set, the cuts and the spectrum and
differ in one controlled thing, so their coefficients are directly comparable
even though the runs differ in history count: GRAS normalises to one incident
particle/cm2. The runs and the steps compared are defined per campaign in
Plotting/CRDSCampaigns.py.

HUS is a field-size and table series. At 20 x 20 mm a ~5 MeV electron field is
below lateral scatter equilibrium, so dose on the axis is depressed relative to
a broad field; the clinical applicator was 36 x 36 cm. The 100 mm step measures
what the specified 20 mm field costs, the 200 mm step asks whether 100 mm had
already saturated, and the no-table run bounds the error made by modelling a
carbon-fibre-and-foam couch as solid PMMA.

Kumpula is an epoxy sweep: five filed-LED thicknesses bracketing a residue that
was never measured, each compared against the thickest. The RadFET was not
modified, so its TID is a control that must not move; the LED DDD, and the
active/substrate ratio, are what the sweep is for. A campaign that varies one
numeric parameter also gets a consistency test across the runs and a plot of
every tally against that parameter.

Both instruments are carried through, because they do not have to agree: the
RadFET and the LED sit 10 mm apart at different depths, and displacement damage
weights the spectrum differently from ionising dose.

    python3 Plotting/CRDSSeries.py --campaign HUS
    python3 Plotting/CRDSSeries.py --campaign Kumpula

Aggregation is Dependencies.AggregateRun and Dependencies.TotalNID; nothing
here re-parses a GRAS file.
"""

import argparse
import csv
import math
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from Dependencies.AggregateRun import aggregateRun  # noqa: E402
from Dependencies.HistogramPlots import colour, logLogGrid, savePdf  # noqa: E402
from Dependencies.TotalNID import totalNID  # noqa: E402
from Plotting.CRDSCampaigns import CAMPAIGNS, MODULES, TILES, campaign  # noqa: E402


def ratioWithError(numerator, numeratorError, denominator, denominatorError):
    """Ratio of two independent measurements and its absolute error."""
    ratio = numerator / denominator
    relative = math.sqrt((numeratorError / numerator) ** 2
                         + (denominatorError / denominator) ** 2)
    return ratio, ratio * relative


def weightedMean(values, errors):
    """Inverse-variance mean, its error, and chi2 per degree of freedom about it."""
    weights = [1.0 / error ** 2 for error in errors]
    mean = sum(w * v for w, v in zip(weights, values)) / sum(weights)
    chi2 = sum(w * (v - mean) ** 2 for w, v in zip(weights, values))
    return mean, sum(weights) ** -0.5, chi2 / (len(values) - 1)


def plotSweep(spec, data, output):
    """Plot every tally of a one-parameter campaign against that parameter.

    Three panels on a shared axis: the LED DDD coefficients, the LED
    active/substrate ratio, and the RadFET TID as a deviation from its
    inverse-variance mean over the sweep, which is how a control that should
    not move is best read.
    """
    sweep, runs, particle = spec["sweep"], spec["runs"], spec["particle"]
    labels = sorted(runs, key=lambda label: runs[label][sweep["key"]])
    x = [runs[label][sweep["key"]] for label in labels]

    figure, axes = plt.subplots(3, 1, sharex=True, figsize=(6.4, 9.6))

    for index, (_, _, name) in enumerate(MODULES):
        values = [data[label][("DDD", name)] for label in labels]
        axes[0].errorbar(x, [v[0] for v in values], yerr=[v[1] for v in values],
                         fmt="o-", capsize=3, color=colour(index), label=name)
    axes[0].set_ylabel(f"DDD [MeV/g per {particle}/cm2]")
    axes[0].set_ylim(bottom=0)
    axes[0].set_title(f"CRDS {spec['description']}")

    ratios = [ratioWithError(*data[label][("DDD", MODULES[0][2])][:2],
                             *data[label][("DDD", MODULES[1][2])][:2])
              for label in labels]
    axes[1].errorbar(x, [r[0] for r in ratios], yerr=[r[1] for r in ratios],
                     fmt="o-", capsize=3, color=colour(2))
    axes[1].axhline(1.0, color="k", linewidth=0.8)
    axes[1].set_ylabel(f"LED {MODULES[0][2]} / {MODULES[1][2]}")

    for index, (_, _, name) in enumerate(TILES):
        values = [data[label][("TID", name)] for label in labels]
        mean, _, _ = weightedMean([v[0] for v in values], [v[1] for v in values])
        axes[2].errorbar(x, [100 * (v[0] / mean - 1) for v in values],
                         yerr=[100 * v[1] / mean for v in values],
                         fmt="o-", capsize=3, color=colour(index + 3),
                         label=f"{name}, mean {mean:.4e} kRad cm2/{particle}")
    axes[2].axhline(0.0, color="k", linewidth=0.8)
    axes[2].set_ylabel("RadFET TID vs sweep mean [%]")
    axes[2].set_xlabel(sweep["xlabel"])

    for axis in axes:
        logLogGrid(axis, xscale="linear", yscale="linear")
    savePdf(output, figure)
    plt.close(figure)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--campaign", required=True, choices=list(CAMPAIGNS))
    parser.add_argument("--reference-fluence", type=float,
                        help="particles/cm2; defaults to the campaign's")
    arguments = parser.parse_args()

    spec = campaign(arguments.campaign)
    particle = spec["particle"]
    reference = arguments.reference_fluence or spec["referenceFluence"]
    runs = spec["runs"]

    data, rows = {}, []
    for label, run in runs.items():
        resultPath = spec["path"] / run["folder"] / "Res"
        tid = aggregateRun(str(resultPath), expectedFiles=run["files"],
                           tileCount=len(TILES))
        ddd = totalNID(str(resultPath), expectedFiles=run["files"],
                       expectedModules=len(MODULES))
        entries = {}
        for index, volume, name in TILES:
            entries[("TID", name)] = (float(tid["dose"][index]),
                                      float(tid["error"][index]),
                                      f"kRad cm2/{particle}", volume,
                                      int(tid["entries"][index]),
                                      int(tid["non-zeros"][index]))
        for module, volume, name in MODULES:
            entry = ddd[module]
            entries[("DDD", name)] = (entry["nid"], entry["error"],
                                      f"{entry['unit']} cm2/{particle}", volume,
                                      entry["entries"], entry["non-zeros"])
        data[label] = entries

        for (kind, name), (value, error, unit, volume, n, nz) in entries.items():
            rows.append({
                "run": label, "folder": run["folder"], "slurm_job_id": run["job"],
                "quantity": kind, "tally": name, "volume_name": volume,
                "simulated_primaries": n,
                "nonzero_entries": nz,
                f"coefficient_per_{particle}_cm2": f"{value:.12e}",
                "statistical_error": f"{error:.12e}",
                "relative_error_percent": f"{100 * error / value:.6f}",
                "reference_fluence": f"{reference:.6e}" if reference else "",
                "at_reference_fluence":
                    f"{value * reference:.6e}" if reference else "",
            })

    width = max(len(k) for k in runs) + 2
    atHeader = f"at {reference:.3g}" if reference else ""
    print(f"\n{arguments.campaign}: {spec['description']}")
    print(f"\n{'run':<{width}}{'quantity':<10}{'tally':<16}{'coefficient':>16}"
          f"{'rel err':>10}{atHeader:>16}")
    for label in runs:
        for (kind, name), (value, error, unit, _, _, _) in data[label].items():
            shown = ""
            if reference:
                at = value * reference
                shown = f"{at:.3f} kRad" if kind == "TID" else f"{at:.4g} MeV/g"
            print(f"{label:<{width}}{kind:<10}{name:<16}{value:>16.6e}"
                  f"{100 * error / value:>9.2f}%{shown:>16}")

    # Within-run ratios. Both tallies of a pair see the same primaries, so
    # their errors are correlated and this quadrature error is an upper bound.
    print(f"\n{'run':<{width}}{'oxide / die':>20}{'active / substrate':>24}")
    for label in runs:
        pairs = [(data[label][("TID", TILES[0][2])], data[label][("TID", TILES[1][2])]),
                 (data[label][("DDD", MODULES[0][2])], data[label][("DDD", MODULES[1][2])])]
        shown = []
        for a, b in pairs:
            ratio, error = ratioWithError(a[0], a[1], b[0], b[1])
            shown.append(f"{ratio:.4f} +- {error:.4f}")
        print(f"{label:<{width}}{shown[0]:>20}{shown[1]:>24}")

    print("\nsteps, each against the run named second")
    for numerator, denominator in spec["steps"]:
        print(f"  {numerator}  vs  {denominator}")
        for key in data[numerator]:
            a, ea = data[numerator][key][0], data[numerator][key][1]
            b, eb = data[denominator][key][0], data[denominator][key][1]
            ratio, error = ratioWithError(a, ea, b, eb)
            sigma = abs(ratio - 1) / error if error else float("inf")
            print(f"    {key[0]:<4}{key[1]:<16}{100 * (ratio - 1):>+8.2f}%"
                  f" +- {100 * error:.2f}%   ({sigma:.1f} sigma)")

    if spec["sweep"]:
        # A tally that does not depend on the swept parameter scatters about
        # one mean with chi2/dof near 1; the RadFET is that control here.
        print(f"\nconsistency across the sweep, {len(runs) - 1} degrees of freedom")
        for key in data[next(iter(runs))]:
            values = [data[label][key] for label in runs]
            mean, meanError, reduced = weightedMean([v[0] for v in values],
                                                    [v[1] for v in values])
            spread = 100 * (max(v[0] for v in values)
                            - min(v[0] for v in values)) / mean
            print(f"    {key[0]:<4}{key[1]:<16}mean {mean:.6e} +- {meanError:.2e}"
                  f"   chi2/dof {reduced:6.2f}   max-min {spread:6.2f}%")

    output = spec["path"] / spec["seriesOutput"]
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {output}")
    if spec["sweep"]:
        print(f"Wrote {plotSweep(spec, data, spec['path'] / spec['sweep']['output'])}")


if __name__ == "__main__":
    main()
