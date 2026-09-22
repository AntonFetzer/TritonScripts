"""Compare every run of one CRDS campaign, both instruments side by side.

The runs of a campaign share the detector set, the cuts and the spectrum and
differ in one controlled thing, so their coefficients are directly comparable
even though the runs differ in history count: GRAS normalises to one incident
particle/cm2. The runs and the steps compared are defined per campaign in
Dependencies/CRDSCampaigns.py.

HUS is a field-size and table series. At 20 x 20 mm a ~5 MeV electron field is
below lateral scatter equilibrium, so dose on the axis is depressed relative to
a broad field; the clinical applicator was 36 x 36 cm. The 100 mm step measures
what the specified 20 mm field costs, the 200 mm step asks whether 100 mm had
already saturated, and the no-table run bounds the error made by modelling a
carbon-fibre-and-foam couch as solid PMMA.

Kumpula is an epoxy sweep: eight filed-LED thicknesses bracketing a residue
that was never measured, each compared against 500 um. The RadFET was not
modified, so its TID is a control that must not move; the LED DDD, and the
active/substrate ratio, are what the sweep is for. A campaign that varies one
numeric parameter also gets a consistency test across the runs and a plot of
every tally against that parameter.

Uppsala is two exposures with different fields, neither uniform over the
instrument. The 52 x 54 mm run against the 70 x 74 mm run is a field-size test
on the RadFET; the LED is outside the narrow simulated field and is flagged
there. The campaign's ``delivery`` block then combines the coefficients with the
local fluence at each part from the accelerator-log maps into delivered TID and
DDD, written to its own CSV; see Dependencies/CRDSDelivery.py.

Both instruments are carried through, because they do not have to agree: the
RadFET and the LED sit at different positions and depths, and displacement
damage weights the spectrum differently from ionising dose.

    python3 Plotting/CRDSSeries.py --campaign HUS
    python3 Plotting/CRDSSeries.py --campaign Kumpula
    python3 Plotting/CRDSSeries.py --campaign Uppsala

Aggregation is Dependencies.AggregateRun and Dependencies.TotalNID; nothing
here re-parses a GRAS file.
"""

import argparse
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from Dependencies.HistogramPlots import colour, logLogGrid, savePdf  # noqa: E402
from Dependencies.CRDSCampaigns import CAMPAIGNS, MODULES, TILES, campaign  # noqa: E402
from Dependencies.CRDSSeriesAnalysis import (collectSeries, printSeries,  # noqa: E402
    ratioWithError, reportDelivery, weightedMean, writeSeries)


def dddPanel(axis, spec, data, labels, x):
    """Draw the LED DDD coefficients of a sweep on one axis."""
    for index, (_, _, name) in enumerate(MODULES):
        values = [data[label][("DDD", name)] for label in labels]
        axis.errorbar(x, [v[0] for v in values], yerr=[v[1] for v in values],
                      fmt="o-", capsize=3, color=colour(index), label=name)
    axis.set_ylabel(f"DDD [MeV/g per {spec['particle']}/cm2]")
    axis.set_ylim(bottom=0)
    axis.set_title(f"CRDS {spec['description']}")


def plotSweep(spec, data, output):
    """Plot every tally of a one-parameter campaign against that parameter.

    Three panels on a shared axis: the LED DDD coefficients, the LED
    active/substrate ratio, and the RadFET TID as a deviation from its
    inverse-variance mean over the sweep, which is how a control that should
    not move is best read. The DDD panel is also saved on its own, as
    ``<output stem>_DDD.pdf``.

    Returns:
        list[Path]: the files written.
    """
    sweep, runs, particle = spec["sweep"], spec["runs"], spec["particle"]
    labels = sorted(runs, key=lambda label: runs[label][sweep["key"]])
    x = [runs[label][sweep["key"]] for label in labels]

    figure, axis = plt.subplots(figsize=(6.4, 3.2))
    dddPanel(axis, spec, data, labels, x)
    axis.set_xlabel(sweep["xlabel"])
    logLogGrid(axis, xscale="linear", yscale="linear")
    standalone = output.with_name(f"{output.stem}_DDD.pdf")
    savePdf(standalone, figure)
    plt.close(figure)

    figure, axes = plt.subplots(3, 1, sharex=True, figsize=(6.4, 9.6))
    dddPanel(axes[0], spec, data, labels, x)

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
    return [output, standalone]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--campaign", required=True, choices=list(CAMPAIGNS))
    parser.add_argument("--reference-fluence", type=float,
                        help="particles/cm2; defaults to the campaign's")
    arguments = parser.parse_args()

    spec = campaign(arguments.campaign)
    reference = arguments.reference_fluence or spec["referenceFluence"]
    data, rows = collectSeries(spec, reference)
    printSeries(arguments.campaign, spec, data, reference)
    print(f"\nWrote {writeSeries(spec, rows)}")
    if spec["sweep"]:
        for written in plotSweep(spec, data, spec["path"] / spec["sweep"]["output"]):
            print(f"Wrote {written}")
    if spec.get("delivery"):
        print(f"Wrote {reportDelivery(spec, data, spec.get('notCoefficients', set()))}")


if __name__ == "__main__":
    main()
