"""Inspect measured RadEx-HUS beam scans exported from PTW MEPHYSTO (.mcc)."""

from __future__ import annotations

import csv
from pathlib import Path
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PYTHON_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_ROOT))

from Dependencies.ScanMetrics import scanMetrics
from Read.ReadMCC import readMCC


DEFAULT_BASE = Path("/scratch/work/fetzera1/GRAS/RadEx/RadEx-HUS")

# Depth grid of the measured HUS water PDD, reused as the simulated layer geometry.
DEPTHS_MM = np.array([
    0.0, 1.5, 3.0, 4.5, 6.0, 7.5, 9.0, 10.5, 12.0, 13.5, 15.0,
    16.5, 18.0, 19.5, 21.0, 22.5, 24.0, 25.5, 27.0, 28.5, 30.0,
    31.5, 33.0, 34.5, 36.0, 37.5, 39.0, 40.0, 42.5, 45.0, 47.5,
    50.0, 55.0, 60.0,
])

# Half-width of the uniform source rectangle used by the HUS simulations.
MODELLED_SOURCE_HALFWIDTH_MM = 180.0

AXIS_LABELS = {
    "profile": ("Off-axis position [mm]", "Dose relative to central axis [%]"),
    "depth_dose": ("Depth in water [mm]", "Dose relative to maximum [%]"),
}


def outputStem(path):
    """Build a filesystem-safe output stem from a measurement file name."""
    return re.sub(r"[^A-Za-z0-9]+", "_", Path(path).stem).strip("_")


def metricsRow(metrics):
    """Flatten one scan's metrics into a single CSV row, keyed by metric name."""
    row = {
        "curve_type": metrics["curve_type"],
        "kind": metrics["kind"],
        "depth_mm": metrics["depth_mm"],
        "unit": metrics["unit"],
        "scan_min_mm": metrics["position"].min(),
        "scan_max_mm": metrics["position"].max(),
        "points": len(metrics["position"]),
        "reference_value_arbitrary_units": metrics["reference_value"],
        "truncated": metrics["truncated"],
    }
    if metrics["kind"] == "profile":
        for level in metrics["levels"]:
            row[f"width_{int(level)}_percent_mm"] = metrics["widths"][level]
            row[f"offset_{int(level)}_percent_mm"] = metrics["offsets"][level]
        row.update(
            flatness_percent=metrics["flatness_percent"],
            flatness_span_percent=metrics["flatness_span_percent"],
            asymmetry_percentage_points=metrics["asymmetry_percentage_points"],
            symmetry_reach_mm=metrics["symmetry_reach_mm"],
        )
    else:
        for level in metrics["levels"]:
            row[f"R{int(level)}_mm"] = metrics["ranges"][level]
        row.update(
            dmax_sampled_mm=metrics["dmax_sampled_mm"],
            surface_percent=metrics["surface_percent"],
            tail_depth_mm=metrics["tail_depth_mm"],
            tail_percent=metrics["tail_percent"],
        )
    return row


def curveLabel(curve_type):
    """Render a SCAN_CURVETYPE for display, keeping acronyms such as PDD uppercase."""
    return " ".join(word if word.isupper() and len(word) <= 4 else word.title()
                    for word in curve_type.split("_"))


def plotScans(all_metrics, scans, title, path):
    """Draw one panel per scan, annotated with the metrics of its curve type."""
    figure, axes = plt.subplots(len(all_metrics), 1,
                                figsize=(9, 4 * len(all_metrics)), squeeze=False)
    axes = axes[:, 0]
    for axis, metrics in zip(axes, all_metrics):
        label = curveLabel(metrics["curve_type"])
        axis.plot(metrics["position"], metrics["normalized"], marker=".",
                  color="tab:blue", label=label)

        if metrics["kind"] == "profile":
            annotations = [(level, f"{int(level)}%: {metrics['widths'][level]:.1f} mm")
                           for level in metrics["levels"]
                           if not np.isnan(metrics["widths"][level])]
            for sign in (-1, 1):
                axis.axvline(sign * MODELLED_SOURCE_HALFWIDTH_MM, color="tab:red",
                             linestyle="--", linewidth=1,
                             label="Modelled source edge" if sign < 0 else None)
            summary = (f"flatness ±{metrics['flatness_percent']:.2f}%, "
                       f"asymmetry {metrics['asymmetry_percentage_points']:.2f} pp")
        else:
            annotations = [(level, f"R{int(level)}: {metrics['ranges'][level]:.1f} mm")
                           for level in metrics["levels"]
                           if not np.isnan(metrics["ranges"][level])]
            axis.axvline(metrics["dmax_sampled_mm"], color="tab:red", linestyle="--",
                         linewidth=1, label=f"dmax {metrics['dmax_sampled_mm']:.1f} mm")
            summary = (f"surface {metrics['surface_percent']:.1f}%, "
                       f"tail {metrics['tail_percent']:.1f}% at "
                       f"{metrics['tail_depth_mm']:.0f} mm")

        for level, text in annotations:
            axis.axhline(level, color="0.8", linewidth=0.8, zorder=0)
            axis.annotate(text, xy=(0.01, level), xycoords=("axes fraction", "data"),
                          xytext=(0, 3), textcoords="offset points",
                          fontsize=8, color="0.35")

        depth = "" if np.isnan(metrics["depth_mm"]) else f" at {metrics['depth_mm']:.0f} mm depth"
        axis.set_title(f"{label}{depth}: {summary}"
                       + (", scan truncated" if metrics["truncated"] else ""))
        axis.set_ylim(0, 105)
        axis.set_ylabel(AXIS_LABELS[metrics["kind"]][1])
        axis.grid(True, alpha=0.3)
        axis.legend(loc="best", fontsize=8)
    axes[-1].set_xlabel(AXIS_LABELS[all_metrics[-1]["kind"]][0])
    figure.suptitle(title, fontsize=10)
    figure.tight_layout()
    figure.savefig(path)
    plt.close(figure)


def analyseScanFile(source, base=None):
    """
    Analyse every scan in one .mcc file and write its points, metrics and plot.

    Depth-dose and lateral-profile scans are dispatched by SCAN_CURVETYPE, so the
    same entry point serves the HUS water PDD and the inplane/crossplane exports.
    Derived results are written beside the measurement; the .mcc is not modified.

    Args:
        source (Path): The .mcc measurement file.
        base (Path | None): Output directory. Defaults to the file's own directory.

    Returns:
        list[dict]: The metrics of each scan, in file order.
    """
    source = Path(source)
    base = Path(base) if base is not None else source.parent
    scans = readMCC(source)
    all_metrics = [scanMetrics(scan) for scan in scans]

    stem = outputStem(source)
    points_path = base / f"Scans_{stem}.csv"
    metrics_path = base / f"ScanMetrics_{stem}.csv"
    figure_path = base / f"Scans_{stem}.pdf"

    with points_path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["curve_type", "depth_mm", "position_mm",
                         "value_arbitrary_units", "reference_arbitrary_units",
                         "normalized_percent"])
        for metrics, scan in zip(all_metrics, scans):
            reference = scan["reference"]
            if reference is not None:
                reference = reference[np.argsort(scan["position"])]
            for index, place in enumerate(metrics["position"]):
                writer.writerow([metrics["curve_type"], metrics["depth_mm"], place,
                                 metrics["value"][index],
                                 "" if reference is None else reference[index],
                                 metrics["normalized"][index]])

    rows = [metricsRow(metrics) for metrics in all_metrics]
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with metrics_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, restval="")
        writer.writeheader()
        writer.writerows(rows)

    plotScans(all_metrics, scans, f"HUS 6 eHDTSE measured scans, {source.name}", figure_path)
    print(f"{source.name}: {len(all_metrics)} scan(s)")
    for metrics in all_metrics:
        if metrics["kind"] == "profile":
            widths = ", ".join(f"{int(level)}% {metrics['widths'][level]:.1f} mm"
                               for level in metrics["levels"])
            print(f"  {metrics['curve_type']}: {widths}, "
                  f"flatness ±{metrics['flatness_percent']:.2f}%, "
                  f"asymmetry {metrics['asymmetry_percentage_points']:.2f} pp")
        else:
            ranges = ", ".join(f"R{int(level)} {metrics['ranges'][level]:.2f} mm"
                               for level in metrics["levels"])
            print(f"  {metrics['curve_type']}: dmax {metrics['dmax_sampled_mm']:.1f} mm, "
                  f"{ranges}, surface {metrics['surface_percent']:.1f}%, "
                  f"tail {metrics['tail_percent']:.1f}%")
    print(f"  wrote {points_path.name}, {metrics_path.name}, {figure_path.name}")
    return all_metrics


def main():
    targets = [Path(argument) for argument in sys.argv[1:]]
    if not targets:
        targets = sorted(DEFAULT_BASE.glob("*.mcc"))
    if not targets:
        raise RuntimeError(f"No .mcc files found in {DEFAULT_BASE}")
    for source in targets:
        analyseScanFile(source)


if __name__ == "__main__":
    main()
