#!/usr/bin/env python3
"""Aggregate and plot the RadEx Uppsala mono-energy response kernel."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# The legacy-detector kernel sweep this script aggregates (it used to live there).
ROOT = Path("/scratch/work/fetzera1/GRAS/RadEx/RadEx-Uppsala/EnergyKernel-FinalGeometry")
ENERGY_MANIFEST = ROOT / "energies.csv"
CHANNELS = range(12)
RUNS_PER_ENERGY = 10
PRIMARIES_PER_RUN = 4_000_000
SUMMARY_PATH = ROOT / "energy-kernel-dose-summary.csv"
PLOT_DIR = ROOT / "Plots"
PLOT_STEM = PLOT_DIR / "dose-vs-energy-all-channels"


def read_channel_block(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return dose, error, and entry arrays from the individual-volume block."""
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))

    block_index = next(
        index
        for index, row in enumerate(rows)
        if row and row[0].strip().strip("'") == "VOLUME NAMES"
    )
    block = rows[block_index + 5 : block_index + 17]
    if len(block) != 12:
        raise ValueError(f"{path}: expected 12 channel rows, found {len(block)}")

    dose = np.asarray([float(row[0]) for row in block], dtype=float)
    error = np.asarray([float(row[1]) for row in block], dtype=float)
    entries = np.asarray([int(float(row[2])) for row in block], dtype=int)
    return dose, error, entries


def load_energy_grid() -> list[tuple[float, str]]:
    """Load the ordered energy-to-directory mapping."""
    with ENERGY_MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    grid = [(float(row["energy_MeV"]), row["directory"]) for row in rows]
    if grid != sorted(grid):
        raise RuntimeError(f"{ENERGY_MANIFEST}: energies are not sorted")
    if len({energy for energy, _ in grid}) != len(grid):
        raise RuntimeError(f"{ENERGY_MANIFEST}: duplicate energies")
    return grid


def aggregate() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Aggregate ten independent runs at every energy."""
    energy_grid = load_energy_grid()
    means = []
    propagated_errors = []
    empirical_sems = []

    for energy, directory in energy_grid:
        paths = sorted((ROOT / directory / "Res").glob("TID_*.csv"))
        if len(paths) != RUNS_PER_ENERGY:
            raise RuntimeError(
                f"{energy} MeV: expected {RUNS_PER_ENERGY} results, found {len(paths)}"
            )

        doses = []
        errors = []
        for path in paths:
            dose, error, entries = read_channel_block(path)
            if not np.all(entries == PRIMARIES_PER_RUN):
                raise RuntimeError(f"{path}: unexpected channel entry count")
            doses.append(dose)
            errors.append(error)

        dose_array = np.stack(doses)
        error_array = np.stack(errors)
        means.append(np.mean(dose_array, axis=0))
        propagated_errors.append(
            np.sqrt(np.sum(np.square(error_array), axis=0)) / RUNS_PER_ENERGY
        )
        empirical_sems.append(
            np.std(dose_array, axis=0, ddof=1) / math.sqrt(RUNS_PER_ENERGY)
        )

    return (
        np.asarray([energy for energy, _ in energy_grid], dtype=float),
        np.stack(means),
        np.stack(propagated_errors),
        np.stack(empirical_sems),
    )


def write_summary(
    energies: np.ndarray,
    doses: np.ndarray,
    errors: np.ndarray,
    empirical_sems: np.ndarray,
) -> None:
    """Write one tidy row per energy and channel."""
    with SUMMARY_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "energy_MeV",
                "channel",
                "mean_dose_rad",
                "combined_error_rad",
                "relative_error",
                "empirical_sem_rad",
                "runs",
                "total_primaries",
            ]
        )
        for energy_index, energy in enumerate(energies):
            for channel in CHANNELS:
                dose = doses[energy_index, channel]
                error = errors[energy_index, channel]
                writer.writerow(
                    [
                        f"{energy:g}",
                        channel,
                        f"{dose:.10e}",
                        f"{error:.10e}",
                        f"{error / dose:.10e}" if dose > 0 else "nan",
                        f"{empirical_sems[energy_index, channel]:.10e}",
                        RUNS_PER_ENERGY,
                        RUNS_PER_ENERGY * PRIMARIES_PER_RUN,
                    ]
                )


def plot(
    energies: np.ndarray,
    doses: np.ndarray,
    errors: np.ndarray,
) -> None:
    """Plot all twelve channel response kernels with propagated error bars."""
    PLOT_DIR.mkdir(exist_ok=True)

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 13,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )
    figure, axis = plt.subplots(figsize=(12.0, 7.4), constrained_layout=True)

    color_indices = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 1, 3]
    colors = [plt.get_cmap("tab20")(index) for index in color_indices]
    markers = ["o", "s", "^", "D", "v", "P", "X", "<", ">", "h", "*", "d"]
    line_styles = ["-", "-", "-", "-", "--", "--", "--", "--", "-.", "-.", ":", ":"]

    for channel in CHANNELS:
        positive = doses[:, channel] > 0
        axis.errorbar(
            energies[positive],
            doses[positive, channel],
            yerr=errors[positive, channel],
            label=f"Channel {channel}",
            color=colors[channel],
            marker=markers[channel],
            linestyle=line_styles[channel],
            linewidth=1.25,
            markersize=3.6,
            markeredgewidth=0.5,
            elinewidth=0.65,
            capsize=1.5,
            alpha=0.92,
        )

    axis.set_yscale("log")
    axis.set_xlim(float(energies.min()), float(energies.max()))
    axis.set_xlabel("Proton energy at source plane (MeV)")
    axis.set_ylabel(r"Dose at 1 proton cm$^{-2}$ (rad)")
    axis.set_title("RadEx Uppsala final geometry: mono-energy dose response")
    axis.grid(which="major", color="0.82", linewidth=0.7)
    axis.grid(which="minor", color="0.91", linewidth=0.45)
    axis.legend(ncol=4, loc="upper right", frameon=True)
    axis.tick_params(direction="in", which="both")

    figure.savefig(PLOT_STEM.with_suffix(".png"), dpi=300)
    figure.savefig(PLOT_STEM.with_suffix(".pdf"))
    plt.close(figure)


def plot_refined_windows(
    energies: np.ndarray,
    doses: np.ndarray,
    errors: np.ndarray,
) -> None:
    """Plot the response features resolved by the quarter-MeV refinement."""
    refined_stem = PLOT_DIR / "dose-vs-energy-refined-windows"
    figure, axes = plt.subplots(
        1,
        2,
        figsize=(12.0, 5.2),
        constrained_layout=True,
    )

    panels = [
        (axes[0], 63.0, 66.0, (8, 9), "63-66 MeV: channels 8 and 9", True),
        (axes[1], 83.0, 86.0, (11,), "83-86 MeV: channel 11", False),
    ]
    color_indices = {8: 16, 9: 18, 11: 3}
    markers = {8: ">", 9: "h", 11: "d"}

    for axis, lower, upper, channels, title, logarithmic in panels:
        selected = (energies >= lower) & (energies <= upper)
        for channel in channels:
            axis.errorbar(
                energies[selected],
                doses[selected, channel],
                yerr=errors[selected, channel],
                label=f"Channel {channel}",
                color=plt.get_cmap("tab20")(color_indices[channel]),
                marker=markers[channel],
                linestyle="-",
                linewidth=1.35,
                markersize=4.2,
                markeredgewidth=0.5,
                elinewidth=0.7,
                capsize=1.8,
            )

        if logarithmic:
            axis.set_yscale("log")
        axis.set_xlim(lower, upper)
        axis.set_xticks(np.arange(lower, upper + 0.001, 0.5))
        axis.set_xlabel("Proton energy at source plane (MeV)")
        axis.set_title(title)
        axis.grid(which="major", color="0.82", linewidth=0.7)
        axis.grid(which="minor", color="0.91", linewidth=0.45)
        axis.legend(loc="best", frameon=True)
        axis.tick_params(direction="in", which="both")

    axes[0].set_ylabel(r"Dose at 1 proton cm$^{-2}$ (rad)")
    axes[1].set_ylabel(r"Dose at 1 proton cm$^{-2}$ (rad)")
    figure.suptitle("RadEx Uppsala final geometry: refined response windows")
    figure.savefig(refined_stem.with_suffix(".png"), dpi=300)
    figure.savefig(refined_stem.with_suffix(".pdf"))
    plt.close(figure)


def main() -> None:
    energies, doses, errors, empirical_sems = aggregate()
    write_summary(energies, doses, errors, empirical_sems)
    plot(energies, doses, errors)
    plot_refined_windows(energies, doses, errors)

    print(f"Wrote {SUMMARY_PATH}")
    print(f"Wrote {PLOT_STEM.with_suffix('.png')}")
    print(f"Wrote {PLOT_STEM.with_suffix('.pdf')}")
    refined_stem = PLOT_DIR / "dose-vs-energy-refined-windows"
    print(f"Wrote {refined_stem.with_suffix('.png')}")
    print(f"Wrote {refined_stem.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()


