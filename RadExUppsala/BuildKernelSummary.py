"""Aggregate a mono-energy response-kernel sweep into an InferSpectrum kernel CSV.

Emits the same schema as EnergyKernel-FinalGeometry/energy-kernel-dose-summary.csv,
which is what RadExUppsala/InferSpectrum.py consumes:

    energy_MeV,channel,mean_dose_rad,combined_error_rad,relative_error,
    empirical_sem_rad,runs,total_primaries

Two conventions are load-bearing and easy to get wrong.

Units: the kernel and the measured endpoint doses are in rad, but Read.ReadDose
converts rad to kRad on the way in. The fitted normalisation F in the inference
would absorb a constant factor, but the *inferred fluence* it reports would not,
so the kRad values are converted back to rad here.

Averaging: the documented kernel convention is the arithmetic mean of the N
independently normalised runs with error sqrt(sum(sigma_i^2))/N. Every run in a
sweep has the same entry count, so that is identical to the entry-weighted mean
and pooled error that Dependencies.TotalDose.totalDose computes. This script
uses the per-file reader directly because it also needs the empirical
between-run scatter, which totalDose does not return.

The new sweep scores 24 volumes: gate oxides in columns 0-11 and silicon dies in
columns 12-23. --channel-offset selects which block becomes channels 0-11.
"""

import argparse
import csv
import math
from pathlib import Path
import sys

import numpy as np

PYTHON_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_ROOT))

from Read.ReadDose import readDose  # noqa: E402

KRAD_TO_RAD = 1000.0
N_CHANNELS = 12
FIELDNAMES = [
    "energy_MeV", "channel", "mean_dose_rad", "combined_error_rad",
    "relative_error", "empirical_sem_rad", "runs", "total_primaries",
]


def load_energy_grid(manifest):
    with manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    grid = [(float(r["energy_MeV"]), r["directory"]) for r in rows]
    if grid != sorted(grid):
        raise RuntimeError(f"{manifest}: energies are not sorted")
    if len({e for e, _ in grid}) != len(grid):
        raise RuntimeError(f"{manifest}: duplicate energies")
    return grid


def aggregate(root, offset, runs_per_energy, primaries_per_run, total_volumes):
    grid = load_energy_grid(root / "energies.csv")
    energies, means, errors, sems = [], [], [], []

    for energy, directory in grid:
        paths = sorted((root / directory / "Res").glob("TID_*.csv"))
        if len(paths) != runs_per_energy:
            raise RuntimeError(
                f"{energy} MeV ({directory}): found {len(paths)} result files, "
                f"expected {runs_per_energy}. A seed collision silently drops a "
                f"replicate when two tasks draw the same autoSeed pair."
            )

        doses, errs = [], []
        for path in paths:
            raw = readDose(str(path))
            dose = np.asarray(raw["dose"], dtype=float) * KRAD_TO_RAD
            err = np.asarray(raw["error"], dtype=float) * KRAD_TO_RAD
            entries = np.asarray(raw["entries"], dtype=float)
            if len(dose) != total_volumes:
                raise RuntimeError(
                    f"{path}: {len(dose)} volumes, expected {total_volumes}"
                )
            if not np.all(entries == primaries_per_run):
                raise RuntimeError(f"{path}: unexpected entry count")
            doses.append(dose[offset:offset + N_CHANNELS])
            errs.append(err[offset:offset + N_CHANNELS])

        dose_array = np.stack(doses)
        err_array = np.stack(errs)
        energies.append(energy)
        means.append(np.mean(dose_array, axis=0))
        errors.append(np.sqrt(np.sum(np.square(err_array), axis=0)) / runs_per_energy)
        sems.append(np.std(dose_array, axis=0, ddof=1) / math.sqrt(runs_per_energy))

    return (np.asarray(energies, dtype=float), np.stack(means),
            np.stack(errors), np.stack(sems))


def write_summary(path, energies, means, errors, sems, runs, primaries_per_run):
    total = runs * primaries_per_run
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(FIELDNAMES)
        for i, energy in enumerate(energies):
            for channel in range(N_CHANNELS):
                mean = means[i, channel]
                error = errors[i, channel]
                relative = error / mean if mean > 0 else 0.0
                writer.writerow([
                    f"{energy:g}", channel,
                    f"{mean:.10e}", f"{error:.10e}",
                    f"{relative:.10e}", f"{sems[i, channel]:.10e}",
                    runs, total,
                ])


def validate_against(reference, energies, means, errors):
    """Numerically compare against an existing summary produced by the old script."""
    with reference.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    lookup = {(float(r["energy_MeV"]), int(r["channel"])):
              (float(r["mean_dose_rad"]), float(r["combined_error_rad"]))
              for r in rows}

    worst_dose = 0.0
    worst_error = 0.0
    compared = 0
    for i, energy in enumerate(energies):
        for channel in range(N_CHANNELS):
            key = (float(energy), channel)
            if key not in lookup:
                raise RuntimeError(f"reference has no entry for {key}")
            ref_dose, ref_error = lookup[key]
            if ref_dose > 0:
                worst_dose = max(worst_dose, abs(means[i, channel] / ref_dose - 1.0))
            if ref_error > 0:
                worst_error = max(worst_error, abs(errors[i, channel] / ref_error - 1.0))
            compared += 1
    print(f"Validation against {reference.name}: {compared} points compared")
    print(f"  worst relative difference in mean dose: {worst_dose:.3e}")
    print(f"  worst relative difference in error:     {worst_error:.3e}")
    return worst_dose, worst_error


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sweep", type=Path, required=True,
                        help="sweep root containing energies.csv")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--channel-offset", type=int, default=0,
                        help="0 for gate oxides, 12 for silicon dies")
    parser.add_argument("--runs-per-energy", type=int, default=20)
    parser.add_argument("--primaries-per-run", type=int, default=4_000_000)
    parser.add_argument("--total-volumes", type=int, default=24)
    parser.add_argument("--validate-against", type=Path, default=None)
    args = parser.parse_args()

    energies, means, errors, sems = aggregate(
        args.sweep, args.channel_offset, args.runs_per_energy,
        args.primaries_per_run, args.total_volumes)

    if args.validate_against is not None:
        validate_against(args.validate_against, energies, means, errors)

    write_summary(args.output, energies, means, errors, sems,
                  args.runs_per_energy, args.primaries_per_run)

    relative = np.divide(errors, means, out=np.zeros_like(errors), where=means > 0)
    strong = means > 0.1 * means.max()
    print(f"\nWrote {args.output}")
    print(f"  energies {len(energies)}, channels {N_CHANNELS}")
    print(f"  median relative error on well-populated points: "
          f"{100 * np.median(relative[strong]):.2f} %")
    print(f"  worst relative error on well-populated points:  "
          f"{100 * np.max(relative[strong]):.2f} %")


if __name__ == "__main__":
    main()
