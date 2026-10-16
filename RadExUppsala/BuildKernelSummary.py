"""Aggregate a mono-energy response-kernel sweep into an InferSpectrum kernel CSV.

Emits the kernel schema that RadExUppsala/InferSpectrum.py and FitAllChannels.py
consume:

    energy_MeV,channel,mean_dose_rad,combined_error_rad,relative_error,
    empirical_sem_rad,runs,total_primaries

Two conventions are load-bearing and easy to get wrong.

Units: the kernel and the measured endpoint doses are in rad, but Read.ReadDose
converts rad to kRad on the way in. The fitted normalisation F in the inference
would absorb a constant factor, but the *inferred fluence* it reports would not,
so the kRad values are converted back to rad here.

Averaging: the documented kernel convention is the arithmetic mean of the N
independently normalised runs with error sqrt(sum(sigma_i^2))/N. Every run of one
energy has the same entry count, so that is identical to the entry-weighted mean
and pooled error that Dependencies.TotalDose.totalDose computes. This script
uses the per-file reader directly because it also needs the empirical
between-run scatter, which totalDose does not return.

The kernel sweep scores 24 volumes: gate oxides in columns 0-11 and silicon dies in
columns 12-23. --channel-offset selects which block becomes channels 0-11.

Run counts: by default each energy is expected to hold exactly the number of
result files given by its `array_tasks` column in energies.csv, each with
`primaries_per_task` entries, so energies topped up with extra replicas aggregate
with their own count. --runs-per-energy and --primaries-per-run force one fixed
value for every energy instead.
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


def _manifest_count(row, column, manifest):
    value = row.get(column)
    if value in (None, ""):
        raise RuntimeError(
            f"{manifest}: no '{column}' for {row['energy_MeV']} MeV; add the column "
            f"or pass a fixed value on the command line")
    count = int(value)
    if count < 1:
        raise RuntimeError(f"{manifest}: non-positive '{column}' for {row['energy_MeV']} MeV")
    return count


def load_run_plan(manifest, runs_per_energy=None, primaries_per_run=None):
    """Return [(energy, directory, runs, primaries_per_run)] from energies.csv.

    A fixed runs_per_energy or primaries_per_run overrides the manifest's
    array_tasks or primaries_per_task column for every energy.
    """
    with manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    plan = []
    for r in rows:
        runs = (runs_per_energy if runs_per_energy is not None
                else _manifest_count(r, "array_tasks", manifest))
        primaries = (primaries_per_run if primaries_per_run is not None
                     else _manifest_count(r, "primaries_per_task", manifest))
        plan.append((float(r["energy_MeV"]), r["directory"], runs, primaries))
    grid = [(energy, directory) for energy, directory, _, _ in plan]
    if grid != sorted(grid):
        raise RuntimeError(f"{manifest}: energies are not sorted")
    if len({e for e, _ in grid}) != len(grid):
        raise RuntimeError(f"{manifest}: duplicate energies")
    return plan


def aggregate(root, offset, runs_per_energy=None, primaries_per_run=None, total_volumes=24):
    plan = load_run_plan(root / "energies.csv", runs_per_energy, primaries_per_run)
    energies, means, errors, sems, runs, primaries = [], [], [], [], [], []

    for energy, directory, n_runs, n_primaries in plan:
        paths = sorted((root / directory / "Res").glob("TID_*.csv"))
        if len(paths) != n_runs:
            raise RuntimeError(
                f"{energy} MeV ({directory}): found {len(paths)} result files, "
                f"expected {n_runs}. A seed collision silently drops a "
                f"replicate when two tasks draw the same autoSeed pair."
            )
        if n_runs < 2:
            raise RuntimeError(f"{energy} MeV ({directory}): the empirical scatter needs at least 2 runs")

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
            if not np.all(entries == n_primaries):
                raise RuntimeError(f"{path}: unexpected entry count")
            doses.append(dose[offset:offset + N_CHANNELS])
            errs.append(err[offset:offset + N_CHANNELS])

        dose_array = np.stack(doses)
        err_array = np.stack(errs)
        energies.append(energy)
        means.append(np.mean(dose_array, axis=0))
        errors.append(np.sqrt(np.sum(np.square(err_array), axis=0)) / n_runs)
        sems.append(np.std(dose_array, axis=0, ddof=1) / math.sqrt(n_runs))
        runs.append(n_runs)
        primaries.append(n_primaries)

    return (np.asarray(energies, dtype=float), np.stack(means),
            np.stack(errors), np.stack(sems),
            np.asarray(runs, dtype=int), np.asarray(primaries, dtype=int))


def write_summary(path, energies, means, errors, sems, runs, primaries_per_run):
    """Write the kernel CSV; runs and primaries_per_run are scalars or per-energy arrays."""
    runs = np.broadcast_to(np.asarray(runs, dtype=int), energies.shape)
    totals = runs * np.broadcast_to(np.asarray(primaries_per_run, dtype=int), energies.shape)
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
                    int(runs[i]), int(totals[i]),
                ])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sweep", type=Path, required=True,
                        help="sweep root containing energies.csv")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--channel-offset", type=int, default=0,
                        help="0 for gate oxides, 12 for silicon dies")
    parser.add_argument("--runs-per-energy", type=int, default=None,
                        help="fixed run count for every energy "
                             "(default: array_tasks per energy from energies.csv)")
    parser.add_argument("--primaries-per-run", type=int, default=None,
                        help="fixed primaries per run for every energy "
                             "(default: primaries_per_task per energy from energies.csv)")
    parser.add_argument("--total-volumes", type=int, default=24)
    args = parser.parse_args()

    energies, means, errors, sems, runs, primaries = aggregate(
        args.sweep, args.channel_offset, args.runs_per_energy,
        args.primaries_per_run, args.total_volumes)

    write_summary(args.output, energies, means, errors, sems, runs, primaries)

    relative = np.divide(errors, means, out=np.zeros_like(errors), where=means > 0)
    strong = means > 0.1 * means.max()
    print(f"\nWrote {args.output}")
    print(f"  energies {len(energies)}, channels {N_CHANNELS}, "
          f"runs per energy {runs.min()}-{runs.max()}")
    print(f"  median relative error on well-populated points: "
          f"{100 * np.median(relative[strong]):.2f} %")
    print(f"  worst relative error on well-populated points:  "
          f"{100 * np.max(relative[strong]):.2f} %")


if __name__ == "__main__":
    main()
