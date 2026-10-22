"""Compare the pooled tallies of two GRAS runs, tally by tally.

Each run directory holds Res/ (and log/ for the job IDs). The doses come from
Dependencies.TotalDose.totalDose, which entry-weights the per-file doses, pools
the statistical errors and runs the Birge consistency check, so the two runs may
have any number of result files. For every compared tally the tool reports the
candidate-to-baseline ratio, the change in percent and the significance of the
difference against the two statistical errors combined in quadrature.

By default every tally is compared and both runs must score the same number.
--baseline-offset, --candidate-offset and --tiles compare a block instead, e.g.
the silicon dies of a VT01-module run (tallies 12-23) against the twelve bulk
silicon blocks of a legacy run: --candidate-offset 12 --tiles 12.

Run from the repository root:

    python3 -m RadEx.CompareRuns --baseline <run dir> --candidate <run dir> \\
        --output <stem> [--label "what changed"]

writes <stem>.csv and <stem>.pdf (change per tally).
"""

import argparse
import csv
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

PYTHON_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_ROOT))
from Dependencies.TotalDose import totalDose  # noqa: E402
from RadEx import RunSettings  # noqa: E402

FIELDNAMES = [
    "tile_index", "baseline_tile", "candidate_tile",
    "baseline_slurm_job_id", "candidate_slurm_job_id",
    "baseline_primaries", "candidate_primaries",
    "baseline_dose_coefficient", "baseline_statistical_error",
    "candidate_dose_coefficient", "candidate_statistical_error",
    "candidate_to_baseline_ratio", "change_percent", "combined_significance_sigma",
]


def pooled(run):
    """Pooled dose, error and primaries of one run directory."""
    results = totalDose(str(Path(run) / "Res"))
    dose = np.asarray(results["dose"], dtype=float)
    error = np.asarray(results["error"], dtype=float)
    if not np.all(np.isfinite(dose)) or not np.all(np.isfinite(error)):
        raise RuntimeError(f"Non-finite values in {run}")
    jobs = ";".join(map(str, RunSettings.slurmJobs(run)))
    return dose, error, int(results["entries"][0]), jobs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--baseline", type=Path, required=True, help="baseline run directory")
    parser.add_argument("--candidate", type=Path, required=True, help="candidate run directory")
    parser.add_argument("--baseline-offset", type=int, default=0)
    parser.add_argument("--candidate-offset", type=int, default=0)
    parser.add_argument("--tiles", type=int, help="number of tallies compared (default: all)")
    parser.add_argument("--output", type=Path, required=True, help="output stem")
    parser.add_argument("--label", help="plot legend (default: candidate vs baseline directory names)")
    args = parser.parse_args()

    base_dose, base_error, base_primaries, base_jobs = pooled(args.baseline)
    cand_dose, cand_error, cand_primaries, cand_jobs = pooled(args.candidate)
    tiles = args.tiles
    if tiles is None:
        tiles = len(base_dose) - args.baseline_offset
        if len(cand_dose) - args.candidate_offset != tiles:
            parser.error(f"the runs score {len(base_dose)} and {len(cand_dose)} tallies; "
                         "give --tiles and the offsets to compare a block")
    for name, dose, offset in (("baseline", base_dose, args.baseline_offset),
                               ("candidate", cand_dose, args.candidate_offset)):
        if offset < 0 or offset + tiles > len(dose):
            parser.error(f"{name} tallies {offset}-{offset + tiles - 1} exceed its {len(dose)} tallies")
    b = slice(args.baseline_offset, args.baseline_offset + tiles)
    c = slice(args.candidate_offset, args.candidate_offset + tiles)
    base, sbase, cand, scand = base_dose[b], base_error[b], cand_dose[c], cand_error[c]
    if np.any(base <= 0):
        parser.error(f"baseline has non-positive dose in tallies {np.flatnonzero(base <= 0) + args.baseline_offset}")

    ratio = cand / base
    change = 100.0 * (ratio - 1.0)
    combined = np.hypot(sbase, scand)
    significance = (cand - base) / combined

    rows = [{
        "tile_index": t,
        "baseline_tile": args.baseline_offset + t,
        "candidate_tile": args.candidate_offset + t,
        "baseline_slurm_job_id": base_jobs, "candidate_slurm_job_id": cand_jobs,
        "baseline_primaries": base_primaries, "candidate_primaries": cand_primaries,
        "baseline_dose_coefficient": f"{base[t]:.12e}",
        "baseline_statistical_error": f"{sbase[t]:.12e}",
        "candidate_dose_coefficient": f"{cand[t]:.12e}",
        "candidate_statistical_error": f"{scand[t]:.12e}",
        "candidate_to_baseline_ratio": f"{ratio[t]:.9f}",
        "change_percent": f"{change[t]:+.8f}",
        "combined_significance_sigma": f"{significance[t]:+.5f}",
    } for t in range(tiles)]
    output = args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.with_suffix(".csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    label = args.label or f"{args.candidate.resolve().name} vs {args.baseline.resolve().name}"
    x = np.arange(tiles)
    plt.figure()
    plt.errorbar(x, change, yerr=100.0 * combined / base, fmt="o", capsize=4,
                 elinewidth=1, capthick=1, label=label)
    plt.axhline(y=0, color="r", linestyle="--", label="no change")
    plt.xlabel("Tile number")
    plt.ylabel("Change in dose [%]")
    plt.xticks(x)
    plt.grid(which="both")
    plt.legend(fontsize=8)
    plt.savefig(output.with_suffix(".pdf"), format="pdf", bbox_inches="tight")
    plt.close("all")

    print("\nPer-tally comparison (doses in kRad per unit normalisation fluence)")
    print(f"{'Tile':>4s} {'baseline':>14s} {'candidate':>14s} {'change':>9s} {'sigma':>8s}")
    for t in range(tiles):
        print(f"{t:>4d} {base[t]:>14.5e} {cand[t]:>14.5e} {change[t]:>+8.3f}% {significance[t]:>+8.2f}")
    weights = 1.0 / combined ** 2
    worst = int(np.argmax(np.abs(change)))
    print(f"\nError-weighted mean change: {float(np.sum(change * weights) / np.sum(weights)):+.3f} %")
    print(f"Largest change: tile {worst} at {change[worst]:+.3f} %")
    print(f"Tallies differing by more than 3 sigma: {int(np.sum(np.abs(significance) > 3))} of {tiles}")
    print(f"Wrote {output.with_suffix('.csv')} and {output.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
