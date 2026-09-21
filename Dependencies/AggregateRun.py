"""Pool one GRAS result directory and validate it before a driver reports it.

Every production driver in this repository repeats the same steps before it can
report anything: count the result files against what the job actually produced,
pool them with Dependencies.TotalDose.totalDose, check that the tally list is
the one the detector macro defines, and check that nothing came back NaN or
Inf. Each of those failures yields a plausible-looking number rather than an
error, so the checks live here once instead of being copied into every driver.

The file count is a real guard, not a formality. GRAS extends the output
filename rather than overwriting when a seed file already exists, so leftovers
from an earlier run appear as chained ``TID_a_b_c_d.csv`` names that are not
independent samples. A job that lost tasks leaves a directory that pools
perfectly well and is simply short of histories. Pass the number of files the
run actually produced, which is not always the array size: GRAS ``autoSeed``
can draw a colliding pair and cost one file, in which case the surviving sample
stays valid and only the history count is short.

Aggregation itself is Dependencies.TotalDose.totalDose and is not reimplemented
here. See its docstring for the entry weighting, the rad-to-kRad conversion and
the pooled-error convention, all of which a fresh parser gets wrong silently.
"""

from pathlib import Path

import numpy as np

from Dependencies.TotalDose import totalDose


def aggregateRun(resultPath, expectedFiles=None, tileCount=None):
    """Pool one GRAS result directory, validating completeness and sanity.

    Args:
        resultPath (str or Path): directory holding the run's GRAS CSV files,
            conventionally ``<run>/Res``.
        expectedFiles (int or None): number of CSV files the run actually
            produced. None skips the completeness check, which is appropriate
            only for exploratory use, never for a reported result.
        tileCount (int or None): number of tallies the detector macro defines.
            None skips the check.

    Returns:
        dict: the keys totalDose returns -- 'dose', 'error', 'entries' and
            'non-zeros', all numpy arrays indexed by tally -- plus

            - 'files': number of CSV files pooled;
            - 'relative_error_percent': per-tally ``100 * error / dose``, zero
              where the dose is zero;
            - 'hit_fraction': per-tally non-zero entries over entries.

    Raises:
        FileNotFoundError: the directory does not exist.
        RuntimeError: the directory is empty, or the file count, the tally
            count or the finiteness check fails.
    """
    resultPath = Path(resultPath)
    if not resultPath.is_dir():
        raise FileNotFoundError(f"result directory not found: {resultPath}")

    files = sorted(resultPath.glob("*.csv"))
    if not files:
        raise RuntimeError(f"no CSV files in {resultPath}")
    if expectedFiles is not None and len(files) != expectedFiles:
        raise RuntimeError(
            f"{resultPath} holds {len(files)} CSV files; expected "
            f"{expectedFiles}. The run is incomplete, or Res/ was not empty "
            "when it started.")

    results = totalDose(str(resultPath))

    dose = np.asarray(results["dose"], dtype=float)
    error = np.asarray(results["error"], dtype=float)
    entries = np.asarray(results["entries"], dtype=float)
    nonZeros = np.asarray(results["non-zeros"], dtype=float)

    if tileCount is not None and len(dose) != tileCount:
        raise RuntimeError(
            f"{resultPath} pooled {len(dose)} tallies; expected {tileCount}. "
            "The detector macro and this driver's tally list disagree.")
    if not np.isfinite(dose).all() or not np.isfinite(error).all():
        raise RuntimeError(f"{resultPath} pooled a non-finite dose or error.")

    results["files"] = len(files)
    results["relative_error_percent"] = np.divide(
        100.0 * error, dose, out=np.zeros_like(error), where=dose != 0)
    results["hit_fraction"] = np.divide(
        nonZeros, entries, out=np.zeros_like(nonZeros), where=entries != 0)
    return results


def tileRecords(results, tiles):
    """Flatten aggregateRun output into one plain dict per tally.

    Args:
        results (dict): what aggregateRun returned.
        tiles (iterable): ``(index, volumeName, label)`` triples in the order
            the detector macro defines them. That order is what fixes which
            tally is which; it is not recoverable from the GRAS output, so the
            caller and the detector macro have to be changed together.

    Returns:
        list[dict]: one record per tally, holding the scalars a driver reports.
    """
    records = []
    for index, volumeName, label in tiles:
        records.append({
            "index": index,
            "volume": volumeName,
            "label": label,
            "dose": float(results["dose"][index]),
            "error": float(results["error"][index]),
            "entries": float(results["entries"][index]),
            "nonZeros": float(results["non-zeros"][index]),
            "relativePercent": float(results["relative_error_percent"][index]),
            "hitFraction": float(results["hit_fraction"][index]),
        })
    return records
