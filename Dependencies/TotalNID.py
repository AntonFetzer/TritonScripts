"""Pool GRAS non-ionising dose across a result directory.

The displacement-damage counterpart of Dependencies.TotalDose, built on the
same numerical conventions and on Dependencies.PoolDoseModules, whose shape it
follows: several single-volume modules rather than one multi-volume tally.

Conventions carried over deliberately, because they are what a fresh
implementation gets wrong:

* The per-module mean is weighted by ENTRIES, not by file. GRAS writes
  intermediate snapshots into Res/ while a job is still running, so a folder
  can legitimately hold files with different entry counts and an unweighted
  mean over those is wrong.
* The pooled error is sqrt(sum(N_i^2 * sigma_i^2)) / sum(N_i).
* A Birge ratio compares the between-file scatter with the quoted errors. For
  files differing only by their random seed it should sit near 1; well above 2
  means the files differ by more than a seed, which is a bug signal rather
  than extra statistics.

Units are preserved as GRAS wrote them, normally MeV/g. Nothing here rescales
displacement damage.
"""
import os
from math import sqrt

from Read.ReadNID import readNID


def totalNID(path):
    """Pool every NID module in a GRAS result directory.

    Args:
        path (str): Directory of GRAS CSV output files.

    Returns:
        dict: Mapping from module name to a dict with ``nid``, ``error``,
        ``entries``, ``non-zeros``, ``unit``, ``birge_ratio``, ``inconsistent``,
        ``relative_error_percent`` and ``files``.

    Raises:
        ValueError: on a module or unit mismatch between files, non-finite or negative values, or no NID
            module at all.
    """
    files = sorted(f for f in os.listdir(path) if f.endswith(".csv"))
    if not files:
        raise ValueError("No csv files found in folder: " + str(path))

    data = [readNID(os.path.join(path, f)) for f in files]
    data = [d for d in data if d]
    if not data:
        raise ValueError(
            "No NID module found in any file in " + str(path) +
            ". A run scoring only TID writes no NID block."
        )

    names = set(data[0])
    if any(set(d) != names for d in data):
        raise ValueError("NID module names differ between files")

    results = {}
    for name in sorted(names):
        items = [d[name] for d in data]
        units = {d["unit"] for d in items}
        if len(units) != 1:
            raise ValueError("Unit mismatch in NID module " + name)

        entries = [float(d["entries"]) for d in items]
        values = [float(d["nid"]) for d in items]
        errors = [float(d["error"]) for d in items]
        nonZeros = [float(d["non-zeros"]) for d in items]

        for series in (entries, values, errors, nonZeros):
            if any(v != v or v in (float("inf"), float("-inf")) for v in series):
                raise ValueError("Non-finite value in NID module " + name)
        if any(n <= 0 for n in entries):
            raise ValueError("Non-positive entry count in NID module " + name)
        if any(v < 0 for v in values) or any(e < 0 for e in errors):
            raise ValueError("Negative NID or error in module " + name)
        if any(z < 0 or z > n for z, n in zip(nonZeros, entries)):
            raise ValueError("Impossible non-zero entry count in module " + name)

        totalEntries = sum(entries)
        mean = sum(n * v for n, v in zip(entries, values)) / totalEntries
        internal = sum((n * e) ** 2 for n, e in zip(entries, errors))
        external = sum((n * (v - mean)) ** 2 for n, v in zip(entries, values))
        error = sqrt(internal) / totalEntries
        if internal:
            birge = sqrt(external / internal)
        else:
            birge = 0.0 if external == 0 else float("inf")

        relative = 100.0 * error / mean if mean else float("inf")
        results[name] = {
            "nid": mean,
            "error": error,
            "entries": int(totalEntries),
            "non-zeros": int(sum(nonZeros)),
            "unit": units.pop(),
            "birge_ratio": birge,
            "inconsistent": bool(birge > 2),
            "relative_error_percent": relative,
            "files": len(data),
        }

    return results
