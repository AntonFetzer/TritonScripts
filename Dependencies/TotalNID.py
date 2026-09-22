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


def totalNID(path, expectedFiles=None, expectedModules=None, modulePrefix=None):
    """Pool every NID module in a GRAS result directory.

    Args:
        path (str): Directory of GRAS CSV output files.
        expectedFiles (int or None): If given, refuse a directory that does not
            hold exactly this many CSV files. A short directory means a job is
            still running or tasks failed; a long one means output from another
            run has been mixed in.
        expectedModules (int or None): If given, refuse a file set that does not
            carry exactly this many NID modules.
        modulePrefix (str or None): Optional module-name prefix to retain.

    Returns:
        dict: Mapping from module name to a dict with ``nid``, ``error``,
        ``entries``, ``non-zeros``, ``unit``, ``birge_ratio``, ``inconsistent``,
        ``relative_error_percent``, ``hit_fraction`` and
        ``particles_for_one_percent``.

    Raises:
        ValueError: on a file-count or module-count mismatch, a module or unit
            mismatch between files, non-finite or negative values, or no NID
            module at all.
    """
    files = sorted(f for f in os.listdir(path) if f.endswith(".csv"))
    if not files:
        raise ValueError("No csv files found in folder: " + str(path))
    if expectedFiles is not None and len(files) != expectedFiles:
        raise ValueError(
            "Expected {} csv files in {}, found {}. A short directory means the "
            "job is still running or tasks failed; a long one means another "
            "run's output has been mixed in.".format(expectedFiles, path, len(files))
        )

    data = [readNID(os.path.join(path, f), modulePrefix) for f in files]
    data = [d for d in data if d]
    if not data:
        raise ValueError(
            "No NID module found in any file in " + str(path) +
            ". A run scoring only TID writes no NID block."
        )

    names = set(data[0])
    if any(set(d) != names for d in data):
        raise ValueError("NID module names differ between files")
    if expectedModules is not None and len(names) != expectedModules:
        raise ValueError(
            "Expected {} NID modules, found {}: {}".format(
                expectedModules, len(names), sorted(names))
        )

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
        # How many primaries a 1% result needs, on the 1/sqrt(N) scaling the
        # pooled error follows. Reported rather than applied: a tally with few
        # hits has an unreliable error to extrapolate from in the first place.
        if mean and error:
            needed = (relative / 1.0) ** 2 * totalEntries
        else:
            needed = float("inf")

        results[name] = {
            "nid": mean,
            "error": error,
            "entries": int(totalEntries),
            "non-zeros": int(sum(nonZeros)),
            "unit": units.pop(),
            "birge_ratio": birge,
            "inconsistent": bool(birge > 2),
            "relative_error_percent": relative,
            "hit_fraction": sum(nonZeros) / totalEntries,
            "particles_for_one_percent": needed,
            "files": len(data),
        }

    return results


def moduleRecords(results, modules):
    """Flatten totalNID output into one plain dict per tally, in a fixed order.

    Mirrors Dependencies.AggregateRun.tileRecords.

    Args:
        results (dict): Output of :func:`totalNID`.
        modules (list): ``(moduleName, volumeName, label)`` triples in the order
            the detector macro defines them.

    Returns:
        list: One dict per module, ready for a csv.DictWriter.

    Raises:
        KeyError: if a named module is absent from the results.
    """
    records = []
    for moduleName, volumeName, label in modules:
        if moduleName not in results:
            raise KeyError(
                "NID module {!r} not in results; present: {}".format(
                    moduleName, sorted(results))
            )
        entry = results[moduleName]
        records.append({
            "module": moduleName,
            "volume": volumeName,
            "label": label,
            "nid": entry["nid"],
            "error": entry["error"],
            "unit": entry["unit"],
            "entries": entry["entries"],
            "nonZeros": entry["non-zeros"],
            "relativePercent": entry["relative_error_percent"],
            "hitFraction": entry["hit_fraction"],
            "particlesForOnePercent": entry["particles_for_one_percent"],
            "birgeRatio": entry["birge_ratio"],
        })
    return records