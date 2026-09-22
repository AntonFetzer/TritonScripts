"""CRDS series aggregation, uncertainty calculations and CSV reporting.

GRAS pooling stays in Dependencies.AggregateRun and Dependencies.TotalNID.
Within-run ratio errors use the documented independent quadrature upper bound.
"""

import csv
import math

from Dependencies.AggregateRun import aggregateRun
from Dependencies.TotalNID import totalNID
from Dependencies.CRDSCampaigns import MATERIALS, MODULES, TILES
from Dependencies.CRDSDelivery import deliveredRecords


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


def reportDelivery(spec, data, notCoefficients):
    """Print and write the delivered value of every tally; return the CSV path.

    Coefficients flagged in ``notCoefficients`` are never offered to the
    combination, so a campaign cannot silently route a part through a run in
    which it was outside the simulated field.
    """
    coefficients = {}
    for label, entries in data.items():
        for (kind, name), (value, error, *_rest) in entries.items():
            if (label, kind) not in notCoefficients:
                coefficients.setdefault((kind, name), {})[label] = (value, error)
    records = deliveredRecords(spec["delivery"], coefficients)
    units = {key: entry[2] for entries in data.values() for key, entry in entries.items()}
    volumes = {key: entry[3] for entries in data.values() for key, entry in entries.items()}

    shift = spec["delivery"]["alignmentShift_mm"]
    print(f"\ndelivered, from the local fluence at each part "
          f"(alignment range: board shifted +-{shift:g} mm)")
    rows = []
    for record in records:
        key = (record["quantity"], record["tally"])
        unit = units[key]
        delivered, error = record["delivered"], record["error"]
        print(f"  {record['quantity']:<4}{record['tally']:<16}{delivered:>12.4e}"
              f" +- {error:.2e} {unit}   field-average {record['fieldAverageDelivered']:.4e}"
              f" ({100 * (record['fieldAverageDelivered'] / delivered - 1):+.1f}%)"
              f"   alignment {record['alignmentLow']:.4e} to {record['alignmentHigh']:.4e}")
        common = {"quantity": record["quantity"], "tally": record["tally"],
                  "volume_name": volumes[key], "material": MATERIALS[volumes[key]],
                  "part": record["part"]}
        for exposure in record["exposures"]:
            rows.append({**common,
                         "exposure": exposure["label"],
                         "coefficient_run": exposure["run"],
                         "coefficient": f"{exposure['coefficient']:.12e}",
                         "coefficient_error": f"{exposure['coefficientError']:.12e}",
                         "coefficient_unit": f"{unit} cm2/{spec['particle']}",
                         "local_fluence_per_cm2": f"{exposure['localFluence']:.6e}",
                         "field_average_fluence_per_cm2": f"{exposure['fieldAverageFluence']:.6e}",
                         "delivered": f"{exposure['delivered']:.6e}",
                         "delivered_error": "", "delivered_field_average": "",
                         "alignment_low": "", "alignment_high": "",
                         "unit": unit})
        rows.append({**common,
                     "exposure": "total", "coefficient_run": "", "coefficient": "",
                     "coefficient_error": "", "coefficient_unit": "",
                     "local_fluence_per_cm2": "", "field_average_fluence_per_cm2": "",
                     "delivered": f"{delivered:.6e}",
                     "delivered_error": f"{error:.6e}",
                     "delivered_field_average": f"{record['fieldAverageDelivered']:.6e}",
                     "alignment_low": f"{record['alignmentLow']:.6e}",
                     "alignment_high": f"{record['alignmentHigh']:.6e}",
                     "unit": unit})

    output = spec["path"] / spec["delivery"]["output"]
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output


def collectSeries(spec, reference):
    """Pool campaign runs and return keyed tallies plus established CSV rows."""
    particle = spec["particle"]
    runs = spec["runs"]
    notCoefficients = spec.get("notCoefficients", set())
    notCoefficientNote = "not a coefficient: part outside the simulated field"
    data, rows = {}, []
    for label, run in runs.items():
        resultPath = spec["path"] / run["folder"] / "Res"
        tid = aggregateRun(str(resultPath), expectedFiles=run["files"],
                           tileCount=len(TILES))
        ddd = totalNID(str(resultPath), expectedFiles=run["files"],
                       expectedModules=len(MODULES))
        # Units name the tallied material: GRAS scores per gram of the
        # volume's own material, and Read.ReadDose converts rad to kRad.
        entries = {}
        for index, volume, name in TILES:
            entries[("TID", name)] = (float(tid["dose"][index]),
                                      float(tid["error"][index]),
                                      f"kRad({MATERIALS[volume]})", volume,
                                      int(tid["entries"][index]),
                                      int(tid["non-zeros"][index]))
        for module, volume, name in MODULES:
            entry = ddd[module]
            entries[("DDD", name)] = (entry["nid"], entry["error"],
                                      f"{entry['unit']}({MATERIALS[volume]})", volume,
                                      entry["entries"], entry["non-zeros"])
        data[label] = entries

        for (kind, name), (value, error, unit, volume, n, nz) in entries.items():
            row = {"run": label}
            if spec["sweep"]:
                row[spec["sweep"]["key"]] = run[spec["sweep"]["key"]]
            row.update({
                "folder": run["folder"], "slurm_job_id": run["job"],
                "quantity": kind, "tally": name, "volume_name": volume,
                "material": MATERIALS[volume],
                "simulated_primaries": n,
                "nonzero_entries": nz,
                f"coefficient_per_{particle}_cm2": f"{value:.12e}",
                "statistical_error": f"{error:.12e}",
                "unit": f"{unit} cm2/{particle}",
                "relative_error_percent": f"{100 * error / value:.6f}",
                "reference_fluence": f"{reference:.6e}" if reference else "",
                "at_reference_fluence":
                    f"{value * reference:.6e}" if reference else "",
                "at_reference_unit": unit if reference else "",
            })
            if notCoefficients:
                row["note"] = notCoefficientNote if (label, kind) in notCoefficients else ""
            rows.append(row)

    return data, rows


def printSeries(name, spec, data, reference):
    """Print the campaign comparisons without changing coefficient units."""
    runs = spec["runs"]
    notCoefficients = spec.get("notCoefficients", set())
    width = max(len(k) for k in runs) + 2
    atHeader = f"at {reference:.3g}" if reference else ""
    print(f"\n{name}: {spec['description']}")
    print(f"\n{'run':<{width}}{'quantity':<10}{'tally':<16}{'coefficient':>16}"
          f"{'rel err':>10}{atHeader:>16}")
    for label in runs:
        for (kind, name), (value, error, unit, _, _, _) in data[label].items():
            shown = ""
            if reference:
                at = value * reference
                shown = f"{at:.3f} kRad" if kind == "TID" else f"{at:.4g} MeV/g"
            flag = "   (not a coefficient)" if (label, kind) in notCoefficients else ""
            print(f"{label:<{width}}{kind:<10}{name:<16}{value:>16.6e}"
                  f"{100 * error / value:>9.2f}%{shown:>16}{flag}")

    # Within-run ratios. Both tallies of a pair see the same primaries, so
    # their errors are correlated and this quadrature error is an upper bound.
    print(f"\n{'run':<{width}}{'oxide / die':>20}{'active / substrate':>24}")
    for label in runs:
        pairs = [(data[label][("TID", TILES[0][2])], data[label][("TID", TILES[1][2])]),
                 (data[label][("DDD", MODULES[0][2])], data[label][("DDD", MODULES[1][2])])]
        shown = []
        for kind, (a, b) in zip(("TID", "DDD"), pairs):
            if (label, kind) in notCoefficients:
                shown.append("n/a")
                continue
            ratio, error = ratioWithError(a[0], a[1], b[0], b[1])
            shown.append(f"{ratio:.4f} +- {error:.4f}")
        print(f"{label:<{width}}{shown[0]:>20}{shown[1]:>24}")

    print("\nsteps, each against the run named second")
    for numerator, denominator in spec["steps"]:
        print(f"  {numerator}  vs  {denominator}")
        for key in data[numerator]:
            if {(numerator, key[0]), (denominator, key[0])} & notCoefficients:
                print(f"    {key[0]:<4}{key[1]:<16}  skipped, not a coefficient in both runs")
                continue
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



def writeSeries(spec, rows):
    """Write the established campaign coefficient CSV."""
    output = spec["path"] / spec["seriesOutput"]
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output

