"""CRDS lying-versus-upright LED comparison on the first Uppsala board.

One run scores two identical LEDs 10 mm apart in a flat part of the field: one
lying on its back, one standing upright on its leads. The question is whether
the orientation changes the dose in the GaAs, so the report compares the two
parts tally by tally, then turns each coefficient into a delivered value with
the local fluence at its part, as the carrier-board analysis does.

The two dies see different primaries, so their errors are independent and the
upright/lying ratio error is plain quadrature. The upright active layer stands
edge-on to the beam and collects only tens of hits; it is reported, flagged,
and never carries the comparison, which rests on the substrates.

Pooling is Dependencies.AggregateRun and Dependencies.TotalNID; delivered
values use Dependencies.CRDSDelivery. Nothing here re-parses a GRAS file.
"""

import csv

from Dependencies.AggregateRun import aggregateRun, tileRecords
from Dependencies.CRDSCampaigns import MATERIALS, spectrumName
from Dependencies.CRDSDelivery import combine, exposureFluences
from Dependencies.CRDSSeriesAnalysis import ratioWithError
from Dependencies.TotalNID import moduleRecords, totalNID


def pairTallies(spec):
    """Pool the run and return one record per tally, TID then DDD.

    Returns:
        list[dict]: quantity, label, volume, material, part, value, error,
            unit, relative error, entries, non-zeros, hit fraction and, for
            DDD, the Birge ratio; plus whether the tally is edge-on.
    """
    resultPath = spec["path"] / spec["folder"] / "Res"
    particle = spec["particle"]
    tid = aggregateRun(resultPath, expectedFiles=spec["files"],
                       tileCount=len(spec["tiles"]))
    nid = totalNID(str(resultPath), expectedFiles=spec["files"],
                   expectedModules=len(spec["modules"]))

    records = []
    for record in tileRecords(tid, spec["tiles"]):
        material = MATERIALS[record["volume"]]
        records.append({
            "quantity": "TID", "label": record["label"], "volume": record["volume"],
            "material": material, "part": spec["partOfVolume"][record["volume"]],
            "value": record["dose"], "error": record["error"],
            "unit": f"kRad({material})", "coefficientUnit": f"kRad({material}) cm2/{particle}",
            "relativePercent": record["relativePercent"],
            "entries": int(record["entries"]), "nonZeros": int(record["nonZeros"]),
            "hitFraction": record["hitFraction"], "birgeRatio": None,
            "edgeOn": record["volume"] in spec["edgeOn"],
            "index": record["index"],
        })
    for record in moduleRecords(nid, spec["modules"]):
        material = MATERIALS[record["volume"]]
        records.append({
            "quantity": "DDD", "label": record["label"], "volume": record["volume"],
            "material": material, "part": spec["partOfVolume"][record["volume"]],
            "value": record["nid"], "error": record["error"],
            "unit": f"{record['unit']}({material})",
            "coefficientUnit": f"{record['unit']}({material}) cm2/{particle}",
            "relativePercent": record["relativePercent"],
            "entries": int(record["entries"]), "nonZeros": int(record["nonZeros"]),
            "hitFraction": record["hitFraction"], "birgeRatio": record["birgeRatio"],
            "edgeOn": record["volume"] in spec["edgeOn"],
            "module": record["module"], "rawUnit": record["unit"],
            "particlesForOnePercent": record["particlesForOnePercent"],
        })
    return tid["files"], records


def orientationRatios(records):
    """Upright over lying for every tally that exists on both parts."""
    byKey = {(r["quantity"], r["part"], r["label"].split(" ", 1)[1]): r for r in records}
    ratios = []
    for (quantity, part, layer), upright in byKey.items():
        if part != "upright":
            continue
        lying = byKey[(quantity, "lying", layer)]
        ratio, error = ratioWithError(upright["value"], upright["error"],
                                      lying["value"], lying["error"])
        ratios.append({"quantity": quantity, "layer": layer, "ratio": ratio,
                       "error": error, "edgeOn": upright["edgeOn"]})
    return ratios


def deliver(spec, records):
    """Add delivered values from the local fluence at each part, in place.

    The one run supplies the coefficient for both exposures, so its error
    scales with the summed fluence. The alignment range is the spread over a
    board shift of +-alignmentShift_mm in x and in y.
    """
    delivery = spec["delivery"]
    step = delivery["alignmentShift_mm"]
    nominal = exposureFluences(delivery)
    shifted = [exposureFluences(delivery, shift)
               for shift in ((step, 0.0), (-step, 0.0), (0.0, step), (0.0, -step))]
    runOf = ["run"] * len(nominal)
    for record in records:
        coefficient = {"run": (record["value"], record["error"])}
        local = [exposure["local"][record["part"]] for exposure in nominal]
        record["localFluences"] = local
        record["fieldAverageFluences"] = [exposure["fieldAverage"] for exposure in nominal]
        record["exposureDelivered"] = [record["value"] * fluence for fluence in local]
        record["delivered"], record["deliveredError"] = combine(coefficient, runOf, local)
        record["fieldAverageDelivered"], _ = combine(coefficient, runOf,
                                                     record["fieldAverageFluences"])
        aligned = [combine(coefficient, runOf,
                           [exposure["local"][record["part"]] for exposure in fluences])[0]
                   for fluences in shifted]
        record["alignmentLow"], record["alignmentHigh"] = min(aligned), max(aligned)
    return [exposure["label"] for exposure in nominal]


def writeTotalDose(spec, records, source):
    """Write TotalDose_<folder>.csv in the columns of Dependencies.CRDSTIDReport."""
    particle = spec["particle"]
    runPath = spec["path"] / spec["folder"]
    rows = [{
        "particle": particle,
        "source": source,
        "tile_index": r["index"],
        "volume_name": r["volume"],
        "layer": r["label"],
        "slurm_job_id": spec["job"],
        "simulated_primaries": r["entries"],
        "field_size_mm": spec["field"],
        f"dose_coefficient_kRad_cm2_per_{particle}": f"{r['value']:.12e}",
        f"statistical_error_kRad_cm2_per_{particle}": f"{r['error']:.12e}",
        "relative_error_percent": f"{r['relativePercent']:.8f}",
        "nonzero_entries": r["nonZeros"],
        "hit_fraction": f"{r['hitFraction']:.9e}",
        f"reference_fluence_{particle}s_per_cm2": "",
        "dose_at_reference_kRad": "",
        "dose_error_at_reference_kRad": "",
    } for r in records if r["quantity"] == "TID"]
    output = runPath / f"TotalDose_{spec['folder']}.csv"
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output


def writeTotalNID(spec, records):
    """Write TotalNID_<folder>.csv in the columns of Dependencies.CRDSDDDReport."""
    runPath = spec["path"] / spec["folder"]
    rows = [{
        "module": r["module"],
        "volume": r["volume"],
        "label": r["label"],
        "unit": r["rawUnit"],
        "slurm_job_id": spec["job"],
        "simulated_primaries": r["entries"],
        "nid_coefficient_per_particle_cm2": f"{r['value']:.12e}",
        "statistical_error": f"{r['error']:.12e}",
        "relative_error_percent": f"{r['relativePercent']:.8f}",
        "nonzero_entries": r["nonZeros"],
        "hit_fraction": f"{r['hitFraction']:.9e}",
        "birge_ratio": f"{r['birgeRatio']:.6f}",
        "reference_fluence": "",
        "ddd_at_reference": "",
        "particles_for_one_percent": f"{r['particlesForOnePercent']:.6e}",
    } for r in records if r["quantity"] == "DDD"]
    output = runPath / f"TotalNID_{spec['folder']}.csv"
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output


def writeDelivered(spec, records, labels):
    """Write the delivered values in the format of DeliveredDose_CRDS-Uppsala.csv."""
    rows = []
    for r in records:
        common = {"quantity": r["quantity"], "tally": r["label"],
                  "volume_name": r["volume"], "material": r["material"], "part": r["part"]}
        for label, local, average, delivered in zip(labels, r["localFluences"],
                                                    r["fieldAverageFluences"],
                                                    r["exposureDelivered"]):
            rows.append({**common,
                         "exposure": label,
                         "coefficient_run": spec["runLabel"],
                         "coefficient": f"{r['value']:.12e}",
                         "coefficient_error": f"{r['error']:.12e}",
                         "coefficient_unit": r["coefficientUnit"],
                         "local_fluence_per_cm2": f"{local:.6e}",
                         "field_average_fluence_per_cm2": f"{average:.6e}",
                         "delivered": f"{delivered:.6e}",
                         "delivered_error": "", "delivered_field_average": "",
                         "alignment_low": "", "alignment_high": "",
                         "unit": r["unit"]})
        rows.append({**common,
                     "exposure": "total", "coefficient_run": "", "coefficient": "",
                     "coefficient_error": "", "coefficient_unit": "",
                     "local_fluence_per_cm2": "", "field_average_fluence_per_cm2": "",
                     "delivered": f"{r['delivered']:.6e}",
                     "delivered_error": f"{r['deliveredError']:.6e}",
                     "delivered_field_average": f"{r['fieldAverageDelivered']:.6e}",
                     "alignment_low": f"{r['alignmentLow']:.6e}",
                     "alignment_high": f"{r['alignmentHigh']:.6e}",
                     "unit": r["unit"]})
    output = spec["path"] / spec["deliveredOutput"]
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output


def reportLEDPair(spec):
    """Pool, compare, deliver, print and write every CSV; return their paths.

    Writes the comparison (``output``), the delivered values in the
    DeliveredDose format (``deliveredOutput``), and the run folder's
    TotalDose_<folder>.csv and TotalNID_<folder>.csv.
    """
    files, records = pairTallies(spec)
    labels = deliver(spec, records)

    print(f"\n{spec['description']}")
    print(f"Folder       : {spec['folder']}")
    print(f"Slurm job    : {spec['job']}")
    print(f"Result files : {files}")

    print(f"\n{'quantity':<10}{'tally':<24}{'coefficient':>14}{'rel err':>10}"
          f"{'non-zeros':>11}{'Birge':>8}")
    for r in records:
        birge = f"{r['birgeRatio']:.2f}" if r["birgeRatio"] is not None else ""
        flag = "   edge-on, check only" if r["edgeOn"] else ""
        print(f"{r['quantity']:<10}{r['label']:<24}{r['value']:>14.6e}"
              f"{r['relativePercent']:>9.2f}%{r['nonZeros']:>11}{birge:>8}{flag}")
    print("  units: " + ", ".join(sorted({f"{r['quantity']} {r['coefficientUnit']}"
                                          for r in records})))

    print("\nupright / lying, independent dies, quadrature errors")
    for ratio in orientationRatios(records):
        sigma = abs(ratio["ratio"] - 1) / ratio["error"]
        flag = "   edge-on, check only" if ratio["edgeOn"] else ""
        print(f"  {ratio['quantity']:<5}{ratio['layer']:<16}{ratio['ratio']:.4f} +- "
              f"{ratio['error']:.4f}   ({100 * (ratio['ratio'] - 1):+.2f}%, "
              f"{sigma:.1f} sigma){flag}")

    lying = {r["label"]: r for r in records if r["part"] == "lying"}
    ratio, error = ratioWithError(lying["lying active layer"]["value"],
                                  lying["lying active layer"]["error"],
                                  lying["lying substrate"]["value"],
                                  lying["lying substrate"]["error"])
    print(f"\nlying active / lying substrate, DDD: {ratio:.4f} +- {error:.4f}"
          "   (the substrate stands for the active layer if this is 1)")

    print("\ndelivered over " + " and ".join(labels) + ", from the local fluence")
    for r in records:
        flag = "   edge-on, check only" if r["edgeOn"] else ""
        print(f"  {r['quantity']:<5}{r['label']:<24}{r['delivered']:>12.4e} +- "
              f"{r['deliveredError']:.2e} {r['unit']}   field-average "
              f"{r['fieldAverageDelivered']:.4e}   alignment {r['alignmentLow']:.4e}"
              f" to {r['alignmentHigh']:.4e}{flag}")

    output = spec["path"] / spec["output"]
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "quantity", "tally", "volume_name", "material", "part", "slurm_job_id",
            "simulated_primaries", "nonzero_entries", "hit_fraction", "birge_ratio",
            "coefficient_per_proton_cm2", "statistical_error", "coefficient_unit",
            "relative_error_percent",
            *[f"local_fluence_{label}" for label in labels],
            *[f"field_average_fluence_{label}" for label in labels],
            "delivered", "delivered_error", "delivered_field_average",
            "alignment_low", "alignment_high", "unit", "note"])
        for r in records:
            writer.writerow([
                r["quantity"], r["label"], r["volume"], r["material"], r["part"],
                spec["job"], r["entries"], r["nonZeros"], f"{r['hitFraction']:.9e}",
                f"{r['birgeRatio']:.6f}" if r["birgeRatio"] is not None else "",
                f"{r['value']:.12e}", f"{r['error']:.12e}", r["coefficientUnit"],
                f"{r['relativePercent']:.6f}",
                *[f"{value:.6e}" for value in r["localFluences"]],
                *[f"{value:.6e}" for value in r["fieldAverageFluences"]],
                f"{r['delivered']:.6e}", f"{r['deliveredError']:.6e}",
                f"{r['fieldAverageDelivered']:.6e}",
                f"{r['alignmentLow']:.6e}", f"{r['alignmentHigh']:.6e}", r["unit"],
                "edge-on to the beam, check only" if r["edgeOn"] else ""])
    runPath = spec["path"] / spec["folder"]
    return [output,
            writeDelivered(spec, records, labels),
            writeTotalDose(spec, records, spectrumName(runPath)),
            writeTotalNID(spec, records)]
