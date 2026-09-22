"""CRDS DDD aggregation and reporting for one GRAS run."""

import csv

from Dependencies.TotalNID import totalNID, moduleRecords
from Dependencies.CRDSCampaigns import MODULES, campaign, runForFolder


def reportDDD(arguments, parser):
    """Aggregate one run and write its established CSV report."""
    spec = campaign(arguments.campaign)
    folder = arguments.folder or spec["defaultFolder"]
    label, run = runForFolder(spec, folder)
    jobId = arguments.job_id or (run and run["job"])
    if jobId is None:
        parser.error(f"no Slurm job recorded for {arguments.campaign} run "
                     f"{folder}; pass --job-id")
    particle = spec["particle"]
    reference = arguments.reference_fluence or spec["referenceFluence"]
    notCoefficient = (label, "DDD") in spec.get("notCoefficients", set())

    runPath = spec["path"] / folder
    results = totalNID(
        str(runPath / "Res"),
        expectedFiles=arguments.expected_files,
        expectedModules=len(MODULES),
    )
    records = moduleRecords(results, MODULES)

    print()
    print("Campaign         : {}, {}".format(arguments.campaign, spec["description"]))
    print("Folder           : {}".format(folder))
    print("Slurm job        : {}".format(jobId))
    print("Result files     : {}".format(results[MODULES[0][0]]["files"]))
    if reference:
        print("Reference fluence: {:.3e} {}s/cm2 ({})".format(
            reference, particle,
            "supplied on the command line" if arguments.reference_fluence
            else spec["referenceNote"]))
    else:
        print("Reference fluence: none, coefficients per {}/cm2 only".format(particle))
    if notCoefficient:
        print("WARNING          : the LED is outside this run's simulated field;")
        print("                   its DDD here is NOT a usable coefficient.")
    print()

    for record in records:
        print("{} ({})".format(record["label"], record["volume"]))
        print("  pooled primaries    : {}".format(record["entries"]))
        print("  non-zero entries    : {}".format(record["nonZeros"]))
        print("  hit fraction        : {:.6e}".format(record["hitFraction"]))
        print("  DDD coefficient     : {:.6e} {} per {}/cm2".format(
            record["nid"], record["unit"], particle))
        print("  statistical error   : {:.6e} {}".format(record["error"], record["unit"]))
        print("  relative error      : {:.4f}%".format(record["relativePercent"]))
        print("  Birge ratio         : {:.2f}{}".format(
            record["birgeRatio"], "  INCONSISTENT" if record["birgeRatio"] > 2 else ""))
        if reference:
            print("  DDD at reference    : {:.4f} +- {:.4f} {}".format(
                record["nid"] * reference, record["error"] * reference,
                record["unit"]))
        print("  primaries for 1%    : {:.3e}  ({:.3g}x this run)".format(
            record["particlesForOnePercent"],
            record["particlesForOnePercent"] / record["entries"]))
        print()

    active, substrate = records[0], records[1]
    if substrate["nid"]:
        ratio = active["nid"] / substrate["nid"]
        relative = (active["relativePercent"] ** 2 + substrate["relativePercent"] ** 2) ** 0.5
        print("active / substrate  : {:.4f} +- {:.4f}".format(ratio, ratio * relative / 100))
        print("  " + spec["activeSubstrateNote"])
        print()

    outputPath = runPath / "TotalNID_{}.csv".format(folder)
    fieldnames = ["module", "volume", "label", "unit", "slurm_job_id",
                  "simulated_primaries", "nid_coefficient_per_particle_cm2",
                  "statistical_error", "relative_error_percent", "nonzero_entries",
                  "hit_fraction", "birge_ratio", "reference_fluence",
                  "ddd_at_reference", "particles_for_one_percent"]
    with outputPath.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            writer.writerow({
                "module": record["module"],
                "volume": record["volume"],
                "label": record["label"],
                "unit": record["unit"],
                "slurm_job_id": jobId,
                "simulated_primaries": record["entries"],
                "nid_coefficient_per_particle_cm2": "{:.12e}".format(record["nid"]),
                "statistical_error": "{:.12e}".format(record["error"]),
                "relative_error_percent": "{:.8f}".format(record["relativePercent"]),
                "nonzero_entries": record["nonZeros"],
                "hit_fraction": "{:.9e}".format(record["hitFraction"]),
                "birge_ratio": "{:.6f}".format(record["birgeRatio"]),
                "reference_fluence":
                    "{:.6e}".format(reference) if reference else "",
                "ddd_at_reference":
                    "{:.6e}".format(record["nid"] * reference) if reference else "",
                "particles_for_one_percent": "{:.6e}".format(
                    record["particlesForOnePercent"]),
            })
    print("Wrote {}".format(outputPath))

