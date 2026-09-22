"""Aggregate the LED displacement-damage tallies of one CRDS run, any campaign.

The DDD sibling of Plotting/CRDSTID.py, and it reads the SAME result files.
CRDS scores both instruments from one set of primaries, so a run directory
holds TID and DDD blocks in the same CSVs; CRDSTID.py reports the RadFET side
and this reports the LED side. Campaigns are defined in Plotting/CRDSCampaigns.py.

    python3 Plotting/CRDSDDD.py --campaign HUS --expected-files 100
    python3 Plotting/CRDSDDD.py --campaign Kumpula --expected-files 100 \
        --folder 10MeVProton-CRDS-100umEpoxy
    python3 Plotting/CRDSDDD.py --campaign HUS --expected-files 20 \
        --job-id 20392554 --folder PDD-Electron-CRDS-100mmField-pilot

Displacement damage is quoted in MeV/g, as GRAS writes it. The coefficient is
per incident particle/cm2 because the run normalises in FLUENCE/CURRENT mode,
so a delivered DDD is the coefficient times the delivered fluence. Kumpula has
no delivered fluence in the model yet; pass ``--reference-fluence`` to quote one.

The two tallies are the 350 nm active layer and the substrate under it. DDD is
intensive, MeV per gram, so a field uniform through the die gives the same
answer in both despite the 500x thickness difference. They are reported side by
side for exactly that reason: at HUS they agree, and at Kumpula, where the
protons stop in or near the die, their ratio is what the epoxy sweep tracks.
"""

import argparse
import csv
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.TotalNID import totalNID, moduleRecords  # noqa: E402
from Plotting.CRDSCampaigns import (CAMPAIGNS, MODULES, campaign,  # noqa: E402
                                    runForFolder)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--campaign", required=True, choices=list(CAMPAIGNS))
    parser.add_argument("--expected-files", type=int, required=True)
    parser.add_argument("--folder",
                        help="run directory; defaults to the campaign's "
                             "production run")
    parser.add_argument("--job-id", type=int)
    parser.add_argument("--reference-fluence", type=float,
                        help="particles/cm2; defaults to the campaign's")
    arguments = parser.parse_args()

    spec = campaign(arguments.campaign)
    folder = arguments.folder or spec["defaultFolder"]
    _, run = runForFolder(spec, folder)
    jobId = arguments.job_id or (run and run["job"])
    if jobId is None:
        parser.error(f"{folder} is not a listed {arguments.campaign} run; "
                     "pass --job-id")
    particle = spec["particle"]
    reference = arguments.reference_fluence or spec["referenceFluence"]

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


if __name__ == "__main__":
    main()
