"""Aggregate the displacement-damage tallies of a CRDS HUS production run.

The DDD sibling of Plotting/CRDSHUSTID.py, and it reads the SAME result files.
CRDS scores both instruments from one set of primaries, so a run directory
holds TID and DDD blocks in the same CSVs; CRDSHUSTID.py reports the RadFET
side and this reports the LED side.

    python3 Plotting/CRDSHUSDDD.py --expected-files 100 --job-id 20392554
    python3 Plotting/CRDSHUSDDD.py --expected-files 20 --job-id 20392554 \
        --folder PDD-Electron-CRDS-100mmField-pilot

Displacement damage is quoted in MeV/g, as GRAS writes it. The coefficient is
per incident particle/cm2 because the run normalises in FLUENCE/CURRENT mode,
so a delivered DDD is the coefficient times the delivered fluence.

The two tallies are the 350 nm active layer and the substrate under it. DDD is
intensive, MeV per gram, so a field uniform through the die gives the same
answer in both despite the 500x thickness difference. They are reported side by
side for exactly that reason: a divergence means the field has stopped being
uniform across the die and the depth of the active layer has begun to matter.
"""

import argparse
import csv
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.TotalNID import totalNID, moduleRecords  # noqa: E402

BASE_PATH = Path("/scratch/work/fetzera1/GRAS/CRDS/CRDS-HUS")
DEFAULT_FOLDER = "PDD-Electron-CRDS-20mmField"

# Module order is fixed by CRDS/CRDSRadFETLEDDetector.mac. Change both together.
MODULES = [
    ("nidActive", "LED_active_0_PV", "350 nm active layer"),
    ("nidSubstrate", "LED_die_0_PV", "GaAs substrate"),
]

# Campaign reference normalisation, the same assumed fluence the TID results
# are quoted at. This is an assumed fluence, not an MU calibration.
REFERENCE_FLUENCE = 2e12


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--expected-files", type=int, required=True)
    parser.add_argument("--job-id", type=int, required=True)
    parser.add_argument("--folder", default=DEFAULT_FOLDER)
    parser.add_argument("--reference-fluence", type=float, default=REFERENCE_FLUENCE)
    arguments = parser.parse_args()

    runPath = BASE_PATH / arguments.folder
    results = totalNID(
        str(runPath / "Res"),
        expectedFiles=arguments.expected_files,
        expectedModules=len(MODULES),
    )
    records = moduleRecords(results, MODULES)

    print()
    print("Folder           : {}".format(arguments.folder))
    print("Slurm job        : {}".format(arguments.job_id))
    print("Result files     : {}".format(results[MODULES[0][0]]["files"]))
    print("Reference fluence: {:.3e} particles/cm2 (assumed)".format(
        arguments.reference_fluence))
    print()

    for record in records:
        delivered = record["nid"] * arguments.reference_fluence
        deliveredError = record["error"] * arguments.reference_fluence
        print("{} ({})".format(record["label"], record["volume"]))
        print("  pooled primaries    : {}".format(record["entries"]))
        print("  non-zero entries    : {}".format(record["nonZeros"]))
        print("  hit fraction        : {:.6e}".format(record["hitFraction"]))
        print("  DDD coefficient     : {:.6e} {} per particle/cm2".format(
            record["nid"], record["unit"]))
        print("  statistical error   : {:.6e} {}".format(record["error"], record["unit"]))
        print("  relative error      : {:.4f}%".format(record["relativePercent"]))
        print("  Birge ratio         : {:.2f}{}".format(
            record["birgeRatio"], "  INCONSISTENT" if record["birgeRatio"] > 2 else ""))
        print("  DDD at reference    : {:.4f} +- {:.4f} {}".format(
            delivered, deliveredError, record["unit"]))
        print("  primaries for 1%    : {:.3e}  ({:.3g}x this run)".format(
            record["particlesForOnePercent"],
            record["particlesForOnePercent"] / record["entries"]))
        print()

    active, substrate = records[0], records[1]
    if substrate["nid"]:
        ratio = active["nid"] / substrate["nid"]
        relative = (active["relativePercent"] ** 2 + substrate["relativePercent"] ** 2) ** 0.5
        print("active / substrate  : {:.4f} +- {:.4f}".format(ratio, ratio * relative / 100))
        print("  These should agree: DDD is intensive and the field is uniform")
        print("  across 0.18 mm of GaAs. A divergence means it no longer is.")
        print()

    outputPath = runPath / "TotalNID_{}.csv".format(arguments.folder)
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
                "slurm_job_id": arguments.job_id,
                "simulated_primaries": record["entries"],
                "nid_coefficient_per_particle_cm2": "{:.12e}".format(record["nid"]),
                "statistical_error": "{:.12e}".format(record["error"]),
                "relative_error_percent": "{:.8f}".format(record["relativePercent"]),
                "nonzero_entries": record["nonZeros"],
                "hit_fraction": "{:.9e}".format(record["hitFraction"]),
                "birge_ratio": "{:.6f}".format(record["birgeRatio"]),
                "reference_fluence": "{:.6e}".format(arguments.reference_fluence),
                "ddd_at_reference": "{:.6e}".format(
                    record["nid"] * arguments.reference_fluence),
                "particles_for_one_percent": "{:.6e}".format(
                    record["particlesForOnePercent"]),
            })
    print("Wrote {}".format(outputPath))


if __name__ == "__main__":
    main()