"""Aggregate the RadFET TID tallies of one CRDS production run, any campaign.

One VT01 RadFET package on the 1.6 mm carrier PCB. Two tallies, the 400 nm
gate oxide and the silicon die behind it, in the order
``CRDS/CRDSRadFETLEDDetector.mac`` fixes. The campaign, its runs and its
reference fluence are defined in Plotting/CRDSCampaigns.py.

    python3 Plotting/CRDSTID.py --campaign HUS --expected-files 100
    python3 Plotting/CRDSTID.py --campaign Kumpula --expected-files 100 \
        --folder 10MeVProton-CRDS-100umEpoxy

``--job-id`` defaults to the job recorded for a listed run and is required for
an unlisted folder such as a pilot.

HUS doses are quoted at 2e12 electrons/cm2, the assumed reference fluence the
RadEx-HUS results use; that is not an MU calibration and not a measured
delivered dose. The 20 x 20 mm HUS field is below lateral scatter equilibrium,
so see Plotting/CRDSSeries.py before quoting it for a broad-field exposure.

Kumpula has no delivered fluence in the model yet, so its results stay per
proton/cm2 unless ``--reference-fluence`` supplies one.
"""

import argparse
import csv
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.AggregateRun import aggregateRun, tileRecords  # noqa: E402
from Plotting.CRDSCampaigns import (CAMPAIGNS, TILES, campaign,  # noqa: E402
                                    runForFolder, spectrumName)


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

    results = aggregateRun(
        runPath / "Res",
        expectedFiles=arguments.expected_files,
        tileCount=len(TILES),
    )

    print(f"\nCampaign         : {arguments.campaign}, {spec['description']}")
    print(f"Folder           : {folder}")
    print(f"Result files     : {results['files']}")
    if reference:
        note = ("supplied on the command line" if arguments.reference_fluence
                else spec["referenceNote"])
        print(f"Reference fluence: {reference:.3e} {particle}s/cm2 ({note})\n")
    else:
        print(f"Reference fluence: none, coefficients per {particle}/cm2 only\n")

    rows = []
    records = tileRecords(results, TILES)
    for record in records:
        dose, error = record["dose"], record["error"]
        print(f"{record['label']} ({record['volume']})")
        print(f"  pooled primaries    : {record['entries']:.0f}")
        print(f"  non-zero entries    : {record['nonZeros']:.0f}")
        print(f"  hit fraction        : {record['hitFraction']:.6e}")
        print(f"  dose coefficient    : {dose:.6e} kRad cm2/{particle}")
        print(f"  statistical error   : {error:.6e} kRad cm2/{particle}")
        print(f"  relative error      : {record['relativePercent']:.4f}%")
        if reference:
            print(f"  dose at reference   : {dose * reference:.4f} "
                  f"+- {error * reference:.4f} kRad")
        print()

        rows.append({
            "particle": particle,
            "source": spectrumName(runPath),
            "tile_index": record["index"],
            "volume_name": record["volume"],
            "layer": record["label"],
            "slurm_job_id": jobId,
            "simulated_primaries": int(record["entries"]),
            "field_size_mm": run["field"] if run else "",
            f"dose_coefficient_kRad_cm2_per_{particle}": f"{dose:.12e}",
            f"statistical_error_kRad_cm2_per_{particle}": f"{error:.12e}",
            "relative_error_percent": f"{record['relativePercent']:.8f}",
            "nonzero_entries": int(record["nonZeros"]),
            "hit_fraction": f"{record['hitFraction']:.9e}",
            f"reference_fluence_{particle}s_per_cm2":
                f"{reference:.12e}" if reference else "",
            "dose_at_reference_kRad":
                f"{dose * reference:.12e}" if reference else "",
            "dose_error_at_reference_kRad":
                f"{error * reference:.12e}" if reference else "",
        })

    oxideDose, dieDose = records[0]["dose"], records[1]["dose"]
    print(f"oxide / die                  : {oxideDose / dieDose:.4f}")
    comparison = spec["oxideComparison"]
    if comparison and reference:
        print(f"oxide vs {comparison['label']} : "
              f"{100.0 * (oxideDose * reference / comparison['kRad'] - 1.0):+.2f} % "
              f"({comparison['note']})")

    output_path = runPath / f"TotalDose_{folder}.csv"
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {output_path}")


if __name__ == "__main__":
    main()
