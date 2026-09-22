"""Aggregate the CRDS HUS PDD electron production run.

CRDS in the HUS TrueBeam electron field: one VT01 RadFET package on a 1.6 mm
carrier PCB, on 20 mm standoffs above the plastic patient table, in the
20 x 20 mm PDD field. Two tallies, the 400 nm gate oxide and the silicon die
behind it, in the order ``CRDS/CRDS1RadFETDetector.mac`` fixes.

Sibling of CRDSUppsalaTID.py and of RadExHUSPDDProduction.py.

    python3 Plotting/CRDSHUSTID.py --expected-files 100 --job-id 20332204

The dose is quoted at 2e12 electrons/cm2, the assumed reference fluence the
RadEx-HUS results use. That is an assumed normalisation, not an MU calibration
and not a measured delivered dose.

The 20 x 20 mm field is below lateral scatter equilibrium and understates a
broad-field exposure by about 16% on the die. If the real exposure used the
clinical applicator, the wider runs are the ones to quote; see
Plotting/CRDSHUSFieldSize.py.
"""

import argparse
import csv
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.AggregateRun import aggregateRun, tileRecords  # noqa: E402

BASE_PATH = Path("/scratch/work/fetzera1/GRAS/CRDS/CRDS-HUS")
DEFAULT_FOLDER = "PDD-Electron-CRDS-20mmField"

# Tally order is fixed by CRDS/CRDS1RadFETDetector.mac. Change both together.
TILES = [
    (0, "VT01_gox_0_PV", "gate oxide"),
    (1, "VT01_die_0_PV", "silicon die"),
]

# Campaign reference normalisation, the same assumed fluence the RadEx-HUS
# results are quoted at. This is an assumed fluence, not an MU calibration.
REFERENCE_FLUENCE = 2e12
SOURCE_AREA_CM2 = 2.0 * 2.0

# RadEx-HUS Ch9, the exposed 0 mm channel, is the nearest thing to a bare
# device in the same field with the same detector model: job 20315277 gave
# 61.435 kRad on its gate oxide at the same reference fluence. CRDS sits at the
# same 200 mm air distance but has no aluminium top plate, no spacer and a
# 1.6 mm board instead of 6.08 mm, so this is context, not a prediction.
RADEX_HUS_CH9_OXIDE_KRAD = 61.435


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--expected-files", type=int, required=True)
    parser.add_argument("--job-id", type=int, required=True)
    # Added so the same driver can report a pilot or a field-size control,
    # matching Plotting/CRDSHUSDDD.py which reads the same result files.
    parser.add_argument("--folder", default=DEFAULT_FOLDER)
    arguments = parser.parse_args()
    folder = arguments.folder

    results = aggregateRun(
        BASE_PATH / folder / "Res",
        expectedFiles=arguments.expected_files,
        tileCount=len(TILES),
    )

    print(f"\nResult files     : {results['files']}")
    print(f"Source area      : {SOURCE_AREA_CM2:.2f} cm2")
    print(f"Reference fluence: {REFERENCE_FLUENCE:.3e} electrons/cm2 "
          "(assumed)\n")

    rows = []
    records = tileRecords(results, TILES)
    for record in records:
        dose, error = record["dose"], record["error"]
        print(f"{record['label']} ({record['volume']})")
        print(f"  pooled primaries    : {record['entries']:.0f}")
        print(f"  non-zero entries    : {record['nonZeros']:.0f}")
        print(f"  hit fraction        : {record['hitFraction']:.6e}")
        print(f"  dose coefficient    : {dose:.6e} kRad cm2/electron")
        print(f"  statistical error   : {error:.6e} kRad cm2/electron")
        print(f"  relative error      : {record['relativePercent']:.4f}%")
        print(f"  dose at reference   : {dose * REFERENCE_FLUENCE:.4f} "
              f"+- {error * REFERENCE_FLUENCE:.4f} kRad")
        print()

        rows.append({
            "particle": "electron",
            "source": "HUS-PDD-SSD100",
            "tile_index": record["index"],
            "volume_name": record["volume"],
            "layer": record["label"],
            "slurm_job_id": arguments.job_id,
            "simulated_primaries": int(record["entries"]),
            "field_size_mm": "20 x 20",
            "source_area_cm2": f"{SOURCE_AREA_CM2:.6f}",
            "dose_coefficient_kRad_cm2_per_electron": f"{dose:.12e}",
            "statistical_error_kRad_cm2_per_electron": f"{error:.12e}",
            "relative_error_percent": f"{record['relativePercent']:.8f}",
            "nonzero_entries": int(record["nonZeros"]),
            "hit_fraction": f"{record['hitFraction']:.9e}",
            "reference_fluence_electrons_per_cm2": f"{REFERENCE_FLUENCE:.12e}",
            "dose_at_reference_kRad": f"{dose * REFERENCE_FLUENCE:.12e}",
            "dose_error_at_reference_kRad": f"{error * REFERENCE_FLUENCE:.12e}",
        })

    oxideDose, dieDose = records[0]["dose"], records[1]["dose"]
    oxideAtReference = oxideDose * REFERENCE_FLUENCE
    print(f"oxide / die                  : {oxideDose / dieDose:.4f}")
    print(f"oxide vs RadEx-HUS Ch9 oxide : "
          f"{100.0 * (oxideAtReference / RADEX_HUS_CH9_OXIDE_KRAD - 1.0):+.2f} % "
          "(different shielding and field size, context only)")

    output_path = BASE_PATH / folder / f"TotalDose_{folder}.csv"
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {output_path}")


if __name__ == "__main__":
    main()
