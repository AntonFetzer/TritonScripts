"""Aggregate the CRDS Uppsala 64 MeV proton production run.

CRDS in the Skandion 64 MeV beam: one VT01 RadFET package behind two 1.6 mm
shielding PCBs, standing on the beam-facing side of a third. Two tallies, the
400 nm gate oxide and the silicon die behind it, in the order
``CRDS/CRDS1RadFETDetector.mac`` fixes.

Sibling of RadExUppsalaProduction.py, which does the same job for the twelve
RadEx channels in the same beam.

    python3 Plotting/CRDSUppsalaProduction.py --expected-files 799 --job-id 20330228

``--expected-files`` guards against aggregating a partial run; pass the number
of files actually produced. See Dependencies.AggregateRun for why that is not
always the array size.

Two scope limits worth knowing before quoting anything this prints:

* The delivered fluence below covers the 5.2 x 5.4 cm aggregate only, which is
  what this simulation models. The campaign also delivered a 7.0 x 7.4 cm
  exposure that is not simulated at all; see ``GRAS/CRDS/README.md``.
* The legacy bulk-silicon run this driver used to report a "vs legacy" column
  against has been removed from the GRAS tree. That comparison is recorded in
  ``GRAS/CRDS/CRDS-Uppsala/README-Uppsala.md`` and is not recomputed here.
"""

import argparse
import csv
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.AggregateRun import aggregateRun, tileRecords  # noqa: E402

BASE_PATH = Path("/scratch/work/fetzera1/GRAS/CRDS/CRDS-Uppsala")
FOLDER = "64MeVProton"

# Tally order is fixed by CRDS/CRDS1RadFETDetector.mac. Change both together.
TILES = [
    (0, "VT01_gox_0_PV", "gate oxide"),
    (1, "VT01_die_0_PV", "silicon die"),
]

# Source plane is 52 x 54 mm; GRAS normalises to one incident proton/cm2, so
# the tallied result is a coefficient in kRad cm2/proton.
SOURCE_AREA_CM2 = (52.0 * 54.0) / 100.0

# Delivered beam for the 5.2 x 5.4 cm CRDS exposures, from the Uppsala
# accelerator logs: aggregate CRDS_5.2x5.4cm, experiments exp_10 to exp_12.
DELIVERED_PROTONS = 5.710568363210e12
SCANNED_AREA_CM2 = 5.2 * 5.4
DELIVERED_FLUENCE = DELIVERED_PROTONS / SCANNED_AREA_CM2


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--expected-files", type=int, required=True)
    parser.add_argument("--job-id", type=int, required=True)
    arguments = parser.parse_args()

    results = aggregateRun(
        BASE_PATH / FOLDER / "Res",
        expectedFiles=arguments.expected_files,
        tileCount=len(TILES),
    )

    print(f"\nResult files     : {results['files']}")
    print(f"Source area      : {SOURCE_AREA_CM2:.2f} cm2")
    print(f"Delivered fluence: {DELIVERED_FLUENCE:.6e} protons/cm2 "
          "(5.2 x 5.4 cm exposures only)\n")

    rows = []
    for record in tileRecords(results, TILES):
        dose, error = record["dose"], record["error"]
        print(f"{record['label']} ({record['volume']})")
        print(f"  pooled primaries    : {record['entries']:.0f}")
        print(f"  non-zero entries    : {record['nonZeros']:.0f}")
        print(f"  hit fraction        : {record['hitFraction']:.6e}")
        print(f"  dose coefficient    : {dose:.6e} kRad cm2/proton")
        print(f"  statistical error   : {error:.6e} kRad cm2/proton")
        print(f"  relative error      : {record['relativePercent']:.4f}%")
        print(f"  delivered dose      : {dose * DELIVERED_FLUENCE:.4f} "
              f"+- {error * DELIVERED_FLUENCE:.4f} kRad")
        print()

        rows.append({
            "energy_MeV": 64,
            "particle": "proton",
            "tile_index": record["index"],
            "volume_name": record["volume"],
            "layer": record["label"],
            "slurm_job_id": arguments.job_id,
            "simulated_primaries": int(record["entries"]),
            "source_area_cm2": f"{SOURCE_AREA_CM2:.6f}",
            "dose_coefficient_kRad_cm2_per_proton": f"{dose:.12e}",
            "statistical_error_kRad_cm2_per_proton": f"{error:.12e}",
            "relative_error_percent": f"{record['relativePercent']:.8f}",
            "nonzero_entries": int(record["nonZeros"]),
            "hit_fraction": f"{record['hitFraction']:.9e}",
            "delivered_fluence_protons_per_cm2": f"{DELIVERED_FLUENCE:.12e}",
            "delivered_dose_kRad": f"{dose * DELIVERED_FLUENCE:.12e}",
            "delivered_dose_error_kRad": f"{error * DELIVERED_FLUENCE:.12e}",
        })

    oxide, die = rows[0], rows[1]
    ratio = float(oxide["dose_coefficient_kRad_cm2_per_proton"]) / float(
        die["dose_coefficient_kRad_cm2_per_proton"])
    print(f"oxide / die : {ratio:.4f}")

    output_path = BASE_PATH / FOLDER / f"TotalDose_{FOLDER}.csv"
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {output_path}")


if __name__ == "__main__":
    main()
