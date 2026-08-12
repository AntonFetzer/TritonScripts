"""Compare a revised RadEx Uppsala GRAS run against the legacy geometry."""

import argparse
import csv
from pathlib import Path
import sys

import numpy as np

PYTHON_ROOT = Path("/home/fetzera1/Desktop/fetzera1/Python")
sys.path.insert(0, str(PYTHON_ROOT))

from Dependencies.TotalDose import totalDose


BASE = Path("/scratch/work/fetzera1/GRAS/RadEx/RadEx-Uppsala")
LEGACY = BASE / "RadEx-Uppsala-dose-results-legacy-noCu-noAir.csv"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("energy", type=int)
    parser.add_argument("folder")
    args = parser.parse_args()

    results = totalDose(str(BASE / args.folder / "Res"))
    with LEGACY.open(newline="", encoding="utf-8") as handle:
        legacy_rows = {
            int(row["tile_index"]): row
            for row in csv.DictReader(handle)
            if int(row["energy_MeV"]) == args.energy
        }

    output = BASE / args.folder / (
        f"Comparison_{args.energy}MeV_CuAir_vs_legacy.csv"
    )
    fields = [
        "tile_index",
        "volume_name",
        "new_dose_coefficient_kRad_cm2_per_proton",
        "new_statistical_error_kRad_cm2_per_proton",
        "new_relative_error_percent",
        "legacy_dose_coefficient_kRad_cm2_per_proton",
        "legacy_statistical_error_kRad_cm2_per_proton",
        "change_percent",
        "difference_sigma",
        "nonzero_entries",
    ]

    rows = []
    for tile in range(len(results["dose"])):
        old = legacy_rows[tile]
        old_dose = float(old["dose_coefficient_kRad_cm2_per_proton"])
        old_error = float(old["statistical_error_kRad_cm2_per_proton"])
        new_dose = float(results["dose"][tile])
        new_error = float(results["error"][tile])
        combined_error = np.hypot(new_error, old_error)
        rows.append(
            {
                "tile_index": tile,
                "volume_name": old["volume_name"],
                "new_dose_coefficient_kRad_cm2_per_proton": f"{new_dose:.12e}",
                "new_statistical_error_kRad_cm2_per_proton": f"{new_error:.12e}",
                "new_relative_error_percent": f"{100 * new_error / new_dose:.8f}",
                "legacy_dose_coefficient_kRad_cm2_per_proton": f"{old_dose:.12e}",
                "legacy_statistical_error_kRad_cm2_per_proton": f"{old_error:.12e}",
                "change_percent": f"{100 * (new_dose / old_dose - 1):.8f}",
                "difference_sigma": f"{(new_dose - old_dose) / combined_error:.6f}",
                "nonzero_entries": int(results["non-zeros"][tile]),
            }
        )

    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    print(",".join(fields))
    for row in rows:
        print(",".join(str(row[field]) for field in fields))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
