"""Aggregate production TID results for the RadEx Uppsala simulations."""

import csv
from pathlib import Path
import sys

import numpy as np

PYTHON_ROOT = Path("/home/fetzera1/Desktop/fetzera1/Python")
sys.path.insert(0, str(PYTHON_ROOT))

from Dependencies.TotalDose import totalDose


BASE_PATH = Path("/scratch/work/fetzera1/GRAS/RadEx/RadEx-Uppsala")
OUTPUT_PATH = BASE_PATH / "RadEx-Uppsala-dose-results.csv"
SOURCE_AREA_CM2 = (66.8 * 154.0) / 100.0
NORMALIZATION_FLUENCE_PROTONS_PER_CM2 = 1.0
VOLUME_NAMES = [f"Si_vol_{tile}_PV" for tile in range(12)]
RUNS = [
    {
        "energy_MeV": 64,
        "folder": "64MeVProton-FinalGeometry",
        "slurm_job_id": 19529191,
        "expected_files": 100,
    },
    {
        "energy_MeV": 85,
        "folder": "85MeVProton-FinalGeometry",
        "slurm_job_id": 19529141,
        "expected_files": 50,
    },
]


def main() -> None:
    fieldnames = [
        "energy_MeV",
        "particle",
        "tile_index",
        "volume_name",
        "slurm_job_id",
        "simulated_primaries",
        "source_area_cm2",
        "normalization_fluence_protons_per_cm2",
        "dose_coefficient_kRad_cm2_per_proton",
        "statistical_error_kRad_cm2_per_proton",
        "relative_error_percent",
        "nonzero_entries",
    ]

    rows = []
    for run in RUNS:
        result_path = BASE_PATH / run["folder"] / "Res"
        result_files = list(result_path.glob("*.csv"))
        if len(result_files) != run["expected_files"]:
            raise RuntimeError(
                f"{run['folder']} has {len(result_files)} CSV files; "
                f"expected {run['expected_files']}. Production is incomplete."
            )

        results = totalDose(str(result_path))
        simulated_primaries = int(results["entries"][0])
        relative_error_percent = np.divide(
            100.0 * results["error"],
            results["dose"],
            out=np.zeros_like(results["error"], dtype=float),
            where=results["dose"] != 0,
        )

        for tile, volume_name in enumerate(VOLUME_NAMES):
            rows.append(
                {
                    "energy_MeV": run["energy_MeV"],
                    "particle": "proton",
                    "tile_index": tile,
                    "volume_name": volume_name,
                    "slurm_job_id": run["slurm_job_id"],
                    "simulated_primaries": simulated_primaries,
                    "source_area_cm2": f"{SOURCE_AREA_CM2:.6f}",
                    "normalization_fluence_protons_per_cm2": (
                        f"{NORMALIZATION_FLUENCE_PROTONS_PER_CM2:.1f}"
                    ),
                    "dose_coefficient_kRad_cm2_per_proton": (
                        f"{results['dose'][tile]:.12e}"
                    ),
                    "statistical_error_kRad_cm2_per_proton": (
                        f"{results['error'][tile]:.12e}"
                    ),
                    "relative_error_percent": (
                        f"{relative_error_percent[tile]:.8f}"
                    ),
                    "nonzero_entries": int(results["non-zeros"][tile]),
                }
            )

        front_max = float(np.max(relative_error_percent[:9]))
        all_max = float(np.max(relative_error_percent))
        print(
            f"{run['energy_MeV']} MeV: {simulated_primaries} primaries; "
            f"max relative error tiles 0-8 = {front_max:.4f}%; "
            f"max relative error all tiles = {all_max:.4f}%"
        )

    with OUTPUT_PATH.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
