"""Consolidated dose table for the proton campaign analysis.

Emits the long-form schema that sw/V2_Proton/Analysis/ProtonAnalysis.py's
load_geant4_doses expects, one row per (energy, channel), built from the
refitted-spectrum NewRadFET production runs rather than the legacy
bulk-silicon, nominal-energy ones.

The scored quantity is the 400 nm gate oxide, which is what the RadFET
calibration refers to. The silicon die coefficient is carried alongside as an
extra column for reference; load_geant4_doses selects columns by name and
ignores it.

Aggregation is Dependencies.TotalDose.totalDose, so the coefficients carry full
precision rather than the two-significant-figure rounding in the per-run
TotalDose CSVs.
"""

import argparse
import csv
from pathlib import Path
import sys

import numpy as np

BASE_PATH = Path("/scratch/work/fetzera1/GRAS/RadEx/RadEx-Uppsala")
PYTHON_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_ROOT))
from Dependencies.TotalDose import totalDose  # noqa: E402

N_CHANNELS = 12
SOURCE_AREA_CM2 = (66.8 * 154.0) / 100.0

RUNS = {
    64: {"folder": "64MeVProton-NewGaussian-MeasuredGeometry-NewRadFET-Cuts",
         "job": 20328954, "primaries": 2_000_000_000,
         "centre_MeV": 62.07, "sigma_MeV": 0.49269},
    85: {"folder": "85MeVProton-NewGaussian-MeasuredGeometry-NewRadFET-Cuts",
         "job": 20328955, "primaries": 1_000_000_000,
         "centre_MeV": 83.47, "sigma_MeV": 0.60875},
}

FIELDNAMES = [
    "energy_MeV", "particle", "tile_index", "volume_name", "slurm_job_id",
    "simulated_primaries", "source_area_cm2",
    "normalization_fluence_protons_per_cm2",
    "dose_coefficient_kRad_cm2_per_proton",
    "statistical_error_kRad_cm2_per_proton", "relative_error_percent",
    "nonzero_entries",
    # Extra provenance and reference columns, ignored by load_geant4_doses.
    "spectrum_centre_MeV", "spectrum_sigma_MeV", "detector_model",
    "die_dose_coefficient_kRad_cm2_per_proton", "die_relative_error_percent",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    rows = []
    for energy, run in sorted(RUNS.items()):
        result_path = BASE_PATH / run["folder"] / "Res"
        files = sorted(result_path.glob("*.csv"))
        results = totalDose(str(result_path))
        dose = np.asarray(results["dose"], dtype=float)
        error = np.asarray(results["error"], dtype=float)
        non_zeros = np.asarray(results["non-zeros"], dtype=float)
        if len(dose) != 2 * N_CHANNELS:
            raise RuntimeError(
                f"{run['folder']} scored {len(dose)} volumes, expected 24")
        if not np.all(np.isfinite(dose)) or not np.all(np.isfinite(error)):
            raise RuntimeError(f"Non-finite values in {run['folder']}")
        if np.any(dose <= 0):
            raise RuntimeError(f"Non-positive dose in {run['folder']}")

        print(f"{energy} MeV: {len(files)} files, {run['primaries']:.3g} primaries")
        for channel in range(N_CHANNELS):
            oxide, oxide_error = dose[channel], error[channel]
            die, die_error = dose[N_CHANNELS + channel], error[N_CHANNELS + channel]
            rows.append({
                "energy_MeV": energy,
                "particle": "proton",
                "tile_index": channel,
                "volume_name": f"VT01_gox_{channel}_PV",
                "slurm_job_id": run["job"],
                "simulated_primaries": run["primaries"],
                "source_area_cm2": f"{SOURCE_AREA_CM2:.6f}",
                "normalization_fluence_protons_per_cm2": "1.0",
                "dose_coefficient_kRad_cm2_per_proton": f"{oxide:.12e}",
                "statistical_error_kRad_cm2_per_proton": f"{oxide_error:.12e}",
                "relative_error_percent": f"{100.0 * oxide_error / oxide:.8f}",
                "nonzero_entries": int(non_zeros[channel]),
                "spectrum_centre_MeV": f"{run['centre_MeV']:.2f}",
                "spectrum_sigma_MeV": f"{run['sigma_MeV']:.5f}",
                "detector_model": "VT01 package module, 400 nm gate oxide",
                "die_dose_coefficient_kRad_cm2_per_proton": f"{die:.12e}",
                "die_relative_error_percent": f"{100.0 * die_error / die:.8f}",
            })

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    with arguments.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {arguments.output} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
