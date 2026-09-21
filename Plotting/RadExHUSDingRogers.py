"""Aggregate and compare the RadEx-HUS Ding-Rogers production TID result."""

import csv
from pathlib import Path
import sys

import numpy as np

PYTHON_ROOT = Path("/home/fetzera1/Desktop/fetzera1/Python")
sys.path.insert(0, str(PYTHON_ROOT))

from Dependencies.TotalDose import totalDose


BASE = Path("/home/fetzera1/Desktop/fetzera1/GRAS/RadEx/RadEx-HUS")
RUN = BASE / "6MeVElectron-DingRogers1995-FinalGeometry"
RESULTS = RUN / "Res"
OUTPUT = BASE / "RadEx-HUS-DingRogers1995-dose-results.csv"
COMPARISON = BASE / "RadEx-HUS-source-model-comparison.csv"
MONO_RESULTS = BASE / "RadEx-HUS-dose-results.csv"
EXPECTED_FILES = 100
EXPECTED_HISTORIES = 100_000_000
NORMALIZATION_FLUENCE = 2.0e12
SLURM_JOB = 19663755

DETECTORS = [
    (0, "VT01", "100% PE"),
    (1, "VT01", "25% Pb / 75% PE"),
    (2, "VT01", "50% Pb / 50% PE"),
    (3, "VT01", "75% Pb / 25% PE"),
    (4, "VT01", "100% Pb"),
    (5, "VT01", "6 mm Al"),
    (6, "VT01", "4 mm Al"),
    (7, "VT05", "2 mm Al"),
    (8, "VT05", "1 mm Al"),
    (9, "VT05", "exposed 0 mm"),
    (10, "VT01", "PCB back behind Ch9 / 0 mm"),
    (11, "VT01", "PCB back centre"),
]


def main() -> None:
    result_files = sorted(RESULTS.glob("*.csv"))
    if len(result_files) != EXPECTED_FILES:
        raise RuntimeError(
            f"Found {len(result_files)} result files; expected {EXPECTED_FILES}."
        )

    result = totalDose(str(RESULTS))
    histories = int(result["entries"][0])
    if histories != EXPECTED_HISTORIES:
        raise RuntimeError(
            f"Aggregated {histories} primary histories; expected {EXPECTED_HISTORIES}."
        )

    relative_error = np.divide(
        100.0 * result["error"],
        result["dose"],
        out=np.full_like(result["error"], np.nan, dtype=float),
        where=result["dose"] != 0,
    )

    fieldnames = [
        "nominal_energy_mev",
        "particle",
        "spectrum_model",
        "spectrum_scoring_plane_ssd_cm",
        "geant4_air_transport_cm",
        "tile_index",
        "channel",
        "radfet_type",
        "shielding_or_location",
        "dose_coefficient_kRad_cm2_per_electron",
        "statistical_error_kRad_cm2_per_electron",
        "relative_error_percent",
        "dose_at_2e12_electrons_cm2_kRad",
        "statistical_error_at_2e12_electrons_cm2_kRad",
        "nonzero_entries",
        "primary_histories",
        "slurm_job",
    ]

    rows = []
    for tile, (channel, radfet_type, location) in enumerate(DETECTORS):
        rows.append(
            {
                "nominal_energy_mev": 6,
                "particle": "electron",
                "spectrum_model": "Ding-Rogers-1995 Clinac-2100C all-electrons",
                "spectrum_scoring_plane_ssd_cm": 100,
                "geant4_air_transport_cm": 20,
                "tile_index": tile,
                "channel": channel,
                "radfet_type": radfet_type,
                "shielding_or_location": location,
                "dose_coefficient_kRad_cm2_per_electron": f"{result['dose'][tile]:.15e}",
                "statistical_error_kRad_cm2_per_electron": f"{result['error'][tile]:.15e}",
                "relative_error_percent": f"{relative_error[tile]:.9f}",
                "dose_at_2e12_electrons_cm2_kRad": f"{result['dose'][tile] * NORMALIZATION_FLUENCE:.9f}",
                "statistical_error_at_2e12_electrons_cm2_kRad": f"{result['error'][tile] * NORMALIZATION_FLUENCE:.9f}",
                "nonzero_entries": int(result["non-zeros"][tile]),
                "primary_histories": histories,
                "slurm_job": SLURM_JOB,
            }
        )

    with OUTPUT.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    with MONO_RESULTS.open(newline="", encoding="utf-8") as input_file:
        mono_rows = {int(row["channel"]): row for row in csv.DictReader(input_file)}

    comparison_fields = [
        "channel",
        "radfet_type",
        "shielding_or_location",
        "mono_6MeV_dose_at_2e12_kRad",
        "Ding_Rogers_dose_at_2e12_kRad",
        "Ding_Rogers_to_mono_ratio",
        "percent_change",
    ]
    comparison_rows = []
    for row in rows:
        channel = int(row["channel"])
        mono_dose = float(mono_rows[channel]["dose_at_2e12_electrons_cm2_kRad"])
        spectrum_dose = float(row["dose_at_2e12_electrons_cm2_kRad"])
        ratio = spectrum_dose / mono_dose
        comparison_rows.append(
            {
                "channel": channel,
                "radfet_type": row["radfet_type"],
                "shielding_or_location": row["shielding_or_location"],
                "mono_6MeV_dose_at_2e12_kRad": f"{mono_dose:.9f}",
                "Ding_Rogers_dose_at_2e12_kRad": f"{spectrum_dose:.9f}",
                "Ding_Rogers_to_mono_ratio": f"{ratio:.6f}",
                "percent_change": f"{100.0 * (ratio - 1.0):.3f}",
            }
        )

    with COMPARISON.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=comparison_fields)
        writer.writeheader()
        writer.writerows(comparison_rows)

    print(f"Wrote {OUTPUT}")
    print(f"Wrote {COMPARISON}")
    print(
        f"Max relative error Ch0-Ch10: {np.max(relative_error[:11]):.4f}% "
        f"(Ch{int(np.argmax(relative_error[:11]))})"
    )
    print(f"Relative error Ch11: {relative_error[11]:.4f}%")
    for row in comparison_rows:
        print(
            f"Ch{row['channel']}: {row['Ding_Rogers_dose_at_2e12_kRad']} kRad, "
            f"ratio to mono = {row['Ding_Rogers_to_mono_ratio']}"
        )


if __name__ == "__main__":
    main()
