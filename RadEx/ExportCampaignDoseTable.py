"""Consolidated RadFET dose table for the proton campaign analysis.

Emits the long-form schema that sw/V2_Proton/Analysis/ProtonAnalysis.py's
load_geant4_doses expects, one row per (energy, channel), from GRAS runs that
score the VT01 package module: gate oxides in the first half of the tallies,
silicon dies in the second, channel by channel.

Each --run names the label written to the energy_MeV column (the measured phase
the run is compared with, e.g. 64) and a run directory holding A.mac, Res/ and
log/. Everything else comes from the run itself: the pooled doses and the
primaries from Res/ through Dependencies.TotalDose.totalDose, the Slurm job from
the log names, and the spectrum centre and width, the source area and the
normalisation fluence from its macros (RadEx.RunSettings). Pointing the
tool at other runs is therefore enough to tabulate them.

The scored quantity is the 400 nm gate oxide, which is what the RadFET
calibration refers to. The silicon die coefficient is carried alongside as an
extra column for reference; load_geant4_doses selects columns by name and
ignores it.

Aggregation is Dependencies.TotalDose.totalDose, so the coefficients carry full
precision rather than the two-significant-figure rounding in the per-run
TotalDose CSVs.

Run from the repository root:

    python3 -m RadEx.ExportCampaignDoseTable --run 64 <64 MeV run dir> \\
        --run 85 <85 MeV run dir> --output <table.csv>
"""

import argparse
import csv
from pathlib import Path
import sys

import numpy as np

PYTHON_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_ROOT))
from Dependencies.TotalDose import totalDose  # noqa: E402
from RadEx import RunSettings  # noqa: E402

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


def run_rows(label, run):
    """Table rows of one run, oxide coefficients first, the die as reference."""
    result_path = run / "Res"
    files = sorted(result_path.glob("*.csv"))
    results = totalDose(str(result_path))
    dose = np.asarray(results["dose"], dtype=float)
    error = np.asarray(results["error"], dtype=float)
    non_zeros = np.asarray(results["non-zeros"], dtype=float)
    if len(dose) % 2:
        raise RuntimeError(
            f"{run} scored {len(dose)} volumes; expected an oxide and a die per channel")
    if not np.all(np.isfinite(dose)) or not np.all(np.isfinite(error)):
        raise RuntimeError(f"Non-finite values in {run}")
    if np.any(dose <= 0):
        raise RuntimeError(f"Non-positive dose in {run}")
    channels = len(dose) // 2

    settings = RunSettings.runSettings(run)
    spectrum = RunSettings.spectrum(run, settings)
    jobs = RunSettings.slurmJobs(run)
    if not jobs:
        raise RuntimeError(f"{run}: no Slurm logs in log/ to name the job")
    primaries = int(results["entries"][0])
    area = RunSettings.sourceAreaCm2(run, settings)
    fluence = RunSettings.normalisationFluence(run, settings)

    print(f"{label}: {run.name}, {len(files)} files, {primaries:.4g} primaries, "
          f"job {';'.join(map(str, jobs))}, {spectrum['centre_MeV']:.2f} MeV "
          f"sigma {spectrum['sigma_MeV']:.5f}")
    rows = []
    for channel in range(channels):
        oxide, oxide_error = dose[channel], error[channel]
        die, die_error = dose[channels + channel], error[channels + channel]
        rows.append({
            "energy_MeV": label,
            "particle": "proton",
            "tile_index": channel,
            "volume_name": f"VT01_gox_{channel}_PV",
            "slurm_job_id": ";".join(map(str, jobs)),
            "simulated_primaries": primaries,
            "source_area_cm2": f"{area:.6f}",
            "normalization_fluence_protons_per_cm2": repr(fluence),
            "dose_coefficient_kRad_cm2_per_proton": f"{oxide:.12e}",
            "statistical_error_kRad_cm2_per_proton": f"{oxide_error:.12e}",
            "relative_error_percent": f"{100.0 * oxide_error / oxide:.8f}",
            "nonzero_entries": int(non_zeros[channel]),
            "spectrum_centre_MeV": f"{spectrum['centre_MeV']:.2f}",
            "spectrum_sigma_MeV": f"{spectrum['sigma_MeV']:.5f}",
            "detector_model": "VT01 package module, 400 nm gate oxide",
            "die_dose_coefficient_kRad_cm2_per_proton": f"{die:.12e}",
            "die_relative_error_percent": f"{100.0 * die_error / die:.8f}",
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run", nargs=2, action="append", required=True,
                        metavar=("LABEL", "RUN_DIR"),
                        help="energy_MeV label and run directory; repeat per run, in table order")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    rows = []
    for label, run in arguments.run:
        rows.extend(run_rows(label, Path(run)))

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    with arguments.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {arguments.output} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
