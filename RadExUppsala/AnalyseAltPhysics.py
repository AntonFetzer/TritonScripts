"""Compare the QGSP_BIC_HP + em_standard_opt4 run against the QBBC + opt3 baseline.

Both runs use the same measured Uppsala geometry, the same Skandion 85 MeV
Gaussian spectrum, the same beam and normalisation, and 1e9 primaries, so the
only difference is the physics list. Aggregation is done by
Dependencies.TotalDose.totalDose, which entry-weights the per-file doses, pools
the statistical errors and runs the Birge consistency check.
"""

import csv
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat, ufloat_fromstr

PYTHON_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_ROOT))

from Dependencies.TotalDose import totalDose  # noqa: E402

BASE_PATH = Path("/scratch/work/fetzera1/GRAS/RadEx/RadEx-Uppsala")
BASELINE_FOLDER = "85MeVProton-Gaussian-MeasuredGeometry"
ALT_FOLDER = "85MeVProton-Gaussian-MeasuredGeometry-AltPhysics"
COMPARISON_PATH = BASE_PATH / "RadEx-Uppsala-AltPhysics-comparison.csv"

ENERGY_MEV = 85
SOURCE_AREA_CM2 = (66.8 * 154.0) / 100.0
BASELINE_JOB_ID = 20217122
ALT_JOB_ID = 20305134
EXPECTED_FILES = {BASELINE_FOLDER: 100, ALT_FOLDER: 500}


def aggregate(folder):
    result_path = BASE_PATH / folder / "Res"
    files = sorted(result_path.glob("*.csv"))
    expected = EXPECTED_FILES[folder]
    if len(files) != expected:
        raise RuntimeError(
            f"{folder} has {len(files)} CSV files, expected {expected}."
        )
    results = totalDose(str(result_path))
    dose = np.asarray(results["dose"], dtype=float)
    error = np.asarray(results["error"], dtype=float)
    if not np.all(np.isfinite(dose)) or not np.all(np.isfinite(error)):
        raise RuntimeError(f"Non-finite values returned for {folder}.")
    if np.any(dose <= 0):
        raise RuntimeError(f"Non-positive dose in {folder}.")
    return results, dose, error


def write_total_dose_csv(folder, dose, error, non_zeros):
    """Per-channel summary in the same layout as the sibling run folders.

    Round through the uncertainties package, then read the rounded value back
    with ufloat_fromstr, exactly as Plotting/RadExTID.py does. That round trip
    is what copes with the factored '(1.3975+/-0.0026)e-10' shape the package
    emits at these magnitudes; splitting the string on '+/-' does not.
    """
    output = BASE_PATH / folder / f"TotalDose_{folder}.csv"
    with output.open("w", encoding="utf-8") as handle:
        handle.write("Tile, Dose [kRad], Error [kRad], Non-Zero Entries\n")
        for tile in range(len(dose)):
            rounded = ufloat_fromstr(f"{ufloat(dose[tile], error[tile]):.2u}")
            handle.write(
                f"{tile},{rounded.n:.4e},{rounded.s:.1e},{int(non_zeros[tile])}\n"
            )
    return output


def write_diagnostic_plots(folder, dose, error, non_zeros):
    tiles = np.arange(len(dose))
    folder_path = BASE_PATH / folder

    plt.figure()
    plt.errorbar(tiles, dose, yerr=error, fmt=" ", capsize=5, elinewidth=1,
                 capthick=1, label="Dose")
    plt.title("Dose per tile " + folder)
    plt.xlabel("Tile number")
    plt.ylabel("Dose [kRad]")
    plt.yscale("log")
    plt.grid(which="both")
    plt.legend()
    plt.savefig(folder_path / f"Dose_{folder}.pdf", format="pdf",
                bbox_inches="tight")

    # Linear axis starting at zero: the point is to spot high-error tiles, not
    # to resolve the small ones.
    plt.figure()
    plt.plot(100.0 * error / dose, ".", label="Relative Error")
    plt.axhline(y=1, color="r", linestyle="--", label="1% error")
    plt.title("Relative Error in %")
    plt.xlabel("Tile number")
    plt.ylabel("Relative Error [%]")
    plt.ylim(bottom=0)
    plt.grid(which="both")
    plt.legend()
    plt.savefig(folder_path / f"Error_{folder}.pdf", format="pdf",
                bbox_inches="tight")

    plt.figure()
    plt.plot(non_zeros, ".", label="Non-zero entries")
    plt.title("Number of non-zero entries")
    plt.xlabel("Tile number")
    plt.ylabel("Number of non-zero entries")
    plt.yscale("log")
    plt.grid(axis="x", which="both")
    plt.grid(axis="y", which="major")
    plt.legend()
    plt.savefig(folder_path / f"NonZeros_{folder}.pdf", format="pdf",
                bbox_inches="tight")
    plt.close("all")


def main() -> None:
    base_results, base_dose, base_error = aggregate(BASELINE_FOLDER)
    alt_results, alt_dose, alt_error = aggregate(ALT_FOLDER)

    if len(base_dose) != len(alt_dose):
        raise RuntimeError("The two runs scored different numbers of volumes.")

    write_total_dose_csv(ALT_FOLDER, alt_dose, alt_error,
                         alt_results["non-zeros"])
    write_diagnostic_plots(ALT_FOLDER, alt_dose, alt_error,
                           alt_results["non-zeros"])

    ratio = alt_dose / base_dose
    change_percent = 100.0 * (ratio - 1.0)
    combined_error = np.hypot(base_error, alt_error)
    significance = (alt_dose - base_dose) / combined_error

    fieldnames = [
        "energy_MeV", "tile_index", "volume_name",
        "baseline_physics", "altphysics_physics",
        "baseline_slurm_job_id", "altphysics_slurm_job_id",
        "simulated_primaries", "source_area_cm2",
        "baseline_dose_coefficient", "baseline_statistical_error",
        "altphysics_dose_coefficient", "altphysics_statistical_error",
        "altphysics_to_baseline_ratio", "change_percent",
        "combined_significance_sigma",
    ]
    rows = []
    for tile in range(len(base_dose)):
        rows.append({
            "energy_MeV": ENERGY_MEV,
            "tile_index": tile,
            "volume_name": f"Si_vol_{tile}_PV",
            "baseline_physics": "QBBC+em_standard_opt3",
            "altphysics_physics": "QGSP_BIC_HP+em_standard_opt4",
            "baseline_slurm_job_id": BASELINE_JOB_ID,
            "altphysics_slurm_job_id": ALT_JOB_ID,
            "simulated_primaries": int(base_results["entries"][0]),
            "source_area_cm2": f"{SOURCE_AREA_CM2:.6f}",
            "baseline_dose_coefficient": f"{base_dose[tile]:.12e}",
            "baseline_statistical_error": f"{base_error[tile]:.12e}",
            "altphysics_dose_coefficient": f"{alt_dose[tile]:.12e}",
            "altphysics_statistical_error": f"{alt_error[tile]:.12e}",
            "altphysics_to_baseline_ratio": f"{ratio[tile]:.9f}",
            "change_percent": f"{change_percent[tile]:+.8f}",
            "combined_significance_sigma": f"{significance[tile]:+.5f}",
        })

    with COMPARISON_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    tiles = np.arange(len(base_dose))
    plt.figure()
    plt.errorbar(tiles, change_percent,
                 yerr=100.0 * combined_error / base_dose,
                 fmt="o", capsize=4, elinewidth=1, capthick=1,
                 label="QGSP_BIC_HP + opt4 vs QBBC + opt3")
    plt.axhline(y=0, color="r", linestyle="--", label="no change")
    plt.title(f"Physics-list effect at {ENERGY_MEV} MeV, measured geometry")
    plt.xlabel("Tile number")
    plt.ylabel("Change in dose [%]")
    plt.xticks(tiles)
    plt.grid(which="both")
    plt.legend()
    plt.savefig(BASE_PATH / "AltPhysics_change_per_channel.pdf", format="pdf",
                bbox_inches="tight")
    plt.close("all")

    header = f"{'Ch':>3s} {'QBBC+opt3':>14s} {'BIC_HP+opt4':>14s} {'change':>9s} {'sigma':>8s}"
    print("\nPer-channel comparison (dose coefficients in kRad cm2 / proton)")
    print(header)
    for tile in range(len(base_dose)):
        print(f"{tile:>3d} {base_dose[tile]:>14.5e} {alt_dose[tile]:>14.5e} "
              f"{change_percent[tile]:>+8.3f}% {significance[tile]:>+8.2f}")

    weights = 1.0 / combined_error ** 2
    mean_change = float(np.sum(change_percent * weights) / np.sum(weights))
    worst = int(np.argmax(np.abs(change_percent)))
    print(f"\nError-weighted mean change over all channels: {mean_change:+.3f} %")
    print(f"Largest change: channel {worst} at {change_percent[worst]:+.3f} %")
    print(f"Channels differing by more than 3 sigma: "
          f"{int(np.sum(np.abs(significance) > 3))} of {len(base_dose)}")
    print(f"\nWrote {COMPARISON_PATH}")


if __name__ == "__main__":
    main()
