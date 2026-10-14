"""Compare the VT01 package module + region cuts against the legacy production run.

Both runs use the same measured Uppsala geometry, spectrum, beam and
normalisation. The new run replaces the legacy 1.8 x 1.9 x 0.5 mm bulk silicon
block with the SOT-23-6 package module and confines a 1 um production cut to the
twelve packages. It scores 24 volumes: the gate oxides in columns 0-11 and the
silicon dies in columns 12-23, channel by channel.

Aggregation is done by Dependencies.TotalDose.totalDose, which entry-weights the
per-file doses, pools the statistical errors and runs the Birge check.
"""

import argparse
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
SOURCE_AREA_CM2 = (66.8 * 154.0) / 100.0
N_CHANNELS = 12

DEFAULTS = {
    85: {"legacy": "85MeVProton-Gaussian-MeasuredGeometry",
         "new": "85MeVProton-Gaussian-MeasuredGeometry-NewRadFET-Cuts",
         "legacy_job": 20217122, "new_job": 20309279,
         "comparison": "RadEx-Uppsala-NewRadFET-comparison.csv",
         "plot": "NewRadFET_change_per_channel.pdf"},
    64: {"legacy": "64MeVProton-Gaussian-MeasuredGeometry",
         "new": "64MeVProton-Gaussian-MeasuredGeometry-NewRadFET-Cuts",
         "legacy_job": 20216430, "new_job": 20313480,
         "comparison": "RadEx-Uppsala-NewRadFET-comparison-64MeV.csv",
         "plot": "NewRadFET_change_per_channel_64MeV.pdf"},
}

# Channels 10 and 11 sit on the back of the readout PCB and are placed with a
# 180 degree rotation, so the beam reaches their oxide through the board, the
# leads, the paddle and the die substrate instead of through 0.35 mm of mould.
BACK_MOUNTED = (10, 11)
# Device type per channel, from README-Uppsala.md. All are modelled with the
# VT01 400 nm stack; see the note in the GDML.
VT05_CHANNELS = (7, 8, 9, 11)


def aggregate(folder, expected_volumes):
    result_path = BASE_PATH / folder / "Res"
    files = sorted(result_path.glob("*.csv"))
    results = totalDose(str(result_path))
    dose = np.asarray(results["dose"], dtype=float)
    error = np.asarray(results["error"], dtype=float)
    if len(dose) != expected_volumes:
        raise RuntimeError(
            f"{folder} scored {len(dose)} volumes, expected {expected_volumes}."
        )
    if not np.all(np.isfinite(dose)) or not np.all(np.isfinite(error)):
        raise RuntimeError(f"Non-finite values returned for {folder}.")
    if np.any(dose <= 0):
        raise RuntimeError(f"Non-positive dose in {folder}.")
    print(f"  {folder}: {len(files)} result files")
    return results, dose, error


def write_total_dose_csv(folder, dose, error, non_zeros):
    """Per-volume summary in the same layout as the sibling run folders."""
    output = BASE_PATH / folder / f"TotalDose_{folder}.csv"
    with output.open("w", encoding="utf-8") as handle:
        handle.write("Tile, Dose [kRad], Error [kRad], Non-Zero Entries\n")
        for index in range(len(dose)):
            rounded = ufloat_fromstr(f"{ufloat(dose[index], error[index]):.2u}")
            handle.write(
                f"{index},{rounded.n},{rounded.s},{int(non_zeros[index])}\n"
            )
    return output


def write_diagnostic_plots(folder, dose, error, non_zeros):
    indices = np.arange(len(dose))
    folder_path = BASE_PATH / folder

    plt.figure()
    plt.errorbar(indices, dose, yerr=error, fmt=" ", capsize=5, elinewidth=1,
                 capthick=1, label="Dose")
    plt.title("Dose per tile " + folder)
    plt.xlabel("Tile number")
    plt.ylabel("Dose [kRad]")
    plt.yscale("log")
    plt.grid(which="both")
    plt.legend()
    plt.savefig(folder_path / f"Dose_{folder}.pdf", format="pdf",
                bbox_inches="tight")

    # Linear axis starting at zero: the point is to spot high-error tiles.
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


def compare(new_dose, new_error, legacy_dose, legacy_error):
    ratio = new_dose / legacy_dose
    change = 100.0 * (ratio - 1.0)
    combined = np.hypot(new_error, legacy_error)
    sigma = (new_dose - legacy_dose) / combined
    return ratio, change, combined, sigma


def weighted_mean(change, combined, legacy_dose, channels):
    idx = list(channels)
    rel = 100.0 * combined[idx] / legacy_dose[idx]
    weights = 1.0 / rel ** 2
    return float(np.sum(change[idx] * weights) / np.sum(weights))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--energy", type=int, choices=sorted(DEFAULTS),
                        default=85, help="nominal beam energy in MeV")
    args = parser.parse_args()
    cfg = DEFAULTS[args.energy]
    energy_mev = args.energy
    legacy_folder = cfg["legacy"]
    new_folder = cfg["new"]
    comparison_path = BASE_PATH / cfg["comparison"]
    legacy_job_id = cfg["legacy_job"]
    new_job_id = cfg["new_job"]
    plot_path = BASE_PATH / cfg["plot"]

    print(f"Aggregating {energy_mev} MeV")
    legacy_results, legacy_dose, legacy_error = aggregate(legacy_folder, N_CHANNELS)
    new_results, new_dose, new_error = aggregate(new_folder, 2 * N_CHANNELS)

    write_total_dose_csv(new_folder, new_dose, new_error,
                         new_results["non-zeros"])
    write_diagnostic_plots(new_folder, new_dose, new_error,
                           new_results["non-zeros"])

    gox_dose, gox_error = new_dose[:N_CHANNELS], new_error[:N_CHANNELS]
    die_dose, die_error = new_dose[N_CHANNELS:], new_error[N_CHANNELS:]
    gox_nz = np.asarray(new_results["non-zeros"], dtype=float)[:N_CHANNELS]
    die_nz = np.asarray(new_results["non-zeros"], dtype=float)[N_CHANNELS:]

    gox_ratio, gox_change, gox_combined, gox_sigma = compare(
        gox_dose, gox_error, legacy_dose, legacy_error)
    die_ratio, die_change, die_combined, die_sigma = compare(
        die_dose, die_error, legacy_dose, legacy_error)

    fieldnames = [
        "energy_MeV", "channel", "device_type", "mounting",
        "legacy_job_id", "new_job_id", "source_area_cm2",
        "legacy_dose_coefficient", "legacy_statistical_error",
        "oxide_dose_coefficient", "oxide_statistical_error",
        "oxide_nonzero_entries", "oxide_to_legacy_ratio",
        "oxide_change_percent", "oxide_significance_sigma",
        "die_dose_coefficient", "die_statistical_error",
        "die_nonzero_entries", "die_to_legacy_ratio",
        "die_change_percent", "die_significance_sigma",
    ]
    rows = []
    for ch in range(N_CHANNELS):
        rows.append({
            "energy_MeV": energy_mev,
            "channel": ch,
            "device_type": "VT05" if ch in VT05_CHANNELS else "VT01",
            "mounting": "readout_pcb_back" if ch in BACK_MOUNTED else "sensor_pcb_front",
            "legacy_job_id": legacy_job_id,
            "new_job_id": new_job_id,
            "source_area_cm2": f"{SOURCE_AREA_CM2:.6f}",
            "legacy_dose_coefficient": f"{legacy_dose[ch]:.12e}",
            "legacy_statistical_error": f"{legacy_error[ch]:.12e}",
            "oxide_dose_coefficient": f"{gox_dose[ch]:.12e}",
            "oxide_statistical_error": f"{gox_error[ch]:.12e}",
            "oxide_nonzero_entries": int(gox_nz[ch]),
            "oxide_to_legacy_ratio": f"{gox_ratio[ch]:.9f}",
            "oxide_change_percent": f"{gox_change[ch]:+.6f}",
            "oxide_significance_sigma": f"{gox_sigma[ch]:+.4f}",
            "die_dose_coefficient": f"{die_dose[ch]:.12e}",
            "die_statistical_error": f"{die_error[ch]:.12e}",
            "die_nonzero_entries": int(die_nz[ch]),
            "die_to_legacy_ratio": f"{die_ratio[ch]:.9f}",
            "die_change_percent": f"{die_change[ch]:+.6f}",
            "die_significance_sigma": f"{die_sigma[ch]:+.4f}",
        })

    with comparison_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    channels = np.arange(N_CHANNELS)
    plt.figure(figsize=(8, 5))
    plt.errorbar(channels - 0.08, gox_change,
                 yerr=100.0 * gox_combined / legacy_dose,
                 fmt="o", capsize=4, elinewidth=1, capthick=1,
                 label="gate oxide (SiO2, 400 nm) vs legacy block")
    plt.errorbar(channels + 0.08, die_change,
                 yerr=100.0 * die_combined / legacy_dose,
                 fmt="s", capsize=4, elinewidth=1, capthick=1,
                 label="silicon die (280 um) vs legacy block")
    plt.axhline(y=0, color="r", linestyle="--", label="no change")
    plt.axvline(x=9.5, color="grey", linestyle=":",
                label="left: sensor PCB front / right: readout PCB back")
    plt.title(f"RadFET package model vs legacy bulk silicon, {energy_mev} MeV")
    plt.xlabel("Channel")
    plt.ylabel("Change in dose [%]")
    plt.xticks(channels)
    plt.grid(which="both")
    plt.legend(fontsize=8)
    plt.savefig(plot_path, format="pdf",
                bbox_inches="tight")
    plt.close("all")

    print("\nDose coefficients [kRad cm2 / proton]")
    head = (f"{'Ch':>3s} {'dev':>5s} {'mount':>6s} {'legacy':>12s} "
            f"{'oxide':>12s} {'ox %':>8s} {'ox err':>7s} "
            f"{'die':>12s} {'die %':>8s}")
    print(head)
    for ch in range(N_CHANNELS):
        mount = "back" if ch in BACK_MOUNTED else "front"
        dev = "VT05" if ch in VT05_CHANNELS else "VT01"
        rel = 100.0 * gox_error[ch] / gox_dose[ch]
        print(f"{ch:>3d} {dev:>5s} {mount:>6s} {legacy_dose[ch]:>12.4e} "
              f"{gox_dose[ch]:>12.4e} {gox_change[ch]:>+7.2f}% {rel:>6.2f}% "
              f"{die_dose[ch]:>12.4e} {die_change[ch]:>+7.2f}%")

    front = [c for c in range(N_CHANNELS) if c not in BACK_MOUNTED]
    print("\nError-weighted mean change vs legacy")
    print(f"  oxide, front-mounted channels 0-9 : "
          f"{weighted_mean(gox_change, gox_combined, legacy_dose, front):+.2f} %")
    print(f"  oxide, back-mounted channels 10-11: "
          f"{weighted_mean(gox_change, gox_combined, legacy_dose, BACK_MOUNTED):+.2f} %")
    print(f"  die,   front-mounted channels 0-9 : "
          f"{weighted_mean(die_change, die_combined, legacy_dose, front):+.2f} %")
    print(f"  die,   back-mounted channels 10-11: "
          f"{weighted_mean(die_change, die_combined, legacy_dose, BACK_MOUNTED):+.2f} %")

    print(f"\nOxide relative errors: min {100.0 * np.min(gox_error / gox_dose):.2f} %, "
          f"max {100.0 * np.max(gox_error / gox_dose):.2f} %")
    print(f"Die relative errors:   min {100.0 * np.min(die_error / die_dose):.2f} %, "
          f"max {100.0 * np.max(die_error / die_dose):.2f} %")
    print(f"\nWrote {comparison_path}")


if __name__ == "__main__":
    main()
