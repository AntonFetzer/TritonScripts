"""Inspect the RadEx-HUS monoenergetic PDD pilot at the measured depths."""

from __future__ import annotations

import csv
from pathlib import Path
import re
import sys

import matplotlib.pyplot as plt
import numpy as np

PYTHON_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_ROOT))

from Read.ReadDose import readDoseModules


DEFAULT_BASE = Path("/scratch/work/fetzera1/GRAS/RadEx/RadEx-HUS")
DEPTHS_MM = np.array([
    0.0, 1.5, 3.0, 4.5, 6.0, 7.5, 9.0, 10.5, 12.0, 13.5, 15.0,
    16.5, 18.0, 19.5, 21.0, 22.5, 24.0, 25.5, 27.0, 28.5, 30.0,
    31.5, 33.0, 34.5, 36.0, 37.5, 39.0, 40.0, 42.5, 45.0, 47.5,
    50.0, 55.0, 60.0,
])


def cell_bounds(depths):
    boundaries = np.empty(len(depths) + 1)
    boundaries[0] = 0.0
    boundaries[1:-1] = (depths[:-1] + depths[1:]) / 2
    boundaries[-1] = depths[-1] + (depths[-1] - depths[-2]) / 2
    return boundaries


def read_mcc(path):
    depths = []
    values = []
    reading = False
    with path.open("r", encoding="latin-1") as stream:
        for line in stream:
            stripped = line.strip()
            if stripped == "BEGIN_DATA":
                reading = True
                continue
            if stripped == "END_DATA":
                break
            if reading and stripped:
                fields = stripped.split()
                depths.append(float(fields[0]))
                values.append(float(fields[1]))
    return np.asarray(depths), np.asarray(values)


def energy_from_name(path):
    match = re.search(r"PDD_(\d+)p(\d+)MeV", path.name)
    if not match:
        raise ValueError(f"Cannot parse energy from {path.name}")
    return float(f"{match.group(1)}.{match.group(2)}")


def falling_crossing(depths, pdd, level):
    peak = int(np.argmax(pdd))
    for right in range(peak + 1, len(pdd)):
        if pdd[right] <= level < pdd[right - 1]:
            left = right - 1
            fraction = (level - pdd[left]) / (pdd[right] - pdd[left])
            return depths[left] + fraction * (depths[right] - depths[left])
    return np.nan


def main():
    base = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_BASE
    run = base / "PDD-Spectrum-Pilot"
    result_files = sorted((run / "Res").glob("PDD_*MeV_*.csv"))
    if len(result_files) != 3:
        raise RuntimeError(f"Found {len(result_files)} result files; expected 3")

    measured_depths, measured = read_mcc(next(base.glob("*eHDTSE_PDD.mcc")))
    if not np.array_equal(measured_depths, DEPTHS_MM):
        raise RuntimeError("Measured MCC depths do not match the pilot layer labels")
    measured_pdd = 100 * measured / np.max(measured)

    boundaries = cell_bounds(DEPTHS_MM)
    rows = []
    curves = {}
    for path in result_files:
        energy = energy_from_name(path)
        modules = readDoseModules(path, module_prefix="doseLayer-")
        if len(modules) != len(DEPTHS_MM):
            raise RuntimeError(f"{path.name}: found {len(modules)} dose layers; expected 34")
        ordered = [modules[f"doseLayer-{layer}"] for layer in range(1, 35)]
        units = {item["unit"] for item in ordered}
        entries = {item["entries"] for item in ordered}
        if units != {"MeV/g"}:
            raise RuntimeError(f"{path.name}: unexpected units {units}")
        if entries != {5_000_000}:
            raise RuntimeError(f"{path.name}: unexpected entry counts {entries}")

        dose = np.array([item["dose"] for item in ordered])
        error = np.array([item["error"] for item in ordered])
        nonzero = np.array([item["non-zeros"] for item in ordered])
        pdd = 100 * dose / np.max(dose)
        relative_error = np.divide(
            100 * error, dose, out=np.full_like(error, np.nan), where=dose != 0
        )
        curves[energy] = {
            "dose": dose, "error": error, "pdd": pdd,
            "relative_error": relative_error, "nonzero": nonzero,
        }

        for layer, depth in enumerate(DEPTHS_MM, 1):
            index = layer - 1
            rows.append({
                "energy_MeV": energy,
                "layer": layer,
                "depth_mm": depth,
                "lower_mm": boundaries[index],
                "upper_mm": boundaries[index + 1],
                "thickness_mm": boundaries[index + 1] - boundaries[index],
                "dose_MeV_per_g_per_cm2": dose[index],
                "error_MeV_per_g_per_cm2": error[index],
                "relative_error_percent": relative_error[index],
                "nonzero_entries": nonzero[index],
                "primary_histories": ordered[index]["entries"],
                "normalized_PDD_percent": pdd[index],
                "measured_PDD_percent": measured_pdd[index],
                "residual_percentage_points": pdd[index] - measured_pdd[index],
            })

    output = run / "pilot_depth_dose.csv"
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    metrics = []
    electron_mask = DEPTHS_MM <= 40.0
    for energy, curve in sorted(curves.items()):
        pdd = curve["pdd"]
        residual = pdd - measured_pdd
        metrics.append({
            "energy_MeV": energy,
            "dmax_sample_mm": DEPTHS_MM[np.argmax(pdd)],
            "surface_PDD_percent": pdd[0],
            "R80_mm_linear_diagnostic": falling_crossing(DEPTHS_MM, pdd, 80),
            "R50_mm_linear_diagnostic": falling_crossing(DEPTHS_MM, pdd, 50),
            "RMSE_0_to_40mm_percentage_points": np.sqrt(np.mean(residual[electron_mask] ** 2)),
            "max_relative_error_0_to_40mm_percent": np.nanmax(curve["relative_error"][electron_mask]),
            "relative_error_60mm_percent": curve["relative_error"][-1],
            "nonzero_entries_60mm": int(curve["nonzero"][-1]),
        })

    metrics_output = run / "pilot_metrics.csv"
    with metrics_output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=metrics[0].keys())
        writer.writeheader()
        writer.writerows(metrics)

    colors = {5.0: "tab:blue", 5.5: "tab:orange", 6.0: "tab:green"}
    plt.figure(figsize=(7.2, 4.8))
    plt.plot(DEPTHS_MM, measured_pdd, "ko", ms=3.5, label="HUS measured")
    for energy, curve in sorted(curves.items()):
        # Scale the marginal GRAS errors by a fixed observed maximum.
        # A full ratio uncertainty needs inter-layer covariance, not recorded here.
        pdd_error = 100 * curve["error"] / np.max(curve["dose"])
        plt.errorbar(DEPTHS_MM, curve["pdd"], yerr=pdd_error,
                     fmt="o-", ms=3, lw=1.2, capsize=2.5, elinewidth=0.9,
                     color=colors[energy], label=f"{energy:.1f} MeV mono")
    plt.xlabel("Nominal measurement depth [mm]")
    plt.ylabel("Dose / maximum [%]")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.figtext(
        0.5, 0.015,
        "Bars: GRAS statistical errors scaled by a fixed maximum.\n"
        "Shared normalization uncertainty excluded; measured errors unavailable.",
        ha="center", va="bottom", fontsize=8,
    )
    plt.tight_layout(rect=(0, 0.09, 1, 1))
    plt.savefig(run / "PDD-pilot.png", dpi=180)
    plt.savefig(run / "PDD-pilot.pdf")
    plt.close()

    plt.figure(figsize=(7.2, 4.8))
    for energy, curve in sorted(curves.items()):
        plt.plot(DEPTHS_MM, curve["relative_error"], "o-", ms=3,
                     color=colors[energy], label=f"{energy:.1f} MeV")
    plt.axhline(1, color="black", ls="--", lw=1, label="1%")
    plt.xlabel("Nominal measurement depth [mm]")
    plt.ylabel("GRAS relative statistical error [%]")
    plt.ylim(bottom=0)
    plt.grid(which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(run / "PDD-pilot-relative-error.png", dpi=180)
    plt.savefig(run / "PDD-pilot-relative-error.pdf")
    plt.close()

    plt.figure(figsize=(7.2, 4.8))
    for energy, curve in sorted(curves.items()):
        plt.semilogy(DEPTHS_MM, curve["nonzero"], "o-", ms=3,
                     color=colors[energy], label=f"{energy:.1f} MeV")
    plt.xlabel("Nominal measurement depth [mm]")
    plt.ylabel("Histories with non-zero layer dose")
    plt.grid(which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(run / "PDD-pilot-nonzero.png", dpi=180)
    plt.savefig(run / "PDD-pilot-nonzero.pdf")
    plt.close()

    print(f"Wrote {output}")
    print(f"Wrote {metrics_output}")
    for metric in metrics:
        print(
            f"{metric['energy_MeV']:.1f} MeV: dmax={metric['dmax_sample_mm']:.1f} mm, "
            f"R50={metric['R50_mm_linear_diagnostic']:.2f} mm, "
            f"RMSE(0-40)={metric['RMSE_0_to_40mm_percentage_points']:.2f} pp, "
            f"error@60={metric['relative_error_60mm_percent']:.2f}%"
        )


if __name__ == "__main__":
    main()
