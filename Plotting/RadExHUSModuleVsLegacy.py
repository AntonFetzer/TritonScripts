"""Compare the VT01 package module against the legacy bulk-silicon RadFET
for the HUS PDD electron exposure.

Both runs use the same measured HUS geometry, the same PDD-derived source at
SSD100, the same beam macro and the same 3.98e9 total histories. The only
differences are the RadFET representation and the 1 um production cut confined
to the RadFET region, so the per-channel ratio isolates the detector model.

The module scores 24 volumes: tiles 0-11 are the gate oxides, tiles 12-23 the
silicon dies of the same twelve packages. The legacy model scores one bulk
silicon tile per channel, so the legacy number is compared against both.

Aggregation is Dependencies.TotalDose.totalDose in both cases, which applies the
entry-weighted pooling, the rad-to-kRad conversion and the Birge-ratio
consistency check. Nothing here re-parses a GRAS file.

The module coefficient table is rewritten without the first-irradiation MU
columns; rerun Calibration.PropagateHUSMU afterwards to restore them.
"""

import argparse
import csv
import os
from pathlib import Path
import sys

import numpy as np
from uncertainties import ufloat

PYTHON_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_ROOT))
from Dependencies.TotalDose import totalDose  # noqa: E402

HUS = "/scratch/work/fetzera1/GRAS/RadEx/RadEx-HUS"
LEGACY = os.path.join(HUS, "PDD-Electron-MeasuredGeometry", "Res")
MODULE = os.path.join(HUS, "PDD-Electron-MeasuredGeometry-NewRadFET-Cuts", "Res")

# Channel labels as used throughout the campaign records.
LABEL = {
    0: "100% PE", 1: "25% Pb / 75% PE", 2: "50% Pb / 50% PE",
    3: "75% Pb / 25% PE", 4: "100% Pb", 5: "6 mm Al", 6: "4 mm Al",
    7: "2 mm Al", 8: "1 mm Al", 9: "exposed 0 mm",
    10: "PCB back behind Ch9", 11: "PCB back centre",
}
DEVICE = {c: ("VT05" if c in (7, 8, 9) else "VT01") for c in range(12)}
FRONT_VT01 = [0, 1, 2, 3, 4, 5, 6]
FLUENCE = 2e12  # electrons/cm2, the campaign reference normalisation


def collect(path, count):
    r = totalDose(path)
    if len(r["dose"]) != count:
        raise RuntimeError(f"{path}: expected {count} tiles, got {len(r['dose'])}")
    return [ufloat(float(d), float(e)) for d, e in zip(r["dose"], r["error"])]


def order(values, channels):
    return sorted(channels, key=lambda c: values[c].n)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=os.path.join(
        HUS, "PDD-Electron-MeasuredGeometry-NewRadFET-Cuts",
        "ModuleVsLegacy_PDD-Electron.csv"))
    parser.add_argument("--coefficients", default=os.path.join(
        HUS, "PDD-Electron-MeasuredGeometry-NewRadFET-Cuts",
        "dose_results_full_precision.csv"))
    args = parser.parse_args()

    legacy = collect(LEGACY, 12)
    module = collect(MODULE, 24)
    oxide = module[:12]
    die = module[12:]

    print(f"{'ch':>3s} {'sample':<22s}{'legacy':>12s}{'oxide':>12s}"
          f"{'ox/leg %':>10s}{'die':>12s}{'ox/die %':>10s}")
    rows = []
    for c in range(12):
        lg, ox, di = legacy[c] * FLUENCE, oxide[c] * FLUENCE, die[c] * FLUENCE
        r_leg = 100.0 * (ox / lg - 1.0)
        r_die = 100.0 * (ox / di - 1.0)
        print(f"{c:>3d} {LABEL[c]:<22s}{lg.n:>12.3f}{ox.n:>12.3f}"
              f"{r_leg.n:>+8.2f}{'':2s}{di.n:>12.3f}{r_die.n:>+8.2f}")
        rows.append(dict(
            channel=c, radfet_type=DEVICE[c], shielding_or_location=LABEL[c],
            legacy_dose_kRad=lg.n, legacy_error_kRad=lg.s,
            oxide_dose_kRad=ox.n, oxide_error_kRad=ox.s,
            die_dose_kRad=di.n, die_error_kRad=di.s,
            oxide_over_legacy=(ox / lg).n, oxide_over_legacy_error=(ox / lg).s,
            oxide_over_die=(ox / di).n, oxide_over_die_error=(ox / di).s,
        ))

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {args.output}")

    # Coefficient table in the campaign's comparison format so the existing
    # ordering analysis can consume the module result unchanged.
    with open(args.coefficients, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "energy_mev", "particle", "tile_index", "channel", "radfet_type",
            "shielding_or_location", "dose_coefficient_kRad_cm2_per_electron",
            "statistical_error_kRad_cm2_per_electron", "relative_error_percent",
            "dose_at_2e12_electrons_cm2_kRad",
            "statistical_error_at_2e12_electrons_cm2_kRad",
            "die_dose_coefficient_kRad_cm2_per_electron"])
        for c in range(12):
            ox = oxide[c]
            writer.writerow([
                6, "electron", c, c, DEVICE[c], LABEL[c],
                repr(ox.n), repr(ox.s), repr(100.0 * ox.s / ox.n),
                repr(ox.n * FLUENCE), repr(ox.s * FLUENCE), repr(die[c].n)])
    print(f"Wrote {args.coefficients}")

    print("\nFront VT01 ordering, low to high")
    for name, values in (("legacy", legacy), ("oxide", oxide), ("die", die)):
        seq = order(values, FRONT_VT01)
        print(f"  {name:<8s}" + " < ".join(f"Ch{c}" for c in seq))

    print("\nFront VT01 spread (max/min)")
    for name, values in (("legacy", legacy), ("oxide", oxide), ("die", die)):
        v = [values[c].n for c in FRONT_VT01]
        print(f"  {name:<8s}{max(v) / min(v):.3f}")

    print("\nBack-mounted contrast, Ch/Ch0")
    for name, values in (("legacy", legacy), ("oxide", oxide), ("die", die)):
        print(f"  {name:<8s}Ch10/Ch0 {values[10].n / values[0].n:.4f}"
              f"   Ch11/Ch0 {values[11].n / values[0].n:.4f}")

    front_shift = np.mean([100.0 * (oxide[c].n / legacy[c].n - 1.0)
                           for c in FRONT_VT01])
    print(f"\nMean oxide/legacy shift over front VT01 Ch0-6: {front_shift:+.2f} %")


if __name__ == "__main__":
    main()
