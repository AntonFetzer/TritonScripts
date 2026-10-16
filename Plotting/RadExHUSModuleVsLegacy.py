"""Tabulate the HUS VT01 package-module run and write its coefficient table.

The module scores 24 volumes: tiles 0-11 are the gate oxides, tiles 12-23 the
silicon dies of the same twelve packages. The script writes a per-channel table
of oxide and die doses at the 2e12 electrons/cm2 campaign reference, and the
coefficient table in the campaign's comparison format that ElectronOrdering.py
and Calibration.PropagateHUSMU consume. The coefficient table is rewritten
without the first-irradiation MU columns; rerun Calibration.PropagateHUSMU
afterwards to restore them.

Two optional comparisons extend the table channel by channel:

    --legacy            a run with one bulk-silicon tile per channel (12 tiles),
                        compared against the oxide and the die
    --reference-module  a second module run (24 tiles), compared oxide to oxide

The default module run is the milled-absorber run, and the channel labels give
the aluminium as milled (Exterior_Top_Plate.step).

Aggregation is Dependencies.TotalDose.totalDose throughout, which applies the
entry-weighted pooling, the rad-to-kRad conversion and the Birge-ratio
consistency check. Nothing here re-parses a GRAS file.
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
RUN = os.path.join(HUS, "PDD-Electron-MeasuredGeometry-NewRadFET-Cuts-MilledAl")

# Channel labels as used throughout the campaign records, aluminium as milled.
LABEL = {
    0: "100% PE", 1: "25% Pb / 75% PE", 2: "50% Pb / 50% PE",
    3: "75% Pb / 25% PE", 4: "100% Pb", 5: "6 mm Al", 6: "5 mm Al",
    7: "4 mm Al", 8: "2 mm Al", 9: "exposed 0 mm",
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
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--module", default=os.path.join(RUN, "Res"),
                        help="Res directory of the module run")
    parser.add_argument("--output", default=os.path.join(RUN, "ModuleVsLegacy_PDD-Electron.csv"))
    parser.add_argument("--coefficients", default=os.path.join(RUN, "dose_results_full_precision.csv"))
    parser.add_argument("--legacy", default=None,
                        help="Res directory of a bulk-silicon run to compare against")
    parser.add_argument("--reference-module", default=None,
                        help="Res directory of a second module run to compare against")
    args = parser.parse_args()

    module = collect(args.module, 24)
    oxide = module[:12]
    die = module[12:]
    legacy = collect(args.legacy, 12) if args.legacy else None
    reference = collect(args.reference_module, 24)[:12] if args.reference_module else None

    print(f"{'ch':>3s} {'sample':<22s}{'oxide':>12s}{'die':>12s}{'ox/die %':>10s}"
          + (f"{'legacy':>12s}{'ox/leg %':>10s}" if legacy else ""))
    rows = []
    for c in range(12):
        ox, di = oxide[c] * FLUENCE, die[c] * FLUENCE
        row = dict(
            channel=c, radfet_type=DEVICE[c], shielding_or_location=LABEL[c],
            oxide_dose_kRad=ox.n, oxide_error_kRad=ox.s,
            die_dose_kRad=di.n, die_error_kRad=di.s,
            oxide_over_die=(ox / di).n, oxide_over_die_error=(ox / di).s,
        )
        line = f"{c:>3d} {LABEL[c]:<22s}{ox.n:>12.3f}{di.n:>12.3f}{100.0 * (ox / di - 1.0).n:>+8.2f}  "
        if legacy:
            lg = legacy[c] * FLUENCE
            row.update(legacy_dose_kRad=lg.n, legacy_error_kRad=lg.s,
                       oxide_over_legacy=(ox / lg).n, oxide_over_legacy_error=(ox / lg).s)
            line += f"{lg.n:>12.3f}{100.0 * (ox / lg - 1.0).n:>+8.2f}"
        if reference:
            rf = reference[c] * FLUENCE
            row.update(reference_oxide_dose_kRad=rf.n, reference_oxide_error_kRad=rf.s,
                       oxide_over_reference=(ox / rf).n, oxide_over_reference_error=(ox / rf).s)
        rows.append(row)
        print(line)

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

    series = [("oxide", oxide), ("die", die)] + ([("legacy", legacy)] if legacy else [])
    print("\nFront VT01 ordering, low to high")
    for name, values in series:
        print(f"  {name:<8s}" + " < ".join(f"Ch{c}" for c in order(values, FRONT_VT01)))

    print("\nFront VT01 spread (max/min)")
    for name, values in series:
        v = [values[c].n for c in FRONT_VT01]
        print(f"  {name:<8s}{max(v) / min(v):.3f}")

    print("\nBack-mounted contrast, Ch/Ch0")
    for name, values in series:
        print(f"  {name:<8s}Ch10/Ch0 {values[10].n / values[0].n:.4f}"
              f"   Ch11/Ch0 {values[11].n / values[0].n:.4f}")

    if legacy:
        front_shift = np.mean([100.0 * (oxide[c].n / legacy[c].n - 1.0) for c in FRONT_VT01])
        print(f"\nMean oxide/legacy shift over front VT01 Ch0-6: {front_shift:+.2f} %")

    if reference:
        print("\nModule run over reference module, gate oxide")
        for c in range(12):
            q = oxide[c] / reference[c]
            print(f"  Ch{c:<3d}{LABEL[c]:<22s}{100.0 * (q.n - 1.0):+7.2f} +- {100.0 * q.s:.2f} %")
        print("  reference order " + " < ".join(f"Ch{c}" for c in order(reference, FRONT_VT01)))


if __name__ == "__main__":
    main()
