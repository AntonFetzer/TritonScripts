"""Field-size series for CRDS at HUS: 20 x 20, 100 x 100 and 200 x 200 mm.

All three runs share the GDML, the VT01 package detector, the RadFET-local 1 um
cut, the PDD spectrum, the 200 mm air path and 1e9 histories. Only the source
rectangle differs. GRAS normalises to one incident electron/cm2, so the
coefficients are directly comparable across fields.

At 20 x 20 mm a ~5 MeV electron field is below lateral scatter equilibrium: the
lateral spread of the secondaries is comparable to the field half-width, so dose
on the axis is depressed relative to a broad field. The clinical applicator was
36 x 36 cm. The 100 mm step measures how much the specified 20 mm field costs;
the 200 mm step asks whether 100 mm had already saturated, which decides whether
the 100 mm coefficient can stand in for a clinical-applicator exposure.

    python3 Plotting/CRDSHUSFieldSize.py

Aggregation is Dependencies.AggregateRun.aggregateRun, which wraps
Dependencies.TotalDose.totalDose; nothing here re-parses a GRAS file.
"""

import argparse
import csv
import math
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.AggregateRun import aggregateRun  # noqa: E402

BASE = Path("/scratch/work/fetzera1/GRAS/CRDS/CRDS-HUS")
RUNS = {
    "20 x 20 mm": {"folder": "PDD-Electron-CRDS-NewRadFET-Cuts",
                   "job": 20332204, "area_cm2": 4.0, "files": 100},
    "100 x 100 mm": {"folder": "PDD-Electron-CRDS-NewRadFET-Cuts-100mmField",
                     "job": 20335084, "area_cm2": 100.0, "files": 100},
    "200 x 200 mm": {"folder": "PDD-Electron-CRDS-NewRadFET-Cuts-200mmField",
                     "job": 20335552, "area_cm2": 400.0, "files": 100},
}
# Successive widenings, each compared against the field one step narrower.
STEPS = [("20 x 20 mm", "100 x 100 mm"), ("100 x 100 mm", "200 x 200 mm")]
TILES = [(0, "VT01_gox_0_PV", "gate oxide"), (1, "VT01_die_0_PV", "silicon die")]
REFERENCE_FLUENCE = 2e12

# RadEx-HUS Ch9, the exposed 0 mm channel: same field, same distance, same
# detector model, 61.435 kRad on its gate oxide at the same reference fluence.
# It sits behind the instrument's aluminium, so this is context, not a target.
RADEX_HUS_CH9_OXIDE_KRAD = 61.435


def ratioWithError(numerator, numeratorError, denominator, denominatorError):
    """Ratio of two independent measurements and its absolute error.

    Relative errors combine in quadrature, which is the repository convention
    for independent uncertainties.
    """
    ratio = numerator / denominator
    return ratio, ratio * math.hypot(numeratorError / numerator,
                                     denominatorError / denominator)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path,
                        default=BASE / "FieldSizeControl_CRDS-HUS.csv")
    args = parser.parse_args()

    data = {}
    for label, spec in RUNS.items():
        data[label] = aggregateRun(
            BASE / spec["folder"] / "Res",
            expectedFiles=spec["files"],
            tileCount=len(TILES),
        )

    print(f"\n{'field':<14s}{'layer':<12s}{'coefficient':>16s}{'rel err':>9s}"
          f"{'hit fraction':>14s}{'dose at 2e12':>14s}")
    rows = []
    for label, spec in RUNS.items():
        results = data[label]
        for index, volume, layer in TILES:
            dose = float(results["dose"][index])
            error = float(results["error"][index])
            entries = float(results["entries"][index])
            nonZeros = float(results["non-zeros"][index])
            relative = float(results["relative_error_percent"][index])
            hitFraction = float(results["hit_fraction"][index])
            print(f"{label:<14s}{layer:<12s}{dose:>16.6e}"
                  f"{relative:>8.3f}%{hitFraction:>14.4e}"
                  f"{dose * REFERENCE_FLUENCE:>10.3f} kRad")
            rows.append({
                "field_size_mm": label, "slurm_job_id": spec["job"],
                "source_area_cm2": f"{spec['area_cm2']:.1f}",
                "tile_index": index, "volume_name": volume, "layer": layer,
                "simulated_primaries": int(entries),
                "dose_coefficient_kRad_cm2_per_electron": f"{dose:.12e}",
                "statistical_error_kRad_cm2_per_electron": f"{error:.12e}",
                "relative_error_percent": f"{relative:.8f}",
                "nonzero_entries": int(nonZeros),
                "hit_fraction": f"{hitFraction:.9e}",
                "dose_at_reference_kRad": f"{dose * REFERENCE_FLUENCE:.9e}",
            })

    print("\nField-size steps, each relative to the next narrower field")
    for narrow, wide in STEPS:
        print(f"  {narrow} -> {wide}")
        small, large = data[narrow], data[wide]
        for index, _volume, layer in TILES:
            ratio, error = ratioWithError(
                float(large["dose"][index]), float(large["error"][index]),
                float(small["dose"][index]), float(small["error"][index]))
            print(f"    {layer:<12s}{100 * (ratio - 1):+7.2f} +- {100 * error:5.2f} %"
                  f"   ({(ratio - 1) / error:+6.1f} sigma)   deficit at "
                  f"{narrow.split()[0]} mm: {100 * (1 - 1 / ratio):+.2f} %")

    widest = "200 x 200 mm"
    print(f"\nEvery field against the widest ({widest})")
    for label in RUNS:
        if label == widest:
            continue
        for index, _volume, layer in TILES:
            ratio, error = ratioWithError(
                float(data[label]["dose"][index]),
                float(data[label]["error"][index]),
                float(data[widest]["dose"][index]),
                float(data[widest]["error"][index]))
            print(f"  {label:<14s}{layer:<12s}{100 * (ratio - 1):+7.2f} "
                  f"+- {100 * error:5.2f} % of the widest-field dose")

    for label in RUNS:
        results = data[label]
        ratio, error = ratioWithError(
            float(results["dose"][0]), float(results["error"][0]),
            float(results["dose"][1]), float(results["error"][1]))
        print(f"\n{label}: oxide/die = {ratio:.4f} +- {error:.4f}")
        oxide = float(results["dose"][0]) * REFERENCE_FLUENCE
        print(f"{'':{len(label)}s}  vs RadEx-HUS Ch9 oxide "
              f"{RADEX_HUS_CH9_OXIDE_KRAD} kRad: "
              f"{100 * (oxide / RADEX_HUS_CH9_OXIDE_KRAD - 1):+.2f} %")

    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
