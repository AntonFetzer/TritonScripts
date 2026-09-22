"""Field-size and table series for CRDS at HUS, both instruments.

Four runs sharing the GDML, the VT01 + TSSS2600 detector set, the 1 um cut in
both fine-cut regions, the PDD spectrum and the 200 mm air gap. GRAS normalises
to one incident electron/cm2, so coefficients are directly comparable across
fields even though the runs differ in history count.

At 20 x 20 mm a ~5 MeV electron field is below lateral scatter equilibrium: the
lateral spread of the secondaries is comparable to the field half-width, so
dose on the axis is depressed relative to a broad field. The clinical
applicator was 36 x 36 cm. The 100 mm step measures what the specified 20 mm
field costs; the 200 mm step asks whether 100 mm had already saturated, which
decides whether the 100 mm coefficient can stand in for a clinical exposure.
The no-table run repeats 20 x 20 mm with the patient table absent, bounding the
error made by modelling a carbon-fibre-and-foam couch as solid PMMA.

Both instruments are carried through, because they do not have to agree: the
RadFET and the LED sit 10 mm apart at different depths, and displacement damage
weights the spectrum differently from ionising dose.

    python3 Plotting/CRDSHUSFieldSize.py

Aggregation is Dependencies.AggregateRun and Dependencies.TotalNID; nothing
here re-parses a GRAS file.

There is no source area anywhere in this analysis, and there does not need to
be. GRAS scales every result by fluence * source_surface / n_primaries before
writing it, so the coefficients below are already per unit fluence. The field
each run used is named by its folder and its beam macro.
"""

import csv
import math
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from Dependencies.AggregateRun import aggregateRun  # noqa: E402
from Dependencies.TotalNID import totalNID  # noqa: E402

BASE = Path("/scratch/work/fetzera1/GRAS/CRDS/CRDS-HUS")
RUNS = {
    "20 x 20 mm": {"folder": "PDD-Electron-CRDS-20mmField", "job": 20392947},
    "100 x 100 mm": {"folder": "PDD-Electron-CRDS-100mmField", "job": 20392948},
    "200 x 200 mm": {"folder": "PDD-Electron-CRDS-200mmField", "job": 20392949},
    "20 x 20 mm, no table": {"folder": "PDD-Electron-CRDS-noTable", "job": 20392950},
}
# Successive widenings, then the table control, each against its reference.
STEPS = [("20 x 20 mm", "100 x 100 mm"),
         ("100 x 100 mm", "200 x 200 mm"),
         ("20 x 20 mm, no table", "20 x 20 mm")]
TILES = [(0, "VT01_gox_0_PV", "gate oxide", "TID"),
         (1, "VT01_die_0_PV", "silicon die", "TID")]
MODULES = [("nidActive", "LED_active_0_PV", "active layer", "DDD"),
           ("nidSubstrate", "LED_die_0_PV", "GaAs substrate", "DDD")]
REFERENCE_FLUENCE = 2e12


def ratioWithError(numerator, numeratorError, denominator, denominatorError):
    """Ratio of two independent measurements and its absolute error."""
    ratio = numerator / denominator
    relative = math.sqrt((numeratorError / numerator) ** 2
                         + (denominatorError / denominator) ** 2)
    return ratio, ratio * relative


def main() -> None:
    data, rows = {}, []
    for label, spec in RUNS.items():
        run = BASE / spec["folder"]
        tid = aggregateRun(str(run / "Res"), expectedFiles=100, tileCount=len(TILES))
        ddd = totalNID(str(run / "Res"), expectedFiles=100, expectedModules=len(MODULES))
        entries = {}
        for index, volume, name, kind in TILES:
            entries[(kind, name)] = (float(tid["dose"][index]),
                                     float(tid["error"][index]),
                                     "kRad cm2/electron", volume,
                                     int(tid["entries"][index]),
                                     int(tid["non-zeros"][index]))
        for module, volume, name, kind in MODULES:
            entry = ddd[module]
            entries[(kind, name)] = (entry["nid"], entry["error"],
                                     entry["unit"] + " cm2/electron", volume,
                                     entry["entries"], entry["non-zeros"])
        data[label] = entries

        for (kind, name), (value, error, unit, volume, n, nz) in entries.items():
            rows.append({
                "run": label, "folder": spec["folder"], "slurm_job_id": spec["job"],
                "quantity": kind, "tally": name, "volume_name": volume,
                "simulated_primaries": n,
                "nonzero_entries": nz,
                "coefficient_per_electron_cm2": f"{value:.12e}",
                "statistical_error": f"{error:.12e}",
                "relative_error_percent": f"{100 * error / value:.6f}",
                "reference_fluence": f"{REFERENCE_FLUENCE:.6e}",
                "at_reference_fluence": f"{value * REFERENCE_FLUENCE:.6e}",
            })

    width = max(len(k) for k in RUNS) + 2
    print(f"\n{'run':<{width}}{'quantity':<10}{'tally':<16}{'coefficient':>16}"
          f"{'rel err':>10}{'at 2e12':>14}")
    for label in RUNS:
        for (kind, name), (value, error, unit, _, _, _) in data[label].items():
            at = value * REFERENCE_FLUENCE
            shown = f"{at:.3f} kRad" if kind == "TID" else f"{at:.4g} MeV/g"
            print(f"{label:<{width}}{kind:<10}{name:<16}{value:>16.6e}"
                  f"{100 * error / value:>9.2f}%{shown:>14}")

    print("\nsteps, each against the run named second")
    for numerator, denominator in STEPS:
        print(f"  {numerator}  vs  {denominator}")
        for key in data[numerator]:
            a, ea = data[numerator][key][0], data[numerator][key][1]
            b, eb = data[denominator][key][0], data[denominator][key][1]
            ratio, error = ratioWithError(a, ea, b, eb)
            sigma = abs(ratio - 1) / error if error else float("inf")
            print(f"    {key[0]:<4}{key[1]:<16}{100 * (ratio - 1):>+8.2f}%"
                  f" +- {100 * error:.2f}%   ({sigma:.1f} sigma)")

    output = BASE / "FieldSizeControl_CRDS-HUS.csv"
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {output}")


if __name__ == "__main__":
    main()