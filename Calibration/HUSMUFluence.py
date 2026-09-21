"""Convert HUS water dose to SSD100 electron fluence per MU."""
import argparse
import csv
import json
import math
from pathlib import Path
from Read.ReadDose import readDoseModules
from Dependencies.PoolDoseModules import poolDoseModules

RUN = Path("/scratch/work/fetzera1/GRAS/RadEx/RadEx-HUS/MU-Fluence-Calibration")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", choices=("baseline", "thin", "thick"), default="baseline")
    parser.add_argument("--results", type=Path)
    parser.add_argument("--filename-contains", default="")
    parser.add_argument("--expected-files", type=int, default=5)
    parser.add_argument("--histories-per-file", type=int, default=5_000_000)
    args = parser.parse_args()
    results = args.results or RUN / ("Res" if args.variant == "baseline" else "Res-" + args.variant)
    files = sorted(p for p in results.glob("*.csv") if args.filename_contains in p.name)
    if len(files) != args.expected_files:
        raise ValueError(f"Expected {args.expected_files} files; found {len(files)} in {results}")
    for path in files:
        modules = readDoseModules(path)
        if set(modules) != {"doseWater13mm"}:
            raise ValueError(f"Unexpected scoring modules: {path}")
        value = modules["doseWater13mm"]
        if value["unit"] != "MeV/g" or value["entries"] != args.histories_per_file:
            raise ValueError(f"Unexpected units/history count: {path}")
    pooled = poolDoseModules(files, module_prefix="doseWater13mm")["doseWater13mm"]
    c = pooled["dose"] * 1.602176634e-10
    error = pooled["error"] * 1.602176634e-10
    if not math.isfinite(c) or not math.isfinite(error) or c <= 0 or pooled["inconsistent"]:
        raise ValueError(f"Invalid/inconsistent coefficient: {pooled}")
    phi = (1.025428 / 100) / c
    phi_error = phi * error / c
    row = dict(variant=args.variant, scorer_thickness_mm={"baseline":0.5,"thin":0.25,"thick":1.0}[args.variant],
               source_SSD_cm=100, histories=pooled["entries"],
               Cw_Gy_cm2_per_electron=c, Cw_stat_error_Gy_cm2_per_electron=error,
               relative_stat_error_percent=100*error/c, birge_ratio=pooled["birge_ratio"],
               fluence_SSD100_electrons_cm2_per_MU=phi,
               fluence_stat_error_electrons_cm2_per_MU=phi_error,
               delivered_MU=86573.0, delivered_fluence_SSD100_electrons_cm2=86573*phi,
               delivered_fluence_stat_error_electrons_cm2=86573*phi_error)
    output = results.parent / ("Calibration_" + args.variant)
    with output.with_suffix(".csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=row.keys())
        writer.writeheader(); writer.writerow(row)
    report = dict(result=row, files=[str(p) for p in files],
                  delivered_MU_evidence="User witnessed machine stop at 86573 MU; photo predates irradiation",
                  assumptions=["Calibration SSD100 inferred", "13 mm Plastic Water treated as 13 mm liquid water", "Parallel reconstructed electron-only source"],
                  uncertainty="Monte Carlo only; calibration and model uncertainties not supplied")
    output.with_suffix(".json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
