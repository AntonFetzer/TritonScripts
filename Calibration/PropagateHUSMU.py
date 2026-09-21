"""Add first-irradiation MU-normalized doses to existing HUS coefficient tables."""
import csv
import json
import math
from pathlib import Path

BASE=Path("/scratch/work/fetzera1/GRAS/RadEx/RadEx-HUS")

def main():
    calibration=json.loads((BASE/"MU-Fluence-Calibration/Completed-20370158/Calibration_baseline.json").read_text())
    result=calibration["result"]
    phi=result["delivered_fluence_SSD100_electrons_cm2"]
    sigma=result["delivered_fluence_stat_error_electrons_cm2"]
    calibration["exposure_scope"]="first irradiation only; later MU inferred separately from timing"
    calibration["status"]="provisional: scorer-thickness discrepancy unresolved"
    calibration["relative_normalization_stat_error"]=sigma/phi
    (BASE/"HUS-MU-calibration.json").write_text(json.dumps(calibration,indent=2)+"\n")
    for folder,model in [("PDD-Electron-MeasuredGeometry","legacy_bulk_silicon"),("PDD-Electron-MeasuredGeometry-NewRadFET-Cuts","gate_oxide_module")]:
        path=BASE/folder/"dose_results_full_precision.csv"
        rows=list(csv.DictReader(path.open()))
        for r in rows:
            c=float(r["dose_coefficient_kRad_cm2_per_electron"])
            e=float(r["statistical_error_kRad_cm2_per_electron"])
            r.update(normalization_status="provisional",dose_material_model=model,source_SSD_cm=100,
                     first_irradiation_delivered_MU=86573,
                     first_irradiation_source_fluence_electrons_cm2=phi,
                     first_irradiation_fluence_stat_error_electrons_cm2=sigma,
                     dose_per_MU_kRad=c*phi/86573,
                     dose_first_irradiation_kRad=c*phi,
                     dose_first_irradiation_detector_stat_error_kRad=e*phi,
                     dose_first_irradiation_normalization_stat_error_kRad=c*sigma,
                     dose_first_irradiation_combined_stat_error_kRad=math.hypot(e*phi,c*sigma))
        with path.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
        if model=="gate_oxide_module":
            with (BASE/"RadEx-HUS-MU-normalized-dose-results.csv").open("w",newline="") as f:
                w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
        print(model,[(r["channel"],round(r["dose_first_irradiation_kRad"],3)) for r in rows])
    print("Reference coefficients and 2e12 columns preserved. Common normalization errors are correlated across channels.")

if __name__=="__main__":
    main()
