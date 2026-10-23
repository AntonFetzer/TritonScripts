"""Add first-irradiation MU-normalized doses to existing HUS coefficient tables.

The water calibration gives fluence per MU at the SSD100 reference plane. The
instrument simulations start a parallel beam at SSD100 and transport it through
20 cm of air to the device at SSD120; air scatter and energy loss are in that
transport, geometric divergence is not. The inverse-square factor from SSD100 to
SSD120 is therefore applied here, exactly once, and every consumer must use the
SSD120 fluences written below rather than scaling the SSD100 values itself.
"""
import csv
import json
import math
from pathlib import Path

BASE=Path("/scratch/work/fetzera1/GRAS/RadEx/RadEx-HUS")
CALIBRATION_SSD_CM=100
DEVICE_SSD_CM=120
INVERSE_SQUARE=(CALIBRATION_SSD_CM/DEVICE_SSD_CM)**2
DELIVERED_MU=86573

# Columns this script owns; removed before rewriting so stale values cannot survive.
MU_COLUMNS=("normalization_status","dose_material_model","source_SSD_cm","device_SSD_cm",
            "inverse_square_factor","first_irradiation_delivered_MU",
            "first_irradiation_source_fluence_electrons_cm2",
            "first_irradiation_fluence_SSD100_electrons_cm2",
            "first_irradiation_fluence_SSD120_electrons_cm2",
            "first_irradiation_fluence_stat_error_electrons_cm2",
            "first_irradiation_fluence_SSD120_stat_error_electrons_cm2",
            "dose_per_MU_kRad","dose_first_irradiation_kRad",
            "dose_first_irradiation_detector_stat_error_kRad",
            "dose_first_irradiation_normalization_stat_error_kRad",
            "dose_first_irradiation_combined_stat_error_kRad")

def main():
    # 50 x 50 mm scorers, jobs 20964818-20964820; the three thicknesses agree.
    calibration=json.loads((BASE/"MU-Fluence-Calibration/Calibration_baseline.json").read_text())
    result=calibration["result"]
    if result["source_SSD_cm"]!=CALIBRATION_SSD_CM:
        raise ValueError(f"Expected an SSD{CALIBRATION_SSD_CM} calibration, got SSD{result['source_SSD_cm']}")
    phi100=result["delivered_fluence_SSD100_electrons_cm2"]
    sigma100=result["delivered_fluence_stat_error_electrons_cm2"]
    phi=phi100*INVERSE_SQUARE
    sigma=sigma100*INVERSE_SQUARE
    result.update(device_SSD_cm=DEVICE_SSD_CM,inverse_square_factor_SSD100_to_SSD120=INVERSE_SQUARE,
                  fluence_SSD120_electrons_cm2_per_MU=result["fluence_SSD100_electrons_cm2_per_MU"]*INVERSE_SQUARE,
                  fluence_SSD120_stat_error_electrons_cm2_per_MU=result["fluence_stat_error_electrons_cm2_per_MU"]*INVERSE_SQUARE,
                  delivered_fluence_SSD120_electrons_cm2=phi,
                  delivered_fluence_SSD120_stat_error_electrons_cm2=sigma)
    calibration["exposure_scope"]="first irradiation only; later MU inferred separately from timing"
    calibration["status"]=("provisional: 0.25/0.5/1.0 mm scorers consistent (chi2 3.9 for 2 dof); "
                           "SSD100 calibration geometry, Plastic Water equivalence and the "
                           "100 cm virtual source distance assumed")
    calibration["device_normalization"]=("Use the SSD120 fluences: inverse-square factor applied once here; "
                                         "the parallel instrument beam already includes the 20 cm air transport")
    calibration["relative_normalization_stat_error"]=sigma/phi
    (BASE/"Results"/"HUS-MU-calibration.json").write_text(json.dumps(calibration,indent=2)+"\n")
    # The PDD module run is the current prediction and feeds the shared table; the
    # mono-energetic companion gets its own table, suffixed with the run name.
    for folder,model,table in [("PDD-Electron","gate_oxide_module","RadEx-HUS-MU-normalized-dose-results.csv"),
                               ("Mono6MeV-Electron","gate_oxide_module","RadEx-HUS-MU-normalized-dose-results-Mono6MeV.csv")]:
        if not (BASE/folder/"dose_results_full_precision.csv").exists():
            continue
        path=BASE/folder/"dose_results_full_precision.csv"
        rows=[{k:v for k,v in r.items() if k not in MU_COLUMNS} for r in csv.DictReader(path.open())]
        for r in rows:
            c=float(r["dose_coefficient_kRad_cm2_per_electron"])
            e=float(r["statistical_error_kRad_cm2_per_electron"])
            r.update(normalization_status="provisional",dose_material_model=model,
                     source_SSD_cm=CALIBRATION_SSD_CM,device_SSD_cm=DEVICE_SSD_CM,
                     inverse_square_factor=INVERSE_SQUARE,
                     first_irradiation_delivered_MU=DELIVERED_MU,
                     first_irradiation_fluence_SSD100_electrons_cm2=phi100,
                     first_irradiation_fluence_SSD120_electrons_cm2=phi,
                     first_irradiation_fluence_SSD120_stat_error_electrons_cm2=sigma,
                     dose_per_MU_kRad=c*phi/DELIVERED_MU,
                     dose_first_irradiation_kRad=c*phi,
                     dose_first_irradiation_detector_stat_error_kRad=e*phi,
                     dose_first_irradiation_normalization_stat_error_kRad=c*sigma,
                     dose_first_irradiation_combined_stat_error_kRad=math.hypot(e*phi,c*sigma))
        with path.open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
        if model=="gate_oxide_module":
            with (BASE/"Results"/table).open("w",newline="") as f:
                w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
        print(folder,model,[(r["channel"],round(r["dose_first_irradiation_kRad"],3)) for r in rows])
    print(f"Device fluence {phi:.6e} electrons/cm2 = SSD100 {phi100:.6e} x inverse-square {INVERSE_SQUARE:.6f}.")
    print("Reference coefficients and 2e12 columns preserved. Common normalization errors are correlated across channels.")

if __name__=="__main__":
    main()
