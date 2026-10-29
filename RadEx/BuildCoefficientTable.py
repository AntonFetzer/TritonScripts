"""Write dose_results_full_precision.csv, the per-channel dose coefficient table of one run.

    python3 -m RadEx.BuildCoefficientTable <run folder> --campaign <json>

The run scores the campaign's layers (for the VT01 package module, gate oxides in
tiles 0-11 and silicon dies in 12-23); the table holds the first layer's coefficients
and the second layer's beside them. The channel labels, device types, particle,
energy label and reference fluence come from the campaign file (RadEx.KernelFit).
The format is the one ElectronOrdering.py and Calibration.PropagateHUSMU consume; it
is written without the first-irradiation MU columns, which Calibration.PropagateHUSMU
adds afterwards.

Aggregation is Dependencies.TotalDose.totalDose; nothing here re-parses a GRAS file.
"""

import argparse
import csv
from pathlib import Path

from Dependencies.TotalDose import totalDose
from RadEx import KernelFit as kf


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run", type=Path)
    parser.add_argument("--campaign", type=Path, required=True)
    args = parser.parse_args()
    campaign = kf.load_campaign(args.campaign)
    run = args.run.resolve()
    res = run / "Res" if (run / "Res").is_dir() else run
    r = totalDose(str(res))
    dose = [float(d) for d in r["dose"]]
    err = [float(e) for e in r["error"]]
    layers = sorted(campaign["layers"], key=campaign["layers"].get)
    expected = len(layers) * len(campaign["channels"])
    if len(dose) != expected:
        raise RuntimeError(f"{res}: expected {expected} tiles, got {len(dose)}")
    first, second = kf.layer_slice(campaign, layers[0]), kf.layer_slice(campaign, layers[1])
    oxide, oxide_err, die = dose[first], err[first], dose[second]
    fluence = campaign["reference_fluence"]
    tag = f"{fluence:.0e}".replace("+", "")

    out = run / "dose_results_full_precision.csv"
    with out.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "energy_mev", "particle", "tile_index", "channel", "radfet_type",
            "shielding_or_location", f"dose_coefficient_kRad_cm2_per_{campaign['particle']}",
            f"statistical_error_kRad_cm2_per_{campaign['particle']}", "relative_error_percent",
            f"dose_at_{tag}_{campaign['particle']}s_cm2_kRad",
            f"statistical_error_at_{tag}_{campaign['particle']}s_cm2_kRad",
            f"die_dose_coefficient_kRad_cm2_per_{campaign['particle']}"])
        for c in campaign["channels"]:
            writer.writerow([
                campaign["coefficient_table_energy_MeV"], campaign["particle"], c, c, campaign["device"][c],
                campaign["labels"][c], repr(oxide[c]), repr(oxide_err[c]), repr(100.0 * oxide_err[c] / oxide[c]),
                repr(oxide[c] * fluence), repr(oxide_err[c] * fluence), repr(die[c])])
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
