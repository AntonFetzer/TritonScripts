"""Figure for the HUS PDD detector-model comparison.

Top panel: dose per channel for the legacy bulk-silicon RadFET, the module gate
oxide and the module silicon die, at the 2e12 electrons/cm2 campaign reference.
Bottom panel: the module-oxide change relative to the legacy model, which is the
quantity that matters for the ordering test because it is not a common scale
factor across channels.

Reads the comparison table written by Plotting.RadExHUSModuleVsLegacy; it does not
re-aggregate GRAS output.
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

FOLDER = ("/scratch/work/fetzera1/GRAS/RadEx/RadEx-HUS/"
          "PDD-Electron-MeasuredGeometry-NewRadFET-Cuts")
TABLE = os.path.join(FOLDER, "ModuleVsLegacy_PDD-Electron.csv")
OUTPUT = os.path.join(FOLDER, "ModuleVsLegacy_PDD-Electron.pdf")

FRONT = list(range(10))


def main():
    data = pd.read_csv(TABLE)
    channels = data["channel"].to_numpy()

    figure, (upper, lower) = plt.subplots(
        2, 1, figsize=(9, 7), sharex=True,
        gridspec_kw={"height_ratios": [2, 1], "hspace": 0.08})

    upper.errorbar(channels - 0.16, data["legacy_dose_kRad"],
                   yerr=data["legacy_error_kRad"], fmt="s", capsize=3,
                   markersize=5, label="legacy bulk silicon")
    upper.errorbar(channels, data["oxide_dose_kRad"],
                   yerr=data["oxide_error_kRad"], fmt="o", capsize=3,
                   markersize=5, label="VT01 module, gate oxide")
    upper.errorbar(channels + 0.16, data["die_dose_kRad"],
                   yerr=data["die_error_kRad"], fmt="^", capsize=3,
                   markersize=5, label="VT01 module, silicon die")
    upper.set_ylabel("Dose at 2e12 electrons/cm$^2$ [kRad]")
    upper.set_yscale("log")
    upper.grid(which="both", alpha=0.3)
    upper.legend()
    upper.set_title("HUS PDD electrons: RadFET detector model comparison\n"
                    "same geometry, source and history count; job 20315277 vs 20256114",
                    fontsize=10)

    percent = 100.0 * (data["oxide_over_legacy"] - 1.0)
    percent_error = 100.0 * data["oxide_over_legacy_error"]
    colours = ["tab:blue" if c in FRONT else "tab:red" for c in channels]
    lower.bar(channels, percent, yerr=percent_error, color=colours,
              capsize=3, width=0.6)
    lower.axhline(0.0, color="k", linewidth=0.8)
    lower.set_ylabel("oxide / legacy $-$ 1 [%]")
    lower.set_xlabel("Channel")
    lower.set_xticks(channels)
    lower.grid(axis="y", alpha=0.3)
    lower.text(0.99, 0.05, "red: PCB-back mounted", transform=lower.transAxes,
               ha="right", va="bottom", fontsize=8, color="tab:red")

    figure.savefig(OUTPUT, format="pdf", bbox_inches="tight")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
