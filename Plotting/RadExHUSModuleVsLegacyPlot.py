"""Figure for the HUS PDD module run and its optional comparisons.

Top panel: gate-oxide and silicon-die dose per channel at the 2e12 electrons/cm2
campaign reference, with the legacy bulk-silicon and reference-module doses when
the table carries them. One further panel per comparison shows the oxide change
relative to that run, which is what matters for the ordering test because it is
not a common scale factor across channels.

Reads the table written by Plotting.RadExHUSModuleVsLegacy; it does not
re-aggregate GRAS output.
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

FOLDER = ("/scratch/work/fetzera1/GRAS/RadEx/RadEx-HUS/"
          "PDD-Electron-MeasuredGeometry-NewRadFET-Cuts-MilledAl")
TITLE = "milled aluminium absorbers, job 20996163"

FRONT = list(range(10))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--folder", default=FOLDER,
                        help="run folder holding ModuleVsLegacy_PDD-Electron.csv")
    parser.add_argument("--title", default=TITLE, help="second title line")
    args = parser.parse_args()
    table = os.path.join(args.folder, "ModuleVsLegacy_PDD-Electron.csv")
    output = os.path.join(args.folder, "ModuleVsLegacy_PDD-Electron.pdf")

    data = pd.read_csv(table)
    channels = data["channel"].to_numpy()
    comparisons = [(column, label) for column, label in (
        ("legacy", "oxide / legacy $-$ 1 [%]"),
        ("reference", "oxide / reference $-$ 1 [%]"))
        if f"oxide_over_{column}" in data.columns]

    figure, axes = plt.subplots(
        1 + len(comparisons), 1, figsize=(9, 4.7 + 2.4 * len(comparisons)), sharex=True,
        gridspec_kw={"height_ratios": [2] + [1] * len(comparisons), "hspace": 0.08},
        squeeze=False)
    upper = axes[0, 0]

    upper.errorbar(channels, data["oxide_dose_kRad"],
                   yerr=data["oxide_error_kRad"], fmt="o", capsize=3,
                   markersize=5, label="VT01 module, gate oxide")
    upper.errorbar(channels + 0.16, data["die_dose_kRad"],
                   yerr=data["die_error_kRad"], fmt="^", capsize=3,
                   markersize=5, label="VT01 module, silicon die")
    if "legacy_dose_kRad" in data.columns:
        upper.errorbar(channels - 0.16, data["legacy_dose_kRad"],
                       yerr=data["legacy_error_kRad"], fmt="s", capsize=3,
                       markersize=5, label="legacy bulk silicon")
    if "reference_oxide_dose_kRad" in data.columns:
        upper.errorbar(channels + 0.32, data["reference_oxide_dose_kRad"],
                       yerr=data["reference_oxide_error_kRad"], fmt="x", capsize=3,
                       markersize=5, label="reference module, gate oxide")
    upper.set_ylabel("Dose at 2e12 electrons/cm$^2$ [kRad]")
    upper.set_yscale("log")
    upper.grid(which="both", alpha=0.3)
    upper.legend()
    upper.set_title("HUS PDD electrons, VT01 package module\n" + args.title, fontsize=10)

    colours = ["tab:blue" if c in FRONT else "tab:red" for c in channels]
    for (column, label), ax in zip(comparisons, axes[1:, 0]):
        change = 100.0 * (data[f"oxide_over_{column}"] - 1.0)
        change_error = 100.0 * data[f"oxide_over_{column}_error"]
        ax.bar(channels, change, yerr=change_error, color=colours, capsize=3, width=0.6)
        ax.axhline(0.0, color="k", linewidth=0.8)
        ax.set_ylabel(label)
        ax.grid(axis="y", alpha=0.3)
        ax.text(0.99, 0.05, "red: PCB-back mounted", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=8, color="tab:red")
    axes[-1, 0].set_xlabel("Channel")
    axes[-1, 0].set_xticks(channels)

    figure.savefig(output, format="pdf", bbox_inches="tight")
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
