"""Drop kernel energies where a fitted channel has exactly zero response.

InferSpectrum.py fits the channel pattern in log space, so it needs every
fitted channel to be strictly positive at every energy it scans. The new
sensitive volumes are about 400 times lighter than the legacy bulk-silicon
block, so a handful of points that the legacy kernel recorded as tiny but
non-zero now come back as exact zeros: no proton out of 8e7 deposited anything
in a 1 mm2 x 400 nm layer.

Those points sit below the channel's range threshold, where the true response is
negligible and the legacy kernel's non-zero value was itself noise-dominated
(64% relative error at the worst of them). A mono-energy hypothesis there is
already excluded by the measured doses, so removing it cannot move the optimum;
it only removes hypotheses the data already reject.

Both kernels are filtered to the SAME energy grid -- the union of the offending
energies across both -- so the oxide and die fits remain directly comparable.
"""

import argparse
from pathlib import Path
import sys

import pandas as pd

# Channels used by either primary fit: VT01 0-6 at 64 MeV, plus channel 10 at
# 85 MeV. Channels 7-9, 11 are diagnostics and are not required to be positive.
FIT_CHANNELS = [0, 1, 2, 3, 4, 5, 6, 10]


def offending_energies(path):
    frame = pd.read_csv(path)
    pivot = frame.pivot(index="energy_MeV", columns="channel", values="mean_dose_rad")
    missing = [c for c in FIT_CHANNELS if c not in pivot.columns]
    if missing:
        raise RuntimeError(f"{path}: kernel lacks fitted channels {missing}")
    bad = pivot.index[(pivot[FIT_CHANNELS] <= 0).any(axis=1)]
    return frame, pivot, set(bad.tolist())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kernel", type=Path, nargs="+", required=True)
    parser.add_argument("--suffix", default="-fitgrid")
    args = parser.parse_args()

    loaded = {}
    drop = set()
    for path in args.kernel:
        frame, pivot, bad = offending_energies(path)
        loaded[path] = (frame, pivot)
        if bad:
            print(f"{path.name}: non-positive fitted channels at {sorted(bad)} MeV")
        drop |= bad

    if not drop:
        print("No energies need dropping; kernels are already usable as they are.")
        return

    print(f"\nCommon grid drops {len(drop)} energies: {sorted(drop)}")
    for path, (frame, pivot) in loaded.items():
        kept = frame[~frame["energy_MeV"].isin(drop)]
        n_energies = kept["energy_MeV"].nunique()
        if n_energies != len(pivot) - len(drop):
            sys.exit(f"{path}: unexpected energy count after filtering")
        if (kept["mean_dose_rad"].isin([0]) & kept["channel"].isin(FIT_CHANNELS)).any():
            sys.exit(f"{path}: a fitted channel is still non-positive after filtering")
        out = path.with_name(path.stem + args.suffix + path.suffix)
        kept.to_csv(out, index=False)
        print(f"  wrote {out.name}: {n_energies} energies, "
              f"{kept['energy_MeV'].min()}-{kept['energy_MeV'].max()} MeV")


if __name__ == "__main__":
    main()
