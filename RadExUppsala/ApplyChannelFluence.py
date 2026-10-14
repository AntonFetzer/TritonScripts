"""Weight a response kernel by the per-channel delivered fluence.

InferSpectrum.py fits

    D_ch = scale * integral(K_ch(E) p(E) dE)

with one free `scale` shared by every channel, which assumes every RadFET saw
the same fluence. The Uppsala accelerator logs say otherwise: bilinear sampling
of the ionisation-chamber heatmap at the twelve RadFET locations gives fluences
spanning about 8% at 85 MeV and 11% at 64 MeV, and all twelve sit below the
instrument-area mean that the fit currently normalises to.

The correct model is

    D_ch = C * Phi_ch * integral(K_ch(E) p(E) dE)

with C a single calibration constant. Folding Phi_ch into the kernel,

    K'_ch(E) = K_ch(E) * Phi_ch / Phi_ref

recovers a model with one shared scale again, so InferSpectrum.py can be used
unmodified. The fitted scale is then C * Phi_ref, and with Phi_ref set to the
phase's log-derived fluence the reported `fluence_ratio_to_log` is exactly C.

This scales the SIMULATED kernel, never the measured doses. The delivered
fluences differ between the 64 and 85 MeV phases, so one weighted kernel serves
one phase only; build both and read each phase from its matching fit.
"""

import argparse
import csv
from pathlib import Path
import sys

import pandas as pd

# Reference fluences are the log-derived, instrument-area means already used by
# InferSpectrum.py, so the ratio_to_log column keeps its meaning.
PHASE_REFERENCE = {"64mev": 6.942e10, "85mev": 9.152e10}


def load_channel_fluence(path):
    frame = pd.read_csv(path)
    required = {"channel", "delivered_fluence_protons_cm2"}
    missing = required.difference(frame.columns)
    if missing:
        raise RuntimeError(f"{path}: missing columns {sorted(missing)}")
    series = frame.set_index("channel")["delivered_fluence_protons_cm2"].astype(float)
    if series.le(0).any():
        raise RuntimeError(f"{path}: non-positive delivered fluence")
    return series


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kernel", type=Path, required=True)
    parser.add_argument("--fluence", type=Path, required=True)
    parser.add_argument("--phase", choices=sorted(PHASE_REFERENCE), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    reference = PHASE_REFERENCE[args.phase]
    fluence = load_channel_fluence(args.fluence)
    kernel = pd.read_csv(args.kernel)

    channels = sorted(kernel["channel"].unique())
    missing = [c for c in channels if c not in fluence.index]
    if missing:
        sys.exit(f"No delivered fluence for kernel channels {missing}")

    factor = (fluence / reference).reindex(channels)
    weighted = kernel.copy()
    weights = weighted["channel"].map(factor)
    for column in ("mean_dose_rad", "combined_error_rad", "empirical_sem_rad"):
        if column in weighted.columns:
            weighted[column] = weighted[column] * weights
    # relative_error is a ratio, so the common factor cancels; leave it alone.
    weighted.to_csv(args.output, index=False)

    print(f"Phase {args.phase}, reference fluence {reference:.4g} protons/cm2")
    print(f"{'ch':>3s} {'delivered':>14s} {'factor':>8s}")
    for channel in channels:
        print(f"{channel:>3d} {fluence[channel]:>14.4e} {factor[channel]:>8.4f}")
    print(f"\nspread over channels: "
          f"{100 * (fluence[channels].max() / fluence[channels].min() - 1):.2f} %")
    print(f"mean factor vs reference: {factor.mean():.4f}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
