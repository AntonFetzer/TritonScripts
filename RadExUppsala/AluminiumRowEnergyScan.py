"""Constrain a phase's centre energy with same-type channel pairs in the aluminium row.

Ch7-9 are VT05 parts whose dose conversion rests on a single X-ray calibration
fit, so their absolute doses are not directly comparable with the VT01 channels.
Ratios WITHIN one device type need much less from the calibration:

    R65 = D6 / D5   VT01 under the 2 mm absorber over VT01 behind the open aperture
    R78 = D7 / D8   VT05 under the 4 mm over VT05 under the 5 mm absorber

Near the end of the proton range both ratios change steeply with energy, so
they test the beam energy and the absorber thicknesses independently of the
VT01 dose pattern.

For each centre energy the kernel is folded with the published Gaussian width
and each channel is weighted by its delivered fluence, through
InferSpectrum.weighted_kernel and InferSpectrum.fold. The measured R78 is quoted
as a band: the tabulated doses use the datasheet conversion
dose = (shift/a)^(1/b) with b = InferSpectrum.VT05_EXPONENT, while a response
linear in dose would give the shift ratio, i.e. the tabulated ratio raised to
the power b. R65 is compared within the VT01 device scatter,
InferSpectrum.SYSTEMATIC.

Run from the repository root:

    python3 -m RadExUppsala.AluminiumRowEnergyScan --kernel <energy-kernel-dose-summary-oxide.csv> \\
        --phase 64mev --dose <phase-end-dose-64mev.csv> \\
        --fluence <RadEx_64MeV/radfet_location_fluence.csv> --output <stem>

writes <stem>.csv and <stem>.pdf.
"""

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PYTHON_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_ROOT))

from RadExUppsala import InferSpectrum as inf  # noqa: E402

ROW_CHANNELS = [5, 6, 7, 8, 9]
#: Spacing of the printed ratio table [MeV].
PRINT_STEP = 0.5


def ratio(means, sigmas, num, den):
    a, b = means[:, ROW_CHANNELS.index(num)], means[:, ROW_CHANNELS.index(den)]
    sa, sb = sigmas[:, ROW_CHANNELS.index(num)], sigmas[:, ROW_CHANNELS.index(den)]
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(b > 0, a / b, np.nan)
        sr = np.abs(r) * np.sqrt((sa / a) ** 2 + (sb / b) ** 2)
    return r, sr


def intervals(mus, inside):
    """Contiguous [lo, hi] centre-energy ranges where `inside` is true."""
    out, start, prev = [], None, None
    for mu, ok in zip(mus, inside):
        if ok and start is None:
            start = mu
        if not ok and start is not None:
            out.append((start, prev))
            start = None
        prev = mu
    if start is not None:
        out.append((start, prev))
    return out


def fmt(ranges):
    return ", ".join(f"{lo:.2f}-{hi:.2f}" for lo, hi in ranges) or "none"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--kernel", type=Path, required=True,
                        help="energy-kernel-dose-summary CSV from BuildKernelSummary")
    parser.add_argument("--phase", choices=sorted(inf.PHASES), default="64mev")
    parser.add_argument("--dose", type=Path, required=True)
    parser.add_argument("--fluence", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="output stem")
    args = parser.parse_args()

    phase = args.phase
    d = pd.read_csv(args.dose).set_index("channel")[inf.PHASES[phase]["dose_column"]]
    r65_meas = float(d[6] / d[5])
    r78_datasheet = float(d[7] / d[8])
    r78_linear = r78_datasheet ** inf.VT05_EXPONENT
    r78_lo, r78_hi = sorted((r78_datasheet, r78_linear))
    tolerance = inf.SYSTEMATIC

    print(f"Measured {phase} phase: R65 = D6/D5 = {r65_meas:.3f} (VT01, +-{100 * tolerance:.0f}% "
          f"device scatter assumed)")
    print(f"                       R78 = D7/D8 = {r78_datasheet:.2f} with the VT05 datasheet exponent, "
          f"{r78_linear:.2f} if linear in dose")
    print(f"                       D9 = {d[9]:.4f} kRad, D8 = {d[8]:.3f} kRad, D7 = {d[7]:.2f} kRad\n")

    lo, hi = inf.PHASES[phase]["mu_range"]
    mus = np.arange(lo, hi + inf.FINE_STEP / 2, inf.FINE_STEP)
    energies, kernel, error = inf.weighted_kernel(args.kernel, args.fluence, phase, ROW_CHANNELS)
    means, sigmas = inf.fold(energies, kernel, error, mus)
    r65, s65 = ratio(means, sigmas, 6, 5)
    r78, s78 = ratio(means, sigmas, 7, 8)
    ok65 = np.abs(r65 / r65_meas - 1.0) <= np.hypot(tolerance, s65 / r65)
    ok78 = (r78 + s78 >= r78_lo) & (r78 - s78 <= r78_hi)
    print(f"R65 consistent for mu in {fmt(intervals(mus, ok65))} MeV")
    print(f"R78 consistent for mu in {fmt(intervals(mus, ok78))} MeV")
    print(f"both consistent for mu in {fmt(intervals(mus, ok65 & ok78))} MeV")
    for mu_show in np.arange(lo, hi + PRINT_STEP / 2, PRINT_STEP):
        i = int(np.argmin(np.abs(mus - mu_show)))
        print(f"    mu {mus[i]:6.2f}: R65 {r65[i]:.3f}+-{s65[i]:.3f}   R78 {r78[i]:8.2f}+-{s78[i]:.2f}")

    rows = [{"mu_MeV": round(float(mu), 4), "sigma_MeV": inf.published_sigma(float(mu)),
             "R65": r65[i], "R65_error": s65[i], "R78": r78[i], "R78_error": s78[i],
             **{f"dose_ch{c}_rad": means[i, k] for k, c in enumerate(ROW_CHANNELS)}}
            for i, mu in enumerate(mus)]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, r, s in ((axes[0], r65, s65), (axes[1], r78, s78)):
        line, = ax.plot(mus, r, label="model")
        ax.fill_between(mus, r - s, r + s, color=line.get_color(), alpha=0.25, lw=0)
    axes[0].axhspan(r65_meas * (1 - tolerance), r65_meas * (1 + tolerance),
                    color="grey", alpha=0.3, label=f"measured +-{100 * tolerance:.0f}%")
    axes[0].axhline(r65_meas, color="black", lw=1)
    axes[0].set_ylabel("D6 / D5 (VT01)")
    axes[1].axhspan(r78_lo, r78_hi, color="grey", alpha=0.3, label="measured band")
    axes[1].set_yscale("log")
    axes[1].set_ylabel("D7 / D8 (VT05)")
    for ax in axes:
        ax.set_xlabel(f"{phase.replace('mev', ' MeV')} phase centre energy [MeV]")
        ax.grid(alpha=0.3)
        ax.legend()
    fig.tight_layout()

    output = args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output.with_suffix(".csv"), index=False)
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    print(f"Wrote {output.with_suffix('.csv')} and {output.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
