"""Constrain a phase's centre energy with same-type channel dose ratios.

Absolute doses of different device types rest on different calibrations, so
they are not directly comparable. Ratios WITHIN one device type need much less
from the calibration. At Uppsala the aluminium row gives two such pairs:

    6/5   VT01 under the 2 mm absorber over VT01 behind the open aperture
    7/8   VT05 under the 4 mm over VT05 under the 5 mm absorber

Near the end of the proton range these ratios change steeply with energy, so
they test the beam energy and the absorber thicknesses independently of the
VT01 dose pattern. Each --pair NUM/DEN names one ratio; both channels must have
the same device type in the campaign's dose table.

For each centre energy the kernel is folded with the published Gaussian width
and each channel is weighted by its delivered fluence, through
InferSpectrum.weighted_kernel and InferSpectrum.fold. A VT05 ratio is quoted as
a band: the tabulated doses use the datasheet conversion
dose = (shift/a)^(1/b) with b = InferSpectrum.VT05_EXPONENT, while a response
linear in dose would give the shift ratio, i.e. the tabulated ratio raised to
the power b. A VT01 ratio is compared within the campaign's device scatter,
`systematic_fraction`.

Run from the repository root:

    python3 -m RadEx.AluminiumRowEnergyScan --kernel <energy-kernel-dose-summary-oxide.csv> \\
        --campaign <campaign.json> --phase 64mev --pair 6/5 --pair 7/8 \\
        --channels 5 6 7 8 9 --output <stem>

writes <stem>.csv (ratios and the folded dose of every --channels channel) and
<stem>.pdf.
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

from RadEx import InferSpectrum as inf  # noqa: E402

#: Spacing of the printed ratio table [MeV].
PRINT_STEP = 0.5


def parse_pair(text):
    try:
        num, den = (int(v) for v in text.split("/"))
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected NUM/DEN channels, got {text!r}") from None
    return num, den


def ratio(means, sigmas, channels, num, den):
    a, b = means[:, channels.index(num)], means[:, channels.index(den)]
    sa, sb = sigmas[:, channels.index(num)], sigmas[:, channels.index(den)]
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
    parser.add_argument("--campaign", type=Path, required=True,
                        help="campaign JSON next to the measured doses")
    parser.add_argument("--phase", required=True)
    parser.add_argument("--pair", type=parse_pair, action="append", required=True,
                        help="NUM/DEN same-type channel ratio, e.g. 6/5; repeatable")
    parser.add_argument("--channels", type=int, nargs="+",
                        help="channels whose folded doses are tabulated (default: the pair channels)")
    parser.add_argument("--output", type=Path, required=True, help="output stem")
    args = parser.parse_args()

    campaign = inf.load_campaign(args.campaign)
    phase = args.phase
    if phase not in campaign["phases"]:
        parser.error(f"unknown phase {phase!r}; the campaign has {list(campaign['phases'])}")
    device = campaign["device"]
    pair_channels = sorted({c for pair in args.pair for c in pair})
    channels = args.channels or pair_channels
    missing = [c for c in pair_channels if c not in channels] + [c for c in channels if c not in device]
    if missing:
        parser.error(f"channels {missing} are not in --channels or have no measured dose")
    for num, den in args.pair:
        if device[num] != device[den]:
            parser.error(f"pair {num}/{den} mixes {device[num]} and {device[den]}")

    d = campaign["phases"][phase]["table"][campaign["phases"][phase]["dose_column"]]
    tolerance = campaign["systematic"]
    pairs = []
    for num, den in args.pair:
        measured = float(d[num] / d[den])
        name = f"R{num}{den}"
        if device[num] == "VT05":
            linear = measured ** inf.VT05_EXPONENT
            band = tuple(sorted((measured, linear)))
            print(f"Measured {phase} phase: {name} = D{num}/D{den} = {measured:.2f} with the VT05 "
                  f"datasheet exponent, {linear:.2f} if linear in dose")
        else:
            band = (measured * (1 - tolerance), measured * (1 + tolerance))
            print(f"Measured {phase} phase: {name} = D{num}/D{den} = {measured:.3f} ({device[num]}, "
                  f"+-{100 * tolerance:.0f}% device scatter assumed)")
        pairs.append({"num": num, "den": den, "name": name, "measured": measured, "band": band,
                      "device": device[num]})
    print("                       " + ", ".join(f"D{c} = {d[c]:.4g} kRad" for c in reversed(channels)) + "\n")

    lo, hi = campaign["phases"][phase]["mu_range"]
    mus = np.arange(lo, hi + inf.FINE_STEP / 2, inf.FINE_STEP)
    energies, kernel, error = inf.weighted_kernel(args.kernel, campaign, phase, channels)
    means, sigmas = inf.fold(energies, kernel, error, mus)
    every = np.ones(len(mus), dtype=bool)
    for p in pairs:
        r, s = ratio(means, sigmas, channels, p["num"], p["den"])
        if p["device"] == "VT05":
            ok = (r + s >= p["band"][0]) & (r - s <= p["band"][1])
        else:
            ok = np.abs(r / p["measured"] - 1.0) <= np.hypot(tolerance, s / r)
        p.update(model=r, error=s, ok=ok)
        every &= ok
        print(f"{p['name']} consistent for mu in {fmt(intervals(mus, ok))} MeV")
    if len(pairs) > 1:
        print(f"all pairs consistent for mu in {fmt(intervals(mus, every))} MeV")
    for mu_show in np.arange(lo, hi + PRINT_STEP / 2, PRINT_STEP):
        i = int(np.argmin(np.abs(mus - mu_show)))
        print(f"    mu {mus[i]:6.2f}: " + "   ".join(
            f"{p['name']} {p['model'][i]:.3f}+-{p['error'][i]:.3f}" for p in pairs))

    rows = []
    for i, mu in enumerate(mus):
        row = {"mu_MeV": round(float(mu), 4), "sigma_MeV": inf.published_sigma(float(mu))}
        for p in pairs:
            row[p["name"]] = p["model"][i]
            row[f"{p['name']}_error"] = p["error"][i]
        row.update({f"dose_ch{c}_rad": means[i, k] for k, c in enumerate(channels)})
        rows.append(row)

    fig, axes = plt.subplots(1, len(pairs), figsize=(5.5 * len(pairs), 4.2), squeeze=False)
    for ax, p in zip(axes[0], pairs):
        line, = ax.plot(mus, p["model"], label="model")
        ax.fill_between(mus, p["model"] - p["error"], p["model"] + p["error"],
                        color=line.get_color(), alpha=0.25, lw=0)
        if p["device"] == "VT05":
            ax.axhspan(*p["band"], color="grey", alpha=0.3, label="measured band")
            ax.set_yscale("log")
        else:
            ax.axhspan(*p["band"], color="grey", alpha=0.3, label=f"measured +-{100 * tolerance:.0f}%")
            ax.axhline(p["measured"], color="black", lw=1)
        ax.set_ylabel(f"D{p['num']} / D{p['den']} ({p['device']})")
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
