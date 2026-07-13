#!/usr/bin/env python3
"""Build a flat-sampling GPS energy bias for a SPENVIS / Geant4 GPS spectrum macro.

A SPENVIS-exported source macro defines an arbitrary energy spectrum
(``/gps/ene/type Arb`` + ``/gps/hist/type arb`` points). Space spectra are steep --
most primaries are emitted at low energy -- so a plain Monte-Carlo run barely
samples the high-energy tail, and SPENVIS's own importance bias over-corrects the
other way (it fires most primaries *into* the tail). This tool writes a NEW bias
that makes the Geant4 GPS sampler draw primaries **uniformly in energy** (a flat
number of primaries per energy bin), so the whole spectrum is sampled evenly --
while leaving the physically reconstructed result identical to the unbiased source.
It is purely a variance-reduction aid; it never changes the physics, only where the
Monte-Carlo effort (and thus the statistical precision) is spent.

Works on any GPS arb macro -- electrons, protons, ions (``/gps/ion``), Log- or
Lin-interpolated. It copies the source verbatim and inserts a ``/gps/hist/type
biase`` block; point your simulation at the output instead of the original.

How GPS energy biasing works (Geant4 ``G4SPSEneDistribution`` /
``G4SPSRandomGenerator``):
  * the arb spectrum is sampled by drawing ``u`` in [0,1] and inverting its
    cumulative distribution -> ``E = CDF^-1(u)``; a uniform ``u`` reproduces the
    physical spectrum.
  * the ``biase`` histogram redefines the *density of u*. Its points ``(x_i, y_i)``
    are NOT energies: ``x = u`` is a cumulative probability in [0,1], ``y`` the
    relative sampling weight. GPS builds the cumulative of ``y``; a sample lands in
    u-bin ``i`` with biased probability ``y_i / sum(y)`` and carries weight
    ``w_i = (u_i - u_{i-1}) / (y_i / sum(y))`` = natural / biased probability, so the
    weighted result is unbiased regardless of the bias.

To flatten the sampling we place the bias points at ``u_k = CDF(E_k)`` for ``E_k``
spaced uniformly across the spectrum, with **constant y**: every u-bin then gets
equal sampling probability (one energy bin each -> a flat number of primaries),
while the weight ``w_k`` proportional to ``(u_k - u_{k-1})`` = the true probability
of that bin reconstructs the spectrum. Because ``sum(biased_prob * weight) =
sum(delta_u) = 1``, the mean weight is 1 -- so the SAME normalisation factor
(``NORM_FACTOR_SPECTRUM`` in the SPENVIS macro) keeps the absolute result unchanged.

``--space log`` (default) spaces ``E_k`` logarithmically -> flat on a log-energy
axis; ``--space linear`` spaces them linearly -> flat on a linear-energy axis (and
rising as E on a log one). Only the spacing changes; the weights and normalisation
are identical. For a steep spectrum log is almost always the right choice -- linear
spacing wastes most bins on the negligible high-energy tail and starves the
dominant low-energy bulk.

Standalone: depends only on the Python standard library.

Usage:
    make_flat_bias.py <input_spectrum.mac> <output_biased.mac> [--nbins N] [--space log|linear]

Example (a flat-log and a flat-linear bias for one spectrum):
    make_flat_bias.py AE9.mac AE9-flat-log-biass.mac --space log
    make_flat_bias.py AE9.mac AE9-flat-lin-biass.mac --space linear
"""

from __future__ import annotations

import argparse
import bisect
import math
from pathlib import Path


def parse_spectrum(macro: Path):
    """Read the arb (E, flux) points and interpolation from a GPS spectrum macro.

    Reads only the ``/gps/hist/type arb`` section; anything after a
    ``/gps/hist/type biase`` line is ignored, so re-running on an already-biased
    macro still reads the underlying spectrum.
    """
    points: list[tuple[float, float]] = []
    inter = "Log"
    in_bias = False
    for raw in Path(macro).read_text(encoding="ascii", errors="replace").splitlines():
        s = raw.strip()
        parts = s.split()
        if s.startswith("/gps/hist/type") and len(parts) > 1 and parts[1].lower() == "biase":
            in_bias = True
        elif in_bias:
            continue
        elif s.startswith("/gps/hist/point") and len(parts) == 3:
            points.append((float(parts[1]), float(parts[2])))
        elif s.startswith("/gps/hist/inter"):
            inter = parts[1]
    if len(points) < 2:
        raise ValueError(f"{macro}: no GPS arb spectrum (/gps/hist/point) found")
    points.sort()
    return points, inter


def _segment_integral(ea: float, fa: float, eb: float, fb: float, inter: str = "Log") -> float:
    """Integral of the GPS interpolation of f over [ea, eb] (Lin trapezoid / Log power-law)."""
    if eb <= ea:                                 # zero-width (e at a node) -> 0
        return 0.0
    if inter.lower().startswith("lin") or fa <= 0 or fb <= 0:
        return 0.5 * (fa + fb) * (eb - ea)       # linear interpolation -> trapezoid
    alpha = math.log(fb / fa) / math.log(eb / ea)
    if abs(alpha + 1.0) < 1e-12:
        return fa * ea * math.log(eb / ea)
    return fa / (alpha + 1.0) * (eb * (eb / ea) ** alpha - ea)


def cumulative(points: list[tuple[float, float]], inter: str = "Log"):
    """Return F(E) (cumulative of the spectrum under its interpolation) and its total."""
    es = [e for e, _ in points]
    fs = [f for _, f in points]
    cum = [0.0]
    for i in range(1, len(points)):
        cum.append(cum[-1] + _segment_integral(es[i - 1], fs[i - 1], es[i], fs[i], inter))
    total = cum[-1]
    lin = inter.lower().startswith("lin")

    def f_at(e: float) -> float:
        if e <= es[0]:
            return 0.0
        if e >= es[-1]:
            return total
        j = bisect.bisect_right(es, e) - 1
        if (not lin) and fs[j] > 0 and fs[j + 1] > 0:
            alpha = math.log(fs[j + 1] / fs[j]) / math.log(es[j + 1] / es[j])
            fe = fs[j] * (e / es[j]) ** alpha
        else:
            fe = fs[j] + (fs[j + 1] - fs[j]) * (e - es[j]) / (es[j + 1] - es[j])
        return cum[j] + _segment_integral(es[j], fs[j], e, fe, inter)

    return f_at, total


def build_bias_points(points, nbins: int, space: str = "log",
                      inter: str = "Log") -> list[float]:
    """u_k = CDF(E_k) for nbins E_k spaced uniformly in log- or linear-energy.

    ``space="log"`` flattens the primary count per log-energy bin (flat on a
    log axis); ``space="linear"`` flattens it per linear-energy interval (flat on a
    linear axis, rising as E on a log one). Both keep constant y, so weights
    reconstruct the spectrum identically -- only where the samples land differs. For
    a steep spectrum, linear spacing spends most bins on the negligible-probability
    tail (many collapse to CDF=1) and compresses the dominant low-E bulk into few
    bins, so log is the better default.
    """
    f_at, total = cumulative(points, inter)
    # Span only the spectrum's non-zero-flux support: many SPENVIS spectra pad the
    # high end with zero-flux points (e.g. trapped protons out to 2000 MeV), and
    # bias points placed over a zero-flux tail just collapse to CDF=1, wasting
    # resolution (and degenerating the linear bias). Clip to the first/last
    # positive-flux energy so every bin lands where there are actually particles.
    nz = [e for e, f in points if f > 0]
    if not nz:
        raise ValueError("spectrum has no positive-flux points to bias")
    emin, emax = nz[0], nz[-1]
    if space == "linear":
        es = [emin + (emax - emin) * (k / nbins) for k in range(nbins + 1)]
    else:
        ratio = emax / emin
        es = [emin * ratio ** (k / nbins) for k in range(nbins + 1)]
    us = [f_at(e) / total for e in es]
    # Strictly increasing (drop tail bins that collapse to CDF=1 so no zero-width
    # bin -> infinite weight survives).
    out = [us[0]]
    for u in us[1:]:
        if u > out[-1]:
            out.append(u)
    return out


def write_biased_macro(src: Path, out: Path, nbins: int, space: str = "log") -> int:
    """Copy the source macro verbatim and insert a flat-bias block into it.

    Copying (rather than rebuilding from parsed pieces) keeps everything the source
    needs -- ``/gps/ion`` for ions, comments, exact formatting -- and just adds the
    ``biase`` histogram, exactly where SPENVIS puts it: right after
    ``/gps/hist/inter``. Only the arb (E, f) points + interpolation are read out, to
    compute the bias points' cumulative positions. Returns the number of trailing
    bins dropped because they collapsed to CDF=1.
    """
    src_lines = Path(src).read_text(encoding="ascii", errors="replace").splitlines()
    points, inter = parse_spectrum(src)
    us = build_bias_points(points, nbins, space, inter)
    dropped = (nbins + 1) - len(us)

    bias_block = [
        f"# --- flat-{space}-E importance bias added by make_flat_bias.py ---",
        f"# uniform {space}-energy sampling; weights reconstruct the unbiased result",
        "# (mean weight 1 -> same NORM_FACTOR). Original source kept verbatim above.",
        "/gps/hist/type biase",
        *[f"/gps/hist/point {u:.14f}  1.0" for u in us],
    ]
    # Insert after the (last) /gps/hist/inter line; fall back to after the last arb
    # /gps/hist/point if a macro has no interpolation line.
    last_inter = last_point = None
    for i, line in enumerate(src_lines):
        t = line.strip()
        if t.startswith("/gps/hist/inter"):
            last_inter = i
        elif t.startswith("/gps/hist/point"):
            last_point = i
    anchor = last_inter if last_inter is not None else last_point
    if anchor is None:
        raise ValueError(f"{src}: no /gps/hist spectrum found to attach a bias to")

    result = src_lines[:anchor + 1] + bias_block + src_lines[anchor + 1:]
    out.write_text("\n".join(result) + "\n", encoding="ascii")
    return dropped


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", type=Path, help="unbiased SPENVIS / GPS spectrum macro")
    ap.add_argument("output", type=Path, help="flat-biased macro to write")
    ap.add_argument("--nbins", type=int, default=120,
                    help="number of energy bins to flatten the sampling over")
    ap.add_argument("--space", choices=("log", "linear"), default="log",
                    help="flatten per log-energy bin (default) or per linear-energy interval")
    args = ap.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    dropped = write_biased_macro(args.input, args.output, args.nbins, args.space)
    msg = f"Wrote {args.output} ({args.space}-flat, {args.nbins} bins"
    if dropped:
        msg += f", {dropped} tail bin(s) collapsed to CDF=1 and dropped"
    print(msg + ")")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
