"""Best-fit centre energy with the width fixed to the published relation.

The free-width Gaussian fit is not usable as a beam measurement: at 85 MeV it
pins sigma at the grid floor of 0.05 MeV, and at 64 MeV it wants 1.8 MeV, a
factor of 44 apart and on opposite sides of the commissioning value. Fixing the
width to the publication and scanning only the centre energy removes that free
parameter and asks the well-posed question: given the beam width Skandion
publishes, which centre energy best reproduces the measured channel pattern?

Width relation, from Ardenfors et al., J. Phys.: Conf. Ser. 860 (2017) 012025 as
interpolated by the production spectrum macros: sigma_E/E falls linearly from
0.8% at 60 MeV to 0.3% at 226 MeV, so

    sigma_E(E) = E * (0.8 - 0.5 * (E - 60) / 166) / 100   [MeV]

which reproduces the 0.50429 MeV and 0.61600 MeV in the existing macros exactly.
Note that sigma tracks the fitted centre energy rather than staying at the
nominal value.

The folding, the log-pattern objective and the scale fit are imported from
InferSpectrum, so the numbers stay comparable with every earlier fit.
"""

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

PYTHON_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_ROOT))

from RadExUppsala import InferSpectrum as inf  # noqa: E402

PHASE_REFERENCE = {"64mev": 6.942e10, "85mev": 9.152e10}


def published_sigma(energy_mev):
    """Ardenfors interpolation used by the production spectrum macros."""
    return energy_mev * (0.8 - 0.5 * (energy_mev - 60.0) / 166.0) / 100.0


def check_relation():
    for energy, want in ((64.0, 0.50429), (85.0, 0.61600)):
        got = published_sigma(energy)
        if abs(got - want) > 1e-5:
            raise RuntimeError(
                f"width relation gives {got:.5f} at {energy} MeV, "
                f"macro says {want:.5f}")


def scan(kernel_path, dose_path, phase, channels, fluence_csv):
    cfg = inf.PHASES[phase]
    energies, kernel, _ = inf.load_kernel(kernel_path)

    table = pd.read_csv(dose_path).set_index("channel")
    measured_rad = table.loc[channels, cfg["dose_column"]].to_numpy(float) * 1000.0

    flu = pd.read_csv(fluence_csv).set_index("channel")["delivered_fluence_protons_cm2"]
    weights = (flu.loc[channels] / PHASE_REFERENCE[phase]).to_numpy(float)

    fine_energy = np.arange(float(energies.min()), float(energies.max()) + 0.0125, 0.025)
    fine_kernel = inf.interpolate_kernel(energies, kernel[:, channels], fine_energy)
    fine_kernel = fine_kernel * weights[None, :]

    # Keep the Gaussian well inside the simulated range so truncation cannot
    # bias the result: at least 4 sigma of headroom at both ends.
    lo = float(energies.min()) + 4.0 * published_sigma(float(energies.min()))
    hi = float(energies.max()) - 4.0 * published_sigma(float(energies.max()))
    mus = np.arange(lo, hi + 0.0125, 0.025)

    models = np.empty((len(mus), len(channels)), dtype=float)
    for i, mu in enumerate(mus):
        z = (fine_energy - mu) / published_sigma(mu)
        w = np.exp(-0.5 * z**2)
        w /= w.sum()
        models[i] = w @ fine_kernel

    obj = inf.centered_log_objective(measured_rad, models)
    best = int(np.argmin(obj))
    mu = float(mus[best])
    scale, predicted, rss = inf.fit_scale_log(measured_rad, models[best])
    n = len(channels)
    residual = 100.0 * (predicted / measured_rad - 1.0)

    # 1-sigma interval on mu from the usual delta-chi2 = 1 style criterion on
    # the log-pattern RSS, quoted only as an indication of how sharp the
    # minimum is; the residual is dominated by systematics, not by statistics.
    threshold = rss * (1.0 + 1.0 / max(n - 2, 1))
    inside = mus[obj <= threshold]
    return {
        "mu": mu, "sigma": published_sigma(mu), "rss": rss,
        "rms_percent": 100.0 * float(np.sqrt(rss / n)),
        "scale": scale, "calibration_ratio": scale / PHASE_REFERENCE[phase],
        "residual_percent": residual, "channels": channels,
        "mu_lo": float(inside.min()), "mu_hi": float(inside.max()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kernel", type=Path, required=True)
    parser.add_argument("--dose-64", type=Path, required=True)
    parser.add_argument("--dose-85", type=Path, required=True)
    parser.add_argument("--fluence-64", type=Path, required=True)
    parser.add_argument("--fluence-85", type=Path, required=True)
    args = parser.parse_args()

    check_relation()
    print("Width relation reproduces both production macros exactly.\n")

    cases = [
        ("64 MeV, ch 0-6", "64mev", [0, 1, 2, 3, 4, 5, 6], args.dose_64, args.fluence_64),
        ("85 MeV, ch 0-6 + 10", "85mev", [0, 1, 2, 3, 4, 5, 6, 10], args.dose_85, args.fluence_85),
        ("85 MeV, ch 0-6 (no Ch10)", "85mev", [0, 1, 2, 3, 4, 5, 6], args.dose_85, args.fluence_85),
    ]

    print(f"{'case':<26s}{'mu [MeV]':>10s}{'sigma':>8s}{'1-sigma range':>18s}"
          f"{'RMS':>8s}{'C':>7s}")
    results = {}
    for label, phase, channels, dose, fluence in cases:
        r = scan(args.kernel, dose, phase, channels, fluence)
        results[label] = r
        span = f"{r['mu_lo']:.2f} - {r['mu_hi']:.2f}"
        print(f"{label:<26s}{r['mu']:>10.2f}{r['sigma']:>8.3f}{span:>18s}"
              f"{r['rms_percent']:>7.2f}%{r['calibration_ratio']:>7.3f}")

    print("\nPer-channel residuals [%]")
    for label, r in results.items():
        pairs = "  ".join(f"ch{c}:{v:+.1f}" for c, v in
                          zip(r["channels"], r["residual_percent"]))
        print(f"  {label:<26s}{pairs}")


if __name__ == "__main__":
    main()
