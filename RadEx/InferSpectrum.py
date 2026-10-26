"""Infer effective proton energies from RadFET phase-end measurements and a response kernel.

This module also holds the measurement model that FitAllChannels and
AluminiumRowEnergyScan share, so that every fit uses the same constants,
uncertainties and kernel folding.

The campaign is described by a JSON file kept next to the measured doses (for
Uppsala, GRAS/RadEx/RadEx-Uppsala/MeasuredDoses/campaign.json), so the same
code fits any campaign of these devices. Per phase it names the measured
phase-end table, the dose column, the phase it increments on (if any), the
accelerator-log per-channel fluence, the log fluence, the nominal energy, the
searched centre-energy range and the calibration-scale models (each of which may
be given per fit space, as {"dose": ..., "shift": ...}); for the campaign
it gives the fit space, the systematic device scatter, the response exponents
and the channels that faded between phases. Relative paths resolve against the
campaign file. The channels and their device types come from the first phase's
table (columns `channel` and `group`).

Fit spaces (`fit_space`, default "dose"; --space overrides it):

    dose    Every channel enters in linear dose space: the measured dose is the
            phase-end shift inverted through the vendor calibration, and the
            model is scale * K, weighted by

                sigma^2 = (systematic_fraction * D)^2 + sigma_noise^2 + (scale * sigma_K)^2

            where sigma_noise is the endpoint-fit RMS voltage `fit_rms_mv`
            turned into dose with the local slope of the datasheet calibration,
            D_tot / (b * V_tot) at the phase end (for a phase that increments on
            a previous one, the noise of both endpoints adds in quadrature). A
            channel whose endpoint shift is below its noise enters with zero
            dose and the noise converted through the lowest-range fit.

    shift   Every channel enters as its phase-end threshold-voltage SHIFT, and
            the model is

                dV = A_g * ((D + D0)^n - D0^n),    D = D_prev + K * F_log

            with one amplitude A per scale group, the response exponent n and
            the knee dose D0 per device type from the campaign's
            `shift_exponent` and `shift_offset_krad`, and D_prev the previous
            phase's unscaled model dose at its fitted energy (so a phase on
            pre-dosed devices is fitted as a total shift, with no dose-space
            subtraction and no inverted calibration). The response is linear
            below D0 and the power law above it; D0 = 0 is the pure power law.
            Defaults: VT01 n = 0.51, the 64 MeV phase's own fitted exponent on
            the sensor board, and D0 = 0.36 krad, which puts the vendor's 1 and
            10 krad calibration points in their published ratio at that
            exponent; VT05 n = 0.64 (the datasheet exponent, all the single
            VT05 channel with 64 MeV dose allows) and D0 = 0. The weights are

                sigma^2 = noise^2 + (d(dV)/dD * D)^2 * (systematic_fraction^2 + (sigma_K / D)^2)

            with noise = `fit_rms_mv` of the phase's own endpoint. Nothing of
            the vendor calibration enters the fit. Amplitudes are always split
            by device type in this space, because the two types' responses
            have different exponents and amplitudes: `shared` becomes
            `per-type` and `fading` becomes `fading-per-type`.

            `shift_response_table` (optional, type -> CSV with dose_krad and
            shift_v, relative to the campaign file) replaces the power law of
            that device type by the tabulated curve, log-log interpolated
            (Dependencies.EnergyKernel.tabulatedResponse). With the vendor
            calibration as the table the amplitude is the ratio of the
            devices to the calibration, and the response carries the curve's
            own falling slope instead of one exponent.

Calibration scales are solved by iteratively reweighted least squares. Each
phase's own per-channel delivered fluence multiplies the kernel.

Calibration-scale models (each phase's `scale_model` is its default):

    shared           one scale for all channels.
    fading           two scales: the campaign's `faded_channels` and all others.
                     At Uppsala these are the pre-dosed VT01, which lost 31-48 mV
                     in the 15 min beam-off pause before the 85 MeV phase,
                     against no drift for VT05 Ch7 (pre-dosed, 208 mV) and the
                     fresh channels. The grouping comes from the pause, not from
                     the fit.
    per-type         one scale per device type (VT01: Co-60, rad(Si); VT05:
                     X-ray, rad(air)).
    fading-per-type  faded and unfaded channels split again by device type: at
                     Uppsala the pre-dosed VT01 (Ch0-6), the fresh VT01 (Ch10)
                     and the VT05 (Ch7-9, 11). Ch10 reproduces the 64 MeV
                     fresh-VT01 scale (0.29 against 0.31), so it belongs with
                     neither the faded VT01 nor the VT05.

Three beam models are compared:

1. a monoenergetic beam;
2. a Gaussian energy distribution over (mu, sigma), truncated to the kernel range;
3. a Gaussian with the published width at each centre energy (as in FitAllChannels).

Optionally, beside them (off by default, so the published width stays the
reference):

4. --fixed-width-mev W: a Gaussian of fixed width W at every centre energy (at
   Skandion the facility's quoted energy-spread ceiling, 0.3 MeV), written under
   `fixed_width` in the summary and as `fixed_width_*` columns;
5. --fixed-width-cut-sigma C (with 4): the same Gaussian cut off above
   mu + C sigma, the facility's description of a spectrum targeted at the
   Bragg-peak position, under `fixed_width_truncated`. Its mu is the Gaussian's
   centre; `spectrum_mean_energy_MeV` is the mean of the cut spectrum.

Like the published width, each fixed-width model is a mu scan of the phase's
default scale model; in shift space a phase on pre-dosed devices starts from the
previous phase's model dose of the SAME fixed-width model, so each model is a
self-consistent beam hypothesis over the run.

The summary reports the default model's fits and, under "diagnostics", those of
the phase's diagnostic models; the figure shows the default fits and marks the
diagnostic Gaussians.

Intervals come from delta-chi2 <= 1, profiled over the other parameter. The
searched centre energies default to each phase's `mu_range_MeV` (at Uppsala the
nominal setting +-3 MeV, 61-67 and 82-88 MeV), which also sets the energy axis
of the misfit maps, and can be changed with --mu-range PHASE LO HI.

The folding and the scale and amplitude solvers are Dependencies.EnergyKernel; this
module holds the campaign file, the measurement models and the fits.

Run from the repository root:

    python3 -m RadEx.InferSpectrum --kernel <energy-kernel-dose-summary-oxide.csv> \\
        --campaign <campaign.json> --output <stem> [--space dose|shift] \\
        [--fixed-width-mev 0.3 [--fixed-width-cut-sigma 1.0]]

writes <stem>-summary.json, <stem>-channel-residuals.csv and <stem>-diagnostics.pdf.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PYTHON_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_ROOT))

from Dependencies import EnergyKernel as ek  # noqa: E402


MARKER = {"VT01": "o", "VT05": "s"}
DIAGNOSTIC_COLOURS = ("grey", "magenta", "green")
#: Colours of the optional fixed-width models (full Gaussian, cut Gaussian).
FIXED_WIDTH_COLOURS = ("darkorange", "saddlebrown")
#: Calibration-scale models; each phase names its default and diagnostics in the campaign.
SCALE_MODELS = ("shared", "fading", "per-type", "fading-per-type")
FIT_SPACES = ("dose", "shift")
#: Exponents B of the Varadis datasheet fits dV = A * dose_rad**B (QF 16 test
#: records, 13/05/2025, as in fs2_radex ProtonAnalysis.DATASHEET_FITS), by the
#: upper end in krad of the dose range each fit covers.
VT01_EXPONENT_BY_RANGE = ((1.0, 0.91193), (5.0, 0.70621), (10.0, 0.61638),
                          (50.0, 0.47515), (100.0, 0.45509))
VT05_EXPONENT = 0.64146
#: Lowest-range datasheet fit (A [V], B) per type, used to turn a noise voltage
#: into a dose for a channel whose shift is below its noise.
NOISE_FIT = {"VT01": (7.83175e-4, 0.91193), "VT05": (4.35124e-4, 0.64146)}
#: Default response exponents for the shift space when the campaign gives none:
#: VT01 from the 64 MeV phase's own power-law fits (mean of Ch0-6), VT05 the
#: datasheet exponent, which is all the one VT05 channel with 64 MeV dose allows.
DEFAULT_SHIFT_EXPONENT = {"VT01": 0.51, "VT05": VT05_EXPONENT}
#: Default knee doses D0 [krad] of the shift-space response, below which it is
#: linear: VT01 from the vendor's 1 krad / 10 krad shift ratio (0.4262 / 1.9979)
#: at n = 0.51; VT05 unknown, pure power law.
DEFAULT_SHIFT_OFFSET_KRAD = {"VT01": 0.36, "VT05": 0.0}
#: Step of the fine energy grid for interpolation, folding and the mu search [MeV].
FINE_STEP = 0.025
SIGMAS = np.arange(0.05, 4.0001, 0.05)


def published_sigma(energy_mev):
    """Skandion energy spread [MeV] at a centre energy [MeV] (Dependencies.EnergyKernel)."""
    return ek.publishedSkandionSigma(energy_mev)


def fixed_width_specs(width_mev: float | None, cut_sigma: float | None = None) -> list[dict]:
    """The optional fixed-width beam models: none, the Gaussian, and the Gaussian cut above mu + cut sigma.

    Each spec has 'key' (summary key and column prefix), 'sigma' [MeV], 'cut'
    (None or sigmas above the centre), 'label' and 'colour'.
    """
    if width_mev is None:
        if cut_sigma is not None:
            raise SystemExit("--fixed-width-cut-sigma needs --fixed-width-mev")
        return []
    if width_mev <= 0 or (cut_sigma is not None and cut_sigma <= 0):
        raise SystemExit("--fixed-width-mev and --fixed-width-cut-sigma must be positive")
    specs = [{"key": "fixed_width", "sigma": float(width_mev), "cut": None,
              "label": f"{width_mev:g} MeV width", "colour": FIXED_WIDTH_COLOURS[0]}]
    if cut_sigma is not None:
        specs.append({"key": "fixed_width_truncated", "sigma": float(width_mev), "cut": float(cut_sigma),
                      "label": f"{width_mev:g} MeV width cut at +{cut_sigma:g} sigma",
                      "colour": FIXED_WIDTH_COLOURS[1]})
    return specs


def spectrum_mean(mu: float, sigma: float, cut: float | None) -> float:
    """Mean energy [MeV] of a Gaussian (mu, sigma), cut off above mu + cut * sigma if cut is given."""
    if cut is None:
        return float(mu)
    phi = np.exp(-0.5 * cut ** 2) / np.sqrt(2.0 * np.pi)
    big_phi = 0.5 * (1.0 + math.erf(cut / np.sqrt(2.0)))
    return float(mu - sigma * phi / big_phi)


def add_width_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--fixed-width-mev", type=float, default=None,
                        help="also fit a Gaussian of this fixed energy spread [MeV] beside the published "
                             "width (e.g. 0.3, the facility's ceiling); off by default")
    parser.add_argument("--fixed-width-cut-sigma", type=float, default=None,
                        help="with --fixed-width-mev, also fit that Gaussian cut off above "
                             "mu + this many sigma (e.g. 1.0)")


def _by_space(setting, space: str):
    """A campaign setting that may be given per fit space, as {"dose": ..., "shift": ...}."""
    return setting[space] if isinstance(setting, dict) else setting


def load_campaign(path: Path, space: str | None = None, systematic: float | None = None,
                  noise_floor_mv: float | None = None) -> dict:
    """Read a campaign file; see the module docstring for its fields.

    Returns a dict with 'path', 'channels' (sorted), 'device' (channel -> type),
    'space', 'systematic', 'noise_floor' [V], 'exponent' (type -> n),
    'faded_channels' and 'phases' (label -> settings, with the table loaded
    under 'table'), phases in file order. `space`, `systematic` and
    `noise_floor_mv` override the file's `fit_space`, `systematic_fraction`
    and `noise_floor_mv`.

    The noise floor is added in quadrature to every channel's endpoint RMS:
    the baseline drift over a phase, which a channel without dose shows
    directly (Ch10 at Uppsala: 4.3 mV at 64 MeV against 1.0 mV RMS) and which
    the RMS of a 300 s window cannot see.
    """
    path = Path(path)
    data = json.loads(path.read_text(encoding="utf-8"))
    space = space or data.get("fit_space", "dose")
    if space not in FIT_SPACES:
        raise ValueError(f"{path}: fit_space {space!r}; use one of {FIT_SPACES}")
    phases = {}
    for label, cfg in data["phases"].items():
        previous = cfg.get("previous_phase")
        if previous is not None and previous not in phases:
            raise ValueError(f"{path}: {label} increments on {previous!r}, which must be listed before it")
        models = (_by_space(cfg["scale_model"], space), *_by_space(cfg.get("diagnostic_models", ()), space))
        unknown = [m for m in models if m not in SCALE_MODELS]
        if unknown:
            raise ValueError(f"{path}: {label} has unknown scale models {unknown}; use {SCALE_MODELS}")
        dose_path = path.parent / cfg["dose"]
        phases[label] = {
            "dose": dose_path,
            "table": pd.read_csv(dose_path).set_index("channel"),
            "dose_column": cfg["dose_column"],
            "previous_phase": previous,
            "fluence": path.parent / cfg["fluence"],
            "log_fluence": float(cfg["log_fluence_protons_cm2"]),
            "nominal_energy": float(cfg["nominal_energy_MeV"]),
            "mu_range": tuple(float(v) for v in cfg["mu_range_MeV"]),
            "scale_model": models[0],
            "diagnostic_models": models[1:],
        }
    if not phases:
        raise ValueError(f"{path}: no phases")
    first = next(iter(phases.values()))["table"]
    channels = sorted(int(c) for c in first.index)
    device = {c: str(first.loc[c, "group"]) for c in channels}
    for label, cfg in phases.items():
        table = cfg["table"]
        if sorted(int(c) for c in table.index) != channels or any(
                str(table.loc[c, "group"]) != device[c] for c in channels):
            raise ValueError(f"{cfg['dose']}: channels or device types differ from the first phase")
    unknown = sorted({d for d in device.values() if d not in MARKER})
    if unknown:
        raise ValueError(f"{path}: no datasheet calibration for device types {unknown}")
    faded = tuple(int(c) for c in data.get("faded_channels", ()))
    if any(c not in channels for c in faded):
        raise ValueError(f"{path}: faded_channels lists channels that have no measured dose")
    exponent = dict(DEFAULT_SHIFT_EXPONENT)
    exponent.update({k: float(v) for k, v in data.get("shift_exponent", {}).items()})
    offset = dict(DEFAULT_SHIFT_OFFSET_KRAD)
    offset.update({k: float(v) for k, v in data.get("shift_offset_krad", {}).items()})
    floor = float(data.get("noise_floor_mv", 0.0)) if noise_floor_mv is None else float(noise_floor_mv)
    tables = {}
    for kind, name in data.get("shift_response_table", {}).items():
        frame = pd.read_csv(path.parent / name)
        if not {"dose_krad", "shift_v"}.issubset(frame.columns):
            raise ValueError(f"{name}: a response table needs dose_krad and shift_v")
        tables[kind] = (frame["dose_krad"].to_numpy(float) * 1000.0, frame["shift_v"].to_numpy(float))
    return {
        "path": path, "channels": channels, "device": device, "space": space,
        "systematic": float(data["systematic_fraction"]) if systematic is None else float(systematic),
        "noise_floor": floor / 1000.0,
        "exponent": exponent, "offset_krad": offset, "response_table": tables,
        "faded_channels": faded, "phases": phases,
    }


def load_kernel(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    frame = pd.read_csv(path)
    required = {"energy_MeV", "channel", "mean_dose_rad", "combined_error_rad"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Kernel is missing columns: {sorted(missing)}")
    dose = frame.pivot(index="energy_MeV", columns="channel", values="mean_dose_rad")
    err = frame.pivot(index="energy_MeV", columns="channel", values="combined_error_rad")
    dose = dose.sort_index()
    err = err.reindex(index=dose.index, columns=dose.columns)
    if dose.isna().any().any():
        raise ValueError("Kernel is not a complete energy-by-channel matrix")
    return dose.index.to_numpy(float), dose.to_numpy(float), err.to_numpy(float)


def channel_weights(campaign: dict, phase: str, channels: list[int]) -> np.ndarray:
    """Delivered fluence per channel relative to the phase's log fluence."""
    cfg = campaign["phases"][phase]
    frame = pd.read_csv(cfg["fluence"]).set_index("channel")
    fluence = frame.loc[channels, "delivered_fluence_protons_cm2"].to_numpy(float)
    if np.any(fluence <= 0):
        raise ValueError(f"{cfg['fluence']}: non-positive delivered fluence")
    return fluence / cfg["log_fluence"]


def local_exponent(device: str, dose_krad: float) -> float:
    if device == "VT05":
        return VT05_EXPONENT
    for top, exponent in VT01_EXPONENT_BY_RANGE:
        if dose_krad <= top:
            return exponent
    return VT01_EXPONENT_BY_RANGE[-1][1]


def measurements(campaign: dict, phase: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Dose space: measured dose [rad], its 1-sigma uncertainty [rad] and the noise part.

    `dose_krad` is the total dose at the phase end, which sets the calibration
    slope; the phase's `dose_column` is the dose fitted. A phase that increments
    on a previous one adds that phase's endpoint noise in quadrature. A channel
    whose endpoint shift is below its noise (Ch9-11 at 64 MeV) has no inverse
    through the calibration slope; it enters with its dose clipped at zero and
    the noise converted through the lowest-range datasheet fit.
    """
    channels = campaign["channels"]
    cfg = campaign["phases"][phase]
    table = cfg["table"].loc[channels]
    dose = np.nan_to_num(table[cfg["dose_column"]].to_numpy(float) * 1000.0)
    total = np.nan_to_num(table["dose_krad"].to_numpy(float) * 1000.0)
    v_total = table["shift_end_v"].to_numpy(float)
    rms_mv = table["fit_rms_mv"].to_numpy(float)
    if cfg["previous_phase"] is not None:
        previous = campaign["phases"][cfg["previous_phase"]]["table"].loc[channels]
        rms_mv = np.hypot(rms_mv, previous["fit_rms_mv"].to_numpy(float))
    noise_v = np.hypot(rms_mv / 1000.0, campaign["noise_floor"])
    noise = np.empty(len(channels))
    for i, c in enumerate(channels):
        device = campaign["device"][c]
        if v_total[i] > noise_v[i] and total[i] > 0:
            noise[i] = total[i] / (local_exponent(device, total[i] / 1000.0) * v_total[i]) * noise_v[i]
        else:
            a, b = NOISE_FIT[device]
            noise[i] = (noise_v[i] / a) ** (1.0 / b)
    dose = np.maximum(dose, 0.0)
    sigma = np.hypot(campaign["systematic"] * dose, noise)
    return dose, sigma, noise


def shift_measurements(campaign: dict, phase: str):
    """Shift space: phase-end shift [V], its endpoint noise [V], the exponent and the knee dose [rad] per channel."""
    channels = campaign["channels"]
    table = campaign["phases"][phase]["table"].loc[channels]
    shift = table["shift_end_v"].to_numpy(float)
    noise = np.hypot(table["fit_rms_mv"].to_numpy(float) / 1000.0, campaign["noise_floor"])
    exponent = np.array([campaign["exponent"][campaign["device"][c]] for c in channels])
    offset = np.array([campaign["offset_krad"][campaign["device"][c]] * 1000.0 for c in channels])
    return shift, noise, exponent, offset


def measurement(campaign: dict, phase: str, previous_model=None) -> dict:
    """The quantity fitted for one phase, in the campaign's fit space.

    Returns a dict with 'space', 'value', 'noise', 'unit', and in dose space
    'sigma' (the total uncertainty), in shift space 'exponent' and 'offset'
    (per channel), 'systematic', 'log_fluence' (which turns the per-fluence
    kernel into a dose) and 'previous' (the previous phase's unscaled model
    dose per channel [rad], zeros for a first phase).
    """
    if campaign["space"] == "dose":
        dose, sigma, noise = measurements(campaign, phase)
        return {"space": "dose", "value": dose, "sigma": sigma, "noise": noise, "unit": "rad"}
    shift, noise, exponent, offset = shift_measurements(campaign, phase)
    previous = np.zeros(len(shift)) if previous_model is None else np.asarray(previous_model, float)
    tables = campaign.get("response_table", {})
    response = [tables.get(campaign["device"][c]) for c in campaign["channels"]] if tables else None
    return {"space": "shift", "value": shift, "noise": noise, "exponent": exponent, "offset": offset,
            "response": response, "systematic": campaign["systematic"], "previous": previous,
            "log_fluence": campaign["phases"][phase]["log_fluence"], "unit": "V"}


def phase_models(campaign: dict, phase: str) -> tuple[str, ...]:
    """The phase's default scale model followed by its diagnostic models."""
    cfg = campaign["phases"][phase]
    return (cfg["scale_model"], *cfg["diagnostic_models"])


def scale_groups(campaign: dict, channels: list[int], model: str = "shared") -> dict[str, np.ndarray]:
    """Indices into `channels` that share one calibration scale, by model.

    In shift space the amplitudes are always split by device type, so
    `shared` is `per-type` there and `fading` is `fading-per-type`.
    """
    if campaign["space"] == "shift":
        model = {"shared": "per-type", "fading": "fading-per-type"}.get(model, model)
    if model == "shared":
        return {"all": np.arange(len(channels))}
    if model == "fading":
        labels = {c: "faded" if c in campaign["faded_channels"] else "not faded" for c in channels}
    elif model == "per-type":
        labels = {c: campaign["device"][c] for c in channels}
    elif model == "fading-per-type":
        labels = {c: ("faded " if c in campaign["faded_channels"] else "fresh ") + campaign["device"][c]
                  for c in channels}
    else:
        raise ValueError(f"unknown scale model {model!r}; use one of {SCALE_MODELS}")
    groups = {}
    for i, c in enumerate(channels):
        groups.setdefault(labels[c], []).append(i)
    return {k: np.array(v) for k, v in groups.items()}


def solve(meas: dict, model, model_err, groups):
    """Fit the group scales of one model in the measurement's space.

    Returns (scales by group, chi2, pulls, predicted per channel); the scales
    are dose scales in dose space and amplitudes A [V / rad^n] in shift space.
    """
    if meas["space"] == "dose":
        scales, chi2, pulls, scale = ek.fitScales(meas["value"], meas["sigma"], model, model_err, groups)
        return scales, chi2, pulls, scale * model
    f = meas["log_fluence"]
    return ek.fitShiftAmplitudes(meas["value"], meas["noise"], meas["exponent"], meas["offset"],
                                 meas["systematic"], model * f, model_err * f, groups,
                                 previous=meas["previous"], response=meas.get("response"))


def scale_report(meas: dict, scales: dict, log_fluence: float) -> dict:
    """The fitted scales as reported: calibration ratios (dose space) or amplitudes (shift)."""
    if meas["space"] == "dose":
        return {"calibration_ratio": {g: s / log_fluence for g, s in scales.items()}}
    return {"amplitude_v_per_rad_pow_n": dict(scales)}


def fold(energies, kernel, error, mus, sigma=None, cut=None):
    """Fold weighted kernel columns with the published width at each centre energy.

    A given `sigma` [MeV] replaces the published width by a fixed one; `cut`
    cuts each Gaussian off above mu + cut * sigma (Dependencies.EnergyKernel.foldKernel).
    """
    mus = np.asarray(mus, dtype=float)
    sigmas = ek.publishedSkandionSigma(mus) if sigma is None else float(sigma)
    return ek.foldKernel(energies, kernel, error, mus, sigmas, step=FINE_STEP, upperCut=cut)


def gaussian_grid(energies, kernel, error, mus, sigmas):
    """Fold weighted kernel columns with Gaussians on a (mu, sigma) grid: (mu, sigma, tally)."""
    mus = np.asarray(mus, dtype=float)
    sigmas = np.asarray(sigmas, dtype=float)
    means, errs = ek.foldKernel(energies, kernel, error, np.repeat(mus, len(sigmas)),
                                np.tile(sigmas, len(mus)), step=FINE_STEP)
    shape = (len(mus), len(sigmas), kernel.shape[1])
    return means.reshape(shape), errs.reshape(shape)


def weighted_kernel(kernel_path: Path, campaign: dict, phase: str, channels: list[int]):
    """Kernel energies, doses and errors for `channels`, times each channel's fluence weight."""
    energies, kernel, error = load_kernel(kernel_path)
    w = channel_weights(campaign, phase, channels)
    return energies, kernel[:, channels] * w, error[:, channels] * w


def fit_model(meas, energies, k, e, mus, mu_range, groups, log_fluence):
    """Mono and Gaussian fits for one calibration-scale model; returns (summary, arrays)."""
    lo, hi = mu_range
    inside = (energies >= lo) & (energies <= hi)
    if not inside.any():
        raise ValueError(f"no kernel energy inside {lo}-{hi} MeV")
    mono_chi2 = np.array([solve(meas, k[i], e[i], groups)[1] if inside[i] else np.inf
                          for i in range(len(energies))])
    m = int(np.argmin(mono_chi2))
    mono_scales, _, mono_pulls, mono_predicted = solve(meas, k[m], e[m], groups)

    means, errs = gaussian_grid(energies, k, e, mus, SIGMAS)
    chi2 = np.empty((len(mus), len(SIGMAS)))
    for i in range(len(mus)):
        for j in range(len(SIGMAS)):
            chi2[i, j] = solve(meas, means[i, j], errs[i, j], groups)[1]
    i, j = np.unravel_index(int(np.argmin(chi2)), chi2.shape)
    scales, best_chi2, pulls, predicted = solve(meas, means[i, j], errs[i, j], groups)
    mu_ok = mus[chi2.min(axis=1) <= best_chi2 + 1.0]
    sigma_ok = SIGMAS[chi2.min(axis=0) <= best_chi2 + 1.0]

    pw_means, pw_errs = fold(energies, k, e, mus)
    pw_chi2 = np.array([solve(meas, pw_means[p], pw_errs[p], groups)[1] for p in range(len(mus))])
    p = int(np.argmin(pw_chi2))
    pw_scales, _, pw_pulls, pw_predicted = solve(meas, pw_means[p], pw_errs[p], groups)
    pw_ok = mus[pw_chi2 <= pw_chi2[p] + 1.0]

    n = len(meas["value"])
    summary = {
        "mono": {
            "energy_MeV": float(energies[m]), "chi2": float(mono_chi2[m]), "dof": n - 1 - len(groups),
            **scale_report(meas, mono_scales, log_fluence),
        },
        "gaussian": {
            "mean_energy_MeV": float(mus[i]), "sigma_energy_MeV": float(SIGMAS[j]),
            "chi2": float(best_chi2), "dof": n - 2 - len(groups),
            "delta_chi2_vs_mono": float(best_chi2 - mono_chi2[m]),
            "mean_energy_interval_MeV": [float(mu_ok.min()), float(mu_ok.max())],
            "sigma_energy_interval_MeV": [float(sigma_ok.min()), float(sigma_ok.max())],
            **scale_report(meas, scales, log_fluence),
        },
        "published_width": {
            "mean_energy_MeV": float(mus[p]), "sigma_energy_MeV": published_sigma(float(mus[p])),
            "chi2": float(pw_chi2[p]), "dof": n - 1 - len(groups),
            "mean_energy_interval_MeV": [float(pw_ok.min()), float(pw_ok.max())],
            **scale_report(meas, pw_scales, log_fluence),
        },
    }
    arrays = {
        "objective": chi2,
        "mono_predicted": mono_predicted, "mono_pull": mono_pulls,
        "gaussian_predicted": predicted, "gaussian_pull": pulls,
        "published_predicted": pw_predicted, "published_pull": pw_pulls,
        # The unscaled model dose [rad] at the published-width best energy: what
        # a later phase on the same devices starts from in shift space.
        "published_model_rad": pw_means[p] * log_fluence,
    }
    return summary, arrays


def fit_fixed_width(meas, energies, k, e, mus, groups, log_fluence, spec):
    """Centre-energy scan of one fixed-width model (fixed_width_specs); returns (summary, arrays)."""
    means, errs = fold(energies, k, e, mus, spec["sigma"], spec["cut"])
    chi2 = np.array([solve(meas, means[p], errs[p], groups)[1] for p in range(len(mus))])
    p = int(np.argmin(chi2))
    scales, _, pulls, predicted = solve(meas, means[p], errs[p], groups)
    ok = mus[chi2 <= chi2[p] + 1.0]
    summary = {
        "label": spec["label"],
        "mean_energy_MeV": float(mus[p]), "sigma_energy_MeV": spec["sigma"],
        "upper_cut_sigma": spec["cut"],
        "spectrum_mean_energy_MeV": spectrum_mean(float(mus[p]), spec["sigma"], spec["cut"]),
        "chi2": float(chi2[p]), "dof": len(meas["value"]) - 1 - len(groups),
        "mean_energy_interval_MeV": [float(ok.min()), float(ok.max())],
        "at_edge": bool(p in (0, len(mus) - 1)),
        **scale_report(meas, scales, log_fluence),
    }
    arrays = {"chi2": chi2, "predicted": predicted, "pull": pulls, "model_rad": means[p] * log_fluence}
    return summary, arrays


def fit_phase(
    campaign: dict,
    label: str,
    kernel_path: Path,
    mu_range: tuple[float, float],
    previous_model=None,
    fixed_widths=(),
    previous_fixed=None,
) -> tuple[dict, pd.DataFrame, dict]:
    """Fits of one phase to all channels: default scale model plus diagnostics.

    `previous_model` is the previous phase's unscaled model dose per channel
    [rad] at its fitted energy; it enters the shift-space model of a phase on
    pre-dosed devices and is ignored in dose space. `fixed_widths` are the
    optional fixed-width models (fixed_width_specs), fitted with the default
    scale model; `previous_fixed` gives, by spec key, the previous phase's
    model dose of that same model (`previous_model` when missing).
    """
    cfg = campaign["phases"][label]
    channels = campaign["channels"]
    meas = measurement(campaign, label, previous_model)
    energies, k, e = weighted_kernel(kernel_path, campaign, label, channels)
    weights = channel_weights(campaign, label, channels)
    lo, hi = mu_range
    mus = np.arange(lo, hi + FINE_STEP / 2, FINE_STEP)
    default, *diagnostic_models = phase_models(campaign, label)
    fits = {model: fit_model(meas, energies, k, e, mus, mu_range,
                             scale_groups(campaign, channels, model), cfg["log_fluence"])
            for model in (default, *diagnostic_models)}
    main, main_arrays = fits[default]
    fixed = {}
    for spec in fixed_widths:
        previous = (previous_fixed or {}).get(spec["key"], previous_model)
        fixed_meas = measurement(campaign, label, previous)
        fixed[spec["key"]] = (spec, fixed_meas, *fit_fixed_width(
            fixed_meas, energies, k, e, mus, scale_groups(campaign, channels, default), cfg["log_fluence"], spec))

    unit = meas["unit"]
    scale = 1000.0 if unit == "rad" else 1.0
    shown = "krad" if unit == "rad" else "V"
    sigma = meas["sigma"] if meas["space"] == "dose" else meas["noise"]
    table = cfg["table"]
    rows = []
    for c in channels:
        row = {
            "phase": label, "channel": c, "group": table.loc[c, "group"], "sample": table.loc[c, "sample"],
            "fit_space": meas["space"],
            "delivered_fluence_protons_cm2": weights[c] * cfg["log_fluence"],
            f"measured_{shown}": meas["value"][c] / scale, f"sigma_{shown}": sigma[c] / scale,
            f"mono_predicted_{shown}": main_arrays["mono_predicted"][c] / scale,
            "mono_pull": main_arrays["mono_pull"][c],
            f"gaussian_predicted_{shown}": main_arrays["gaussian_predicted"][c] / scale,
            "gaussian_pull": main_arrays["gaussian_pull"][c],
            f"published_width_predicted_{shown}": main_arrays["published_predicted"][c] / scale,
            "published_width_pull": main_arrays["published_pull"][c],
        }
        if meas["space"] == "shift":
            row["previous_model_dose_krad"] = meas["previous"][c] / 1000.0
            row["exponent"] = meas["exponent"][c]
            row["offset_krad"] = meas["offset"][c] / 1000.0
        for model in diagnostic_models:
            key = model.replace("-", "_")
            row[f"{key}_gaussian_predicted_{shown}"] = fits[model][1]["gaussian_predicted"][c] / scale
            row[f"{key}_gaussian_pull"] = fits[model][1]["gaussian_pull"][c]
        for key, (spec, fixed_meas, _, fixed_arrays) in fixed.items():
            row[f"{key}_predicted_{shown}"] = fixed_arrays["predicted"][c] / scale
            row[f"{key}_pull"] = fixed_arrays["pull"][c]
            row[f"{key}_model_dose_krad"] = fixed_arrays["model_rad"][c] / 1000.0
            if meas["space"] == "shift":
                row[f"{key}_previous_model_dose_krad"] = fixed_meas["previous"][c] / 1000.0
        if fixed:
            row["published_width_model_dose_krad"] = main_arrays["published_model_rad"][c] / 1000.0
        rows.append(row)
    summary = {
        "phase": label,
        "fit_space": meas["space"],
        "systematic_fraction": campaign["systematic"],
        "noise_floor_mv": 1e3 * campaign["noise_floor"],
        "channels": channels,
        "nominal_energy_MeV": cfg["nominal_energy"],
        "mu_range_MeV": [lo, hi],
        "scale_model": default,
        "scale_groups": {g: [channels[i] for i in idx]
                         for g, idx in scale_groups(campaign, channels, default).items()},
        **main,
        **{key: fit[2] for key, fit in fixed.items()},
        "diagnostics": {model: fits[model][0] for model in diagnostic_models},
    }
    if meas["space"] == "shift":
        types = sorted(set(campaign["device"].values()))
        summary["exponent"] = {t: campaign["exponent"][t] for t in types}
        summary["offset_krad"] = {t: campaign["offset_krad"][t] for t in types}
    arrays = {
        "mus": mus, "objective": main_arrays["objective"], "unit": shown,
        "channels": np.asarray(channels), "device": [campaign["device"][c] for c in channels],
        "measured": meas["value"] / scale, "sigma": sigma / scale,
        "mono_predicted": main_arrays["mono_predicted"] / scale,
        "gaussian_predicted": main_arrays["gaussian_predicted"] / scale,
        "published_predicted": main_arrays["published_predicted"] / scale,
        "published_model_rad": main_arrays["published_model_rad"],
        "diagnostic_gaussian_predicted": {
            model: fits[model][1]["gaussian_predicted"] / scale for model in diagnostic_models},
        # Optional fixed-width models: spec, prediction (shown unit) and the unscaled
        # model dose [rad] a later phase of the same model starts from.
        "fixed": {key: {"spec": spec, "predicted": fixed_arrays["predicted"] / scale,
                        "model_rad": fixed_arrays["model_rad"]}
                  for key, (spec, _, _, fixed_arrays) in fixed.items()},
    }
    return summary, pd.DataFrame(rows), arrays


def plot_results(results: dict[str, tuple[dict, pd.DataFrame, dict]], output: Path) -> None:
    fig, axes = plt.subplots(2, len(results), figsize=(6.0 * len(results), 8.2),
                             constrained_layout=True, squeeze=False)
    for column, (label, (summary, rows, arrays)) in enumerate(results.items()):
        ax = axes[0, column]
        x = arrays["channels"]
        unit = arrays["unit"]
        device = np.asarray(arrays["device"])
        for dev in sorted(set(device)):
            sel = device == dev
            y, s = arrays["measured"][sel], arrays["sigma"][sel]
            ax.errorbar(x[sel], y, yerr=[np.minimum(s, 0.95 * np.abs(y)), s], fmt=MARKER[dev], color="C0",
                        capsize=3, label=f"RadFET measured, {dev}")
        mono, gauss, pw = summary["mono"], summary["gaussian"], summary["published_width"]
        model = summary["scale_model"]
        curves = [
            (arrays["mono_predicted"], "s", "C1", "--",
             f"best monoenergy {mono['energy_MeV']:.2f} MeV, chi2 {mono['chi2']:.1f}/{mono['dof']}"),
            (arrays["gaussian_predicted"], "^", "C2", "--",
             f"best Gaussian {gauss['mean_energy_MeV']:.2f} MeV, sigma {gauss['sigma_energy_MeV']:.2f}, "
             f"chi2 {gauss['chi2']:.1f}/{gauss['dof']}"),
            (arrays["published_predicted"], "D", "C3", "-",
             f"best Gaussian, published width: {pw['mean_energy_MeV']:.2f} MeV, "
             f"sigma {pw['sigma_energy_MeV']:.3f}, chi2 {pw['chi2']:.1f}/{pw['dof']}"),
        ]
        for (diag_model, diag), colour in zip(summary["diagnostics"].items(), DIAGNOSTIC_COLOURS):
            g = diag["gaussian"]
            curves.append((arrays["diagnostic_gaussian_predicted"][diag_model], "x", colour, ":",
                           f"diagnostic, {diag_model} scale: {g['mean_energy_MeV']:.2f} MeV, "
                           f"chi2 {g['chi2']:.1f}/{g['dof']}"))
        for key, fixed in arrays["fixed"].items():
            f = summary[key]
            curves.append((fixed["predicted"], "v", fixed["spec"]["colour"], "-.",
                           f"best Gaussian, {fixed['spec']['label']}: {f['mean_energy_MeV']:.2f} MeV, "
                           f"chi2 {f['chi2']:.1f}/{f['dof']}"))
        for y, marker, colour, style, name in curves:
            ax.plot(x, np.where(y > 0, y, np.nan), marker + style, color=colour, mfc="none", label=name)
        ax.set_yscale("log")
        ax.set_xticks(x, [str(v) for v in x])
        ax.set_xlabel("RadFET channel")
        ax.set_ylabel("Phase dose (krad)" if unit == "krad" else "Phase-end threshold shift (V)")
        ax.set_title(f"{label.replace('mev', ' MeV phase')}: all {len(x)} channels, {model} scale, "
                     f"{summary['fit_space']} space")
        ax.grid(alpha=0.25, which="both")
        ax.legend(fontsize=8)

        ax = axes[1, column]
        delta = arrays["objective"] - arrays["objective"].min()
        image = ax.pcolormesh(arrays["mus"], SIGMAS, delta.T, shading="auto",
                              cmap="viridis_r", vmin=0, vmax=25, rasterized=True)
        ax.plot(arrays["mus"], [published_sigma(mu) for mu in arrays["mus"]], "w--", lw=1,
                label="published width")
        ax.plot(gauss["mean_energy_MeV"], gauss["sigma_energy_MeV"], "r*", markersize=12,
                label="best Gaussian")
        ax.plot(pw["mean_energy_MeV"], pw["sigma_energy_MeV"], "D", color="C3", mfc="none", markersize=8,
                mew=2, label="best Gaussian, published width")
        for (diag_model, diag), colour in zip(summary["diagnostics"].items(), DIAGNOSTIC_COLOURS):
            g = diag["gaussian"]
            ax.plot(g["mean_energy_MeV"], g["sigma_energy_MeV"], "x", color=colour, markersize=9, mew=2,
                    label=f"{diag_model} diagnostic")
        lines = set()
        for key, fixed in arrays["fixed"].items():
            spec, f = fixed["spec"], summary[key]
            if spec["sigma"] not in lines:
                lines.add(spec["sigma"])
                ax.axhline(spec["sigma"], color=spec["colour"], ls=":", lw=1.5,
                           label=f"sigma = {spec['sigma']:g} MeV (fixed width)")
            ax.plot(f["mean_energy_MeV"], spec["sigma"], "v", color=spec["colour"], mfc="none", markersize=9,
                    mew=2, label=f"best Gaussian, {spec['label']}"
                    + (" (map: full Gaussians)" if spec["cut"] is not None else ""))
        ax.set_xlim(*summary["mu_range_MeV"])
        ax.set_xlabel("Gaussian mean energy (MeV)")
        ax.set_ylabel("Gaussian sigma (MeV)")
        ax.set_title(f"Weighted 12-channel misfit, {model} scale")
        ax.legend(fontsize=8, loc="upper right")
        fig.colorbar(image, ax=ax, label="Delta chi2 (clipped at 25)")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def add_mu_range_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--mu-range", nargs=3, action="append", default=[],
                        metavar=("PHASE", "LO", "HI"),
                        help="searched centre energies [MeV] for one phase, overriding the "
                             "campaign's mu_range_MeV; repeatable")


def add_space_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--space", choices=FIT_SPACES, default=None,
                        help="fit space, overriding the campaign's fit_space")
    parser.add_argument("--systematic", type=float, default=None,
                        help="device-spread fraction, overriding the campaign's systematic_fraction")
    parser.add_argument("--noise-floor-mv", type=float, default=None,
                        help="noise floor [mV] added in quadrature to every endpoint RMS, "
                             "overriding the campaign's noise_floor_mv")


def mu_range_overrides(campaign: dict, overrides: list[list[str]]) -> dict[str, tuple[float, float]]:
    """Each phase's searched centre-energy range, with --mu-range overrides applied."""
    ranges = {label: cfg["mu_range"] for label, cfg in campaign["phases"].items()}
    for label, lo, hi in overrides:
        if label not in ranges:
            raise SystemExit(f"--mu-range: unknown phase {label!r}; the campaign has {list(ranges)}")
        ranges[label] = (float(lo), float(hi))
    return ranges


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--kernel", type=Path, required=True,
                        help="unweighted kernel summary from BuildKernelSummary")
    parser.add_argument("--campaign", type=Path, required=True,
                        help="campaign JSON next to the measured doses")
    parser.add_argument("--output", type=Path, required=True,
                        help="output stem, e.g. Results/spectrum-fit-oxide")
    add_mu_range_argument(parser)
    add_space_argument(parser)
    add_width_argument(parser)
    args = parser.parse_args()
    specs = fixed_width_specs(args.fixed_width_mev, args.fixed_width_cut_sigma)
    campaign = load_campaign(args.campaign, args.space, args.systematic, args.noise_floor_mv)
    mu_ranges = mu_range_overrides(campaign, args.mu_range)
    stem = args.output
    stem.parent.mkdir(parents=True, exist_ok=True)

    # Phases in order: a phase on pre-dosed devices starts, in shift space,
    # from the previous phase's unscaled model dose at its fitted energy
    # (each fixed-width model from its own previous model dose).
    results = {}
    model_by_phase = {}
    fixed_by_phase = {}
    for label, cfg in campaign["phases"].items():
        previous = model_by_phase.get(cfg["previous_phase"]) if cfg["previous_phase"] else None
        previous_fixed = fixed_by_phase.get(cfg["previous_phase"]) if cfg["previous_phase"] else None
        results[label] = fit_phase(campaign, label, args.kernel, mu_ranges[label], previous,
                                   specs, previous_fixed)
        model_by_phase[label] = results[label][2]["published_model_rad"] + (
            previous if previous is not None else 0.0)
        fixed_by_phase[label] = {
            key: fixed["model_rad"] + (previous_fixed[key] if previous_fixed else 0.0)
            for key, fixed in results[label][2]["fixed"].items()}

    summary = {label: result[0] for label, result in results.items()}
    stem.with_name(f"{stem.name}-summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    pd.concat([result[1] for result in results.values()], ignore_index=True).to_csv(
        stem.with_name(f"{stem.name}-channel-residuals.csv"), index=False
    )
    plot_results(results, stem.with_name(f"{stem.name}-diagnostics.pdf"))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
