"""Shared core of the RadFET kernel fits: campaign file, measurement models, scale groups, solver.

Every kernel or spectrum fit of a RadEx campaign builds on this module: the proton fits
of Uppsala (RadEx.InferSpectrum, FitAllChannels, AluminiumRowEnergyScan) and the
electron fits of HUS (RadEx.FitScaledSpectrum, FitSpectrumModels, JointSpectrumFit,
CompareVariant, BuildCoefficientTable). A fit tool chooses its beam models; the
campaign, the measured quantity and its errors, the calibration-scale groups and the
scale or amplitude solve are the same for all of them. The folding and the solvers
themselves are Dependencies.EnergyKernel.

A campaign is a JSON file kept next to the measured tables (for Uppsala,
GRAS/3D/RadEx/RadEx-Uppsala/MeasuredDoses/campaign.json; for HUS,
GRAS/3D/RadEx/RadEx-HUS/MeasuredDoses/campaign.json). Relative paths resolve against
it. Per phase it names:

    dose                the measured table (or the campaign's `dose` for all phases),
                        one row per channel; `rows` ({column: value}) selects the
                        phase's rows of a shared table
    dose_column         the dose fitted
    noise_column, conversion_column
                        optional: the table's own per-channel dose errors; sigma^2 =
                        noise^2 + conversion^2 + (systematic_fraction x dose)^2 in
                        the table's unit. Without them the dose and its noise come
                        from the phase-end shift through the vendor calibration
                        (InferSpectrum's dose space).
    log_fluence_protons_cm2 or log_fluence_column
                        the phase fluence (a number, or the table column holding it)
    fluence             optional per-channel delivered fluence; uniform without it
    previous_phase      the phase this one increments on, if any
    nominal_energy_MeV, mu_range_MeV
                        the nominal energy and the searched centre energies
    scale_model, diagnostic_models
                        the default calibration-scale model and the diagnostic ones,
                        optionally per fit space as {"dose": ..., "shift": ...}
    note                free text

and for the campaign: `fit_space`, `systematic_fraction`, `device_column` (the
table column with each channel's device type, default `group`), and for the shift
space `shift_exponent`, `shift_offset_krad`, `shift_response_table`, `noise_floor_mv`
and `faded_channels` (InferSpectrum's docstring describes the two fit spaces and the
scale models). Every other key (for HUS: `labels`, `layers`, `default_layer`,
`measurement_phase`, `spectrum`, `reference_fluence`, `joint_fit`, ...) is passed
through for the tools that use it.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from Dependencies import EnergyKernel as ek

MARKER = {"VT01": "o", "VT05": "s"}
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




def _by_space(setting, space: str):
    """A campaign setting that may be given per fit space, as {"dose": ..., "shift": ...}."""
    return setting[space] if isinstance(setting, dict) else setting


#: Campaign-file keys load_campaign interprets; every other key is passed through unchanged
#: (with the paths among them resolved), for the tools that need it.
_CORE_KEYS = {"phases", "fit_space", "systematic_fraction", "noise_floor_mv", "shift_exponent",
              "shift_offset_krad", "shift_response_table", "faded_channels", "dose", "device_column"}


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
    path = Path(path).resolve()
    data = json.loads(path.read_text(encoding="utf-8"))
    device_column = data.get("device_column", "group")
    space = space or data.get("fit_space", "dose")
    if space not in FIT_SPACES:
        raise ValueError(f"{path}: fit_space {space!r}; use one of {FIT_SPACES}")
    phases = {}
    for label, cfg in data["phases"].items():
        previous = cfg.get("previous_phase")
        if previous is not None and previous not in phases:
            raise ValueError(f"{path}: {label} increments on {previous!r}, which must be listed before it")
        models = (_by_space(cfg.get("scale_model", "shared"), space), *_by_space(cfg.get("diagnostic_models", ()), space))
        unknown = [m for m in models if m not in SCALE_MODELS]
        if unknown:
            raise ValueError(f"{path}: {label} has unknown scale models {unknown}; use {SCALE_MODELS}")
        dose_path = path.parent / cfg.get("dose", data.get("dose"))
        # round_trip parses every value exactly, as float() does; the default parser can be off in the last bit.
        table = pd.read_csv(dose_path, float_precision="round_trip")
        for column, value in cfg.get("rows", {}).items():
            table = table[table[column].astype(str) == str(value)]
        table = table.set_index("channel")
        if "log_fluence_column" in cfg:
            log_fluence = float(table[cfg["log_fluence_column"]].iloc[0])
        else:
            log_fluence = float(cfg["log_fluence_protons_cm2"])
        phases[label] = {
            "dose": dose_path,
            "table": table,
            "dose_column": cfg["dose_column"],
            "noise_column": cfg.get("noise_column"),
            "conversion_column": cfg.get("conversion_column"),
            "previous_phase": previous,
            "fluence": path.parent / cfg["fluence"] if "fluence" in cfg else None,
            "log_fluence": log_fluence,
            "log_fluence_column": cfg.get("log_fluence_column"),
            "nominal_energy": float(cfg["nominal_energy_MeV"]) if "nominal_energy_MeV" in cfg else None,
            "mu_range": tuple(float(v) for v in cfg["mu_range_MeV"]) if "mu_range_MeV" in cfg else None,
            "scale_model": models[0],
            "diagnostic_models": models[1:],
            "note": cfg.get("note", ""),
        }
    if not phases:
        raise ValueError(f"{path}: no phases")
    first = next(iter(phases.values()))["table"]
    channels = sorted(int(c) for c in first.index)
    device = {c: str(first.loc[c, device_column]) for c in channels}
    for label, cfg in phases.items():
        table = cfg["table"]
        if sorted(int(c) for c in table.index) != channels or any(
                str(table.loc[c, device_column]) != device[c] for c in channels):
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
    def resolve(value):
        return None if value is None else (path.parent / value).resolve()

    settings = {k: v for k, v in data.items() if k not in _CORE_KEYS}
    for key in ("spectrum", "reference_run"):
        if key in settings:
            settings[key] = resolve(settings[key])
    if "joint_fit" in settings:
        settings["joint_fit"] = dict(settings["joint_fit"])
        for key in ("response_matrix", "measured_pdd"):
            settings["joint_fit"][key] = resolve(settings["joint_fit"][key])
    if "labels" in settings:
        settings["labels"] = {int(c): v for c, v in settings["labels"].items()}
    return {
        **settings,
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
    if cfg["fluence"] is None:
        return np.ones(len(channels))
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
    if cfg["noise_column"] is not None:
        # The table carries each channel's dose errors itself, in its own unit.
        dose = table[cfg["dose_column"]].to_numpy(float)
        noise = table[cfg["noise_column"]].fillna(0).to_numpy(float)
        conversion = (table[cfg["conversion_column"]].fillna(0).to_numpy(float)
                      if cfg["conversion_column"] is not None else np.zeros(len(channels)))
        sigma = np.sqrt(noise ** 2 + conversion ** 2 + (campaign["systematic"] * dose) ** 2)
        return dose, sigma, np.hypot(noise, conversion)
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
        unit = "table" if campaign["phases"][phase]["noise_column"] is not None else "rad"
        return {"space": "dose", "value": dose, "sigma": sigma, "noise": noise, "unit": unit}
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


def scale_groups(campaign: dict, channels: list[int], model: str = "shared",
                 exclude=()) -> dict[str, np.ndarray]:
    """Indices into `channels` that share one calibration scale, by model.

    In shift space the amplitudes are always split by device type, so
    `shared` is `per-type` there and `fading` is `fading-per-type`.
    """
    if campaign["space"] == "shift":
        model = {"shared": "per-type", "fading": "fading-per-type"}.get(model, model)
    if model == "shared":
        return {"all": np.array([i for i, c in enumerate(channels) if c not in exclude])}
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
        if c not in exclude:
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


def gaussian_grid(energies, kernel, error, mus, sigmas):
    """Fold weighted kernel columns with Gaussians on a (mu, sigma) grid: (mu, sigma, tally)."""
    mus = np.asarray(mus, dtype=float)
    sigmas = np.asarray(sigmas, dtype=float)
    means, errs = ek.foldKernel(energies, kernel, error, np.repeat(mus, len(sigmas)),
                                np.tile(sigmas, len(mus)), step=FINE_STEP)
    shape = (len(mus), len(sigmas), kernel.shape[1])
    return means.reshape(shape), errs.reshape(shape)


def add_space_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--space", choices=FIT_SPACES, default=None,
                        help="fit space, overriding the campaign's fit_space")
    parser.add_argument("--systematic", type=float, default=None,
                        help="device-spread fraction, overriding the campaign's systematic_fraction")
    parser.add_argument("--noise-floor-mv", type=float, default=None,
                        help="noise floor [mV] added in quadrature to every endpoint RMS, "
                             "overriding the campaign's noise_floor_mv")


def add_mu_range_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--mu-range", nargs=3, action="append", default=[],
                        metavar=("PHASE", "LO", "HI"),
                        help="searched centre energies [MeV] for one phase, overriding the "
                             "campaign's mu_range_MeV; repeatable")


def mu_range_overrides(campaign: dict, overrides: list[list[str]]) -> dict[str, tuple[float, float]]:
    """Each phase's searched centre-energy range, with --mu-range overrides applied."""
    ranges = {label: cfg["mu_range"] for label, cfg in campaign["phases"].items()}
    for label, lo, hi in overrides:
        if label not in ranges:
            raise SystemExit(f"--mu-range: unknown phase {label!r}; the campaign has {list(ranges)}")
        ranges[label] = (float(lo), float(hi))
    return ranges


def layer_slice(campaign: dict, layer: str) -> slice:
    """Tally columns of one detector layer of the GRAS runs (the campaign's `layers`)."""
    first = campaign["layers"][layer]
    return slice(first, first + len(campaign["channels"]))


def device_types(campaign: dict) -> list[str]:
    """Device types in order of first appearance."""
    return list(dict.fromkeys(campaign["device"][c] for c in campaign["channels"]))
