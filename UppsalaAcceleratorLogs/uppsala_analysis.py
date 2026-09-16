#!/usr/bin/env python3
"""Streaming analysis of IBA PBS data-recorder logs from the Uppsala campaign.

The program deliberately separates three stages:

* ``inventory`` reads file metadata and builds a reproducible manifest.
* ``parse`` streams assigned map_record files and creates charge-domain products.
* ``reduce`` applies the facility MU calibration and combines parsed products.

No input file is modified.  The charge-domain products preserve enough information
to re-run calibration/reduction without re-reading the multi-GB source CSV files.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np


SCHEMA_VERSION = 1
DEFAULT_RANGE_ENERGY = {
    # Logged water-equivalent range [cm] -> nominal facility energy [MeV].
    # 64 and 85 MeV are directly identified by the facility experiment log.
    # The remaining values are matched to the nearest facility table energy
    # using the standard proton range-energy relation.
    3.2009159491194183: 60.0,
    3.5883448634201853: 64.0,
    4.2048897994367405: 70.0,
    5.923796604781645: 85.0,
    15.029766405465406: 145.0,
    15.933099216775974: 150.0,
    30.94264381823079: 220.0,
    32.15413781909827: 226.0,
}

# The BMS primary-monitor dose channel is calibrated to 3 nC per MU before the
# per-map K_FACTOR correction.  DOSE_PRIM is already an incremental charge in
# coulombs; it must be summed directly and must not be multiplied by the 250 us
# acquisition period.
MONITOR_CHARGE_C_PER_MU = 3.0e-9

# The facility supplied its particle conversion at discrete nominal energies.
# The campaign's "64 MeV" setting uses the facility's 65 MeV conversion row
# (the supplied 7.14E12 total is 125,000 MU × 5.71E7 protons/MU).
FACILITY_CALIBRATION_ENERGY_ALIASES = {64.0: 65.0}

# The facility was explicitly asked to irradiate the complete RadEx instrument
# surface with fluence uniform to within 5%. This fixed area, rather than
# geometric magnification inferred from recorder profiles, is the denominator.
RADEX_INSTRUMENT_WIDTH_X_MM = 154.0
RADEX_INSTRUMENT_WIDTH_Y_MM = 66.8
RADEX_INSTRUMENT_AREA_CM2 = RADEX_INSTRUMENT_WIDTH_X_MM * RADEX_INSTRUMENT_WIDTH_Y_MM / 100.0

# CRDS was irradiated with two requested target fields.  The first three runs
# used the smaller field and the final run used the larger field, so its fluence
# denominator must be selected per irradiation rather than once per aggregate.
CRDS_TARGETS = {
    "exp_10": {"width_x_cm": 5.2, "width_y_cm": 5.4, "area_cm2": 28.08},
    "exp_11": {"width_x_cm": 5.2, "width_y_cm": 5.4, "area_cm2": 28.08},
    "exp_12": {"width_x_cm": 5.2, "width_y_cm": 5.4, "area_cm2": 28.08},
    "exp_13": {"width_x_cm": 7.0, "width_y_cm": 7.4, "area_cm2": 51.8},
}
SESSION_RE = re.compile(r"skandion_irrlogs_(\d{8})_(\d{8})_(\d{6})_\d+\.PBS")
MAP_RE = re.compile(r"\.map_record_(\d{3})_(part|tuning)_(\d{2})\.csv$")


def json_dump(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    tmp.replace(path)


def parse_float(text: str | None) -> float | None:
    if text is None:
        return None
    try:
        value = float(text.strip())
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def parse_number(text: str) -> float:
    return float(text.strip().replace(",", ""))


def parse_experiment_csv(path: Path, date: str, utc_offset_hours: int) -> list[dict[str, Any]]:
    day = datetime.strptime(date, "%Y-%m-%d").date()
    tz = timezone(timedelta(hours=utc_offset_hours))
    experiments: list[dict[str, Any]] = []
    section = "Experiment 1"
    sequence = 0
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.reader(handle):
            row += [""] * (8 - len(row))
            first = row[0].strip()
            if first.lower().startswith("experiment"):
                section = first.strip() or section
                continue
            energy = parse_float(row[1])
            if energy is None or not row[2].strip() or not row[3].strip():
                continue
            try:
                mu = parse_number(row[2])
                start_clock = datetime.strptime(row[3].strip(), "%H:%M").time()
                stop_clock = datetime.strptime(row[4].strip(), "%H:%M").time()
            except (ValueError, IndexError):
                continue
            start = datetime.combine(day, start_clock, tzinfo=tz)
            stop = datetime.combine(day, stop_clock, tzinfo=tz)
            if stop < start:
                stop += timedelta(days=1)
            sequence += 1
            experiments.append(
                {
                    "id": f"exp_{sequence:02d}",
                    "section": section,
                    "label": first,
                    "energy_mev": energy,
                    "prescribed_mu": mu,
                    "start": start.isoformat(),
                    "stop": stop.isoformat(),
                }
            )
    if not experiments:
        raise ValueError(f"No irradiation rows found in {path}")
    return experiments


def parse_conversion_csv(path: Path) -> list[dict[str, float]]:
    values: list[dict[str, float]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            energy = parse_float(row.get("Energy"))
            ppmu = parse_float(row.get("Protons per MU"))
            if energy is not None and ppmu is not None:
                values.append({"energy_mev": energy, "protons_per_mu": ppmu})
    values.sort(key=lambda item: item["energy_mev"])
    if len(values) < 2:
        raise ValueError(f"Insufficient conversion values in {path}")
    return values


def interpolate_ppmu(energy: float, conversion: list[dict[str, float]]) -> float:
    energies = np.asarray([item["energy_mev"] for item in conversion], dtype=float)
    values = np.asarray([item["protons_per_mu"] for item in conversion], dtype=float)
    if energy < energies[0] or energy > energies[-1]:
        raise ValueError(f"Energy {energy} MeV is outside facility conversion table")
    return float(np.interp(energy, energies, values))


def facility_ppmu(
    energy: float, conversion: list[dict[str, float]]
) -> tuple[float, float]:
    calibration_energy = FACILITY_CALIBRATION_ENERGY_ALIASES.get(energy, energy)
    return interpolate_ppmu(calibration_energy, conversion), calibration_energy


def monitor_mu_from_charge(charge_c: float, k_factor: float) -> float:
    if k_factor <= 0:
        raise ValueError(f"Invalid K_FACTOR {k_factor}")
    return charge_c / (MONITOR_CHARGE_C_PER_MU * k_factor)


def monitor_protons_per_c(protons_per_mu: float, k_factor: float) -> float:
    if k_factor <= 0:
        raise ValueError(f"Invalid K_FACTOR {k_factor}")
    return protons_per_mu / (MONITOR_CHARGE_C_PER_MU * k_factor)


def session_datetime(name: str, utc_offset_hours: int) -> datetime:
    match = SESSION_RE.search(name)
    if not match:
        raise ValueError(f"Cannot parse session timestamp from {name}")
    date_text, repeated_date, clock = match.groups()
    if repeated_date != date_text:
        raise ValueError(f"Inconsistent dates in session name {name}")
    dt = datetime.strptime(date_text + clock, "%Y%m%d%H%M%S")
    return dt.replace(tzinfo=timezone(timedelta(hours=utc_offset_hours)))


def closest_experiment(
    session_time: datetime, experiments: list[dict[str, Any]], tolerance_minutes: float
) -> str | None:
    candidates: list[tuple[float, str]] = []
    for item in experiments:
        start = datetime.fromisoformat(item["start"])
        delta = abs((session_time - start).total_seconds())
        candidates.append((delta, item["id"]))
    delta, experiment_id = min(candidates)
    return experiment_id if delta <= tolerance_minutes * 60.0 else None


def read_spec_metadata(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle)
        version = next(reader, [""])[0]
        header = next(reader, [])
        values = next(reader, [])
    lookup = {name.lstrip("#"): idx for idx, name in enumerate(header)}

    def get(name: str) -> str | None:
        idx = lookup.get(name)
        return values[idx] if idx is not None and idx < len(values) else None

    return {
        "version": version,
        "layer_id": int(float(get("LAYER_ID") or 0)),
        "range_cm": parse_float(get("RANGE")),
        "planned_charge_c": parse_float(get("TOTAL_CHARGE")),
        "element_count": int(float(get("ELEMENT_COUNT") or 0)),
        "k_factor": parse_float(get("K_FACTOR")),
    }


def energy_from_range(range_cm: float | None) -> tuple[float | None, str]:
    if range_cm is None:
        return None, "missing"
    best_range = min(DEFAULT_RANGE_ENERGY, key=lambda candidate: abs(candidate - range_cm))
    if abs(best_range - range_cm) > 1e-5:
        return None, "unmapped"
    energy = DEFAULT_RANGE_ENERGY[best_range]
    source = "facility_experiment" if energy in (64.0, 85.0) else "range_inference"
    return energy, source


def build_inventory(args: argparse.Namespace) -> None:
    input_dir = Path(args.input_dir).resolve()
    experiment_csv = Path(args.experiment_csv).resolve()
    conversion_csv = Path(args.conversion_csv).resolve()
    experiments = parse_experiment_csv(
        experiment_csv, args.experiment_date, args.utc_offset_hours
    )
    conversion = parse_conversion_csv(conversion_csv)
    experiment_lookup = {item["id"]: item for item in experiments}
    for item in experiments:
        (
            item["protons_per_mu"],
            item["calibration_energy_mev"],
        ) = facility_ppmu(item["energy_mev"], conversion)
        item["calibration_value"] = item["prescribed_mu"] * item["protons_per_mu"]

    records = sorted(input_dir.glob("skandion_irrlogs_*/*.map_record_*.csv"))
    entries: list[dict[str, Any]] = []
    matched_specs: set[Path] = set()
    sessions: dict[str, dict[str, Any]] = {}
    for record in records:
        match = MAP_RE.search(record.name)
        if not match:
            continue
        map_number, kind, sequence = match.groups()
        spec = record.with_name(record.name.replace(".map_record_", ".map_specif_"))
        spec_metadata: dict[str, Any]
        if spec.exists():
            spec_metadata = read_spec_metadata(spec)
            matched_specs.add(spec.resolve())
        else:
            spec_metadata = {
                "version": None,
                "layer_id": None,
                "range_cm": None,
                "planned_charge_c": None,
                "element_count": None,
                "k_factor": None,
            }
        session = record.parent.name
        session_time = session_datetime(session, args.utc_offset_hours)
        experiment_id = closest_experiment(
            session_time, experiments, args.match_tolerance_minutes
        )
        energy, energy_source = energy_from_range(spec_metadata["range_cm"])
        if experiment_id is not None:
            documented_energy = experiment_lookup[experiment_id]["energy_mev"]
            if energy is not None and abs(energy - documented_energy) > 0.25:
                raise ValueError(
                    f"Range energy {energy} conflicts with experiment energy "
                    f"{documented_energy} for {record}"
                )
            energy = documented_energy
            energy_source = "facility_experiment"
        sessions.setdefault(
            session,
            {
                "session": session,
                "session_time": session_time.isoformat(),
                "experiment_id": experiment_id,
            },
        )
        entries.append(
            {
                "id": record.stem,
                "record_path": str(record.resolve()),
                "spec_path": str(spec.resolve()) if spec.exists() else None,
                "size_bytes": record.stat().st_size,
                "session": session,
                "session_time": session_time.isoformat(),
                "experiment_id": experiment_id,
                "map_number": int(map_number),
                "kind": kind,
                "sequence": int(sequence),
                "energy_mev": energy,
                "energy_source": energy_source,
                **spec_metadata,
            }
        )

    if not entries:
        raise ValueError(f"No map_record files found under {input_dir}")

    # Largest-processing-time-first assignment gives much better array balance than
    # index modulo because record sizes vary by four orders of magnitude.
    task_loads = [0] * args.tasks
    by_size = sorted(range(len(entries)), key=lambda idx: entries[idx]["size_bytes"], reverse=True)
    for idx in by_size:
        task = min(range(args.tasks), key=task_loads.__getitem__)
        entries[idx]["task"] = task
        task_loads[task] += entries[idx]["size_bytes"]
    entries.sort(key=lambda item: (item["session"], item["map_number"], item["kind"], item["sequence"]))

    all_specs = {path.resolve() for path in input_dir.glob("skandion_irrlogs_*/*.map_specif_*.csv")}
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "input_dir": str(input_dir),
        "experiment_csv": str(experiment_csv),
        "conversion_csv": str(conversion_csv),
        "timezone_offset_hours": args.utc_offset_hours,
        "tasks": args.tasks,
        "grid": {
            "minimum_mm": args.grid_min_mm,
            "maximum_mm": args.grid_max_mm,
            "pixel_mm": args.pixel_mm,
            "width_interpretation": args.width_interpretation,
            "width_quantization_mm": args.width_quantization_mm,
        },
        "conversion": conversion,
        "experiments": experiments,
        "sessions": sorted(sessions.values(), key=lambda item: item["session"]),
        "entries": entries,
        "orphan_spec_files": sorted(str(path) for path in all_specs - matched_specs),
        "task_load_bytes": task_loads,
    }
    json_dump(Path(args.output), manifest)
    assigned = sum(item["experiment_id"] is not None for item in entries)
    print(
        f"Wrote {args.output}: {len(entries)} records, {len(all_specs)} specs, "
        f"{assigned} records associated with documented experiments, "
        f"{len(all_specs - matched_specs)} orphan specs"
    )


def safe_float(row: list[str], idx: int) -> float | None:
    if idx >= len(row):
        return None
    return parse_float(row[idx])


def parse_iso_second(text: str) -> tuple[datetime, int]:
    dt = datetime.fromisoformat(text)
    if dt.tzinfo is None:
        raise ValueError(f"Timestamp lacks timezone: {text}")
    return dt, int(dt.timestamp())


def finalize_heatmap(
    width_histograms: dict[tuple[float, float], np.ndarray],
    pixel_mm: float,
    width_interpretation: str,
) -> np.ndarray:
    try:
        from scipy.ndimage import gaussian_filter
    except ImportError:  # Small dependency-free fallback, mainly for local tests.
        def gaussian_filter(
            values: np.ndarray,
            sigma: tuple[float, float],
            mode: str = "constant",
            truncate: float = 5.0,
        ) -> np.ndarray:
            del mode
            result = values
            for axis, axis_sigma in enumerate(sigma):
                radius = max(1, int(math.ceil(truncate * axis_sigma)))
                coordinates = np.arange(-radius, radius + 1, dtype=np.float64)
                kernel = np.exp(-0.5 * (coordinates / axis_sigma) ** 2)
                kernel /= kernel.sum()
                result = np.apply_along_axis(
                    lambda line: np.convolve(
                        np.pad(line, (radius, radius), mode="constant"),
                        kernel,
                        mode="valid",
                    ),
                    axis,
                    result,
                )
            return result

    if not width_histograms:
        return np.zeros((1, 1), dtype=np.float64)
    shape = next(iter(width_histograms.values())).shape
    result = np.zeros(shape, dtype=np.float64)
    for (x_width, y_width), histogram in width_histograms.items():
        input_sum = float(histogram.sum())
        if input_sum <= 0:
            continue
        if width_interpretation == "fwhm":
            x_sigma = x_width / 2.354820045
            y_sigma = y_width / 2.354820045
        else:
            x_sigma = x_width
            y_sigma = y_width
        blurred = gaussian_filter(
            histogram,
            sigma=(max(y_sigma / pixel_mm, 0.01), max(x_sigma / pixel_mm, 0.01)),
            mode="constant",
            truncate=5.0,
        )
        output_sum = float(blurred.sum())
        if output_sum > 0:
            blurred *= input_sum / output_sum
        result += blurred
    return result


def parse_record(entry: dict[str, Any], grid: dict[str, Any]) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    record_path = Path(entry["record_path"])
    grid_min = float(grid["minimum_mm"])
    grid_max = float(grid["maximum_mm"])
    pixel_mm = float(grid["pixel_mm"])
    quantization = float(grid["width_quantization_mm"])
    n_pixels = int(round((grid_max - grid_min) / pixel_mm))
    if n_pixels <= 0:
        raise ValueError("Invalid grid dimensions")

    second_values: dict[int, np.ndarray] = defaultdict(lambda: np.zeros(8, dtype=np.float64))
    width_histograms: dict[tuple[float, float], np.ndarray] = {}
    totals = np.zeros(6, dtype=np.float64)
    rows = 0
    parse_errors = 0
    invalid_spatial_rows = 0
    negative_primary_rows = 0
    sentinel_counts = {"primary": 0, "secondary": 0, "x": 0, "y": 0}
    first_time: datetime | None = None
    last_time: datetime | None = None
    acquisition_period_us: int | None = None
    in_stream = False
    stream_section_count = 0

    with record_path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if not row:
                continue
            first = row[0].strip()
            if first.startswith("acquisition period"):
                try:
                    acquisition_period_us = int(float(first.split("=", 1)[1].strip()))
                except (ValueError, IndexError):
                    pass
                continue
            if first == "SUBMAP_NUMBER":
                stream_section_count += 1
                if stream_section_count > 1:
                    # A canceled layer can be serialized twice in one raw file:
                    # an initial acquisition block, event text, and then a newer,
                    # slightly more complete copy with the same timestamps.  A
                    # map_record file represents one layer part, so retain only
                    # the last numerical section instead of counting both copies.
                    second_values = defaultdict(
                        lambda: np.zeros(8, dtype=np.float64)
                    )
                    width_histograms = {}
                    totals = np.zeros(6, dtype=np.float64)
                    rows = 0
                    parse_errors = 0
                    invalid_spatial_rows = 0
                    negative_primary_rows = 0
                    sentinel_counts = {
                        "primary": 0,
                        "secondary": 0,
                        "x": 0,
                        "y": 0,
                    }
                    first_time = None
                    last_time = None
                in_stream = True
                continue
            if not in_stream or len(row) != 41:
                continue
            try:
                dt, epoch_second = parse_iso_second(row[1].strip())
            except (ValueError, IndexError):
                # Event text after the numerical section is intentionally ignored.
                continue
            values = [safe_float(row, idx) for idx in (6, 7, 26, 27, 4, 5, 2, 3)]
            if any(value is None for value in values):
                parse_errors += 1
                continue
            x_dose, y_dose, primary, secondary, x_pos, y_pos, x_width, y_width = (
                float(value) for value in values
            )
            rows += 1
            first_time = dt if first_time is None or dt < first_time else first_time
            last_time = dt if last_time is None or dt > last_time else last_time
            # The recorder uses -10000 not only for missing positions/widths but
            # occasionally for individual dose channels in otherwise valid rows.
            # A missing primary channel makes the sample unusable for canonical
            # charge, time, and spatial products.
            if primary <= -9999.0:
                sentinel_counts["primary"] += 1
                continue
            if secondary <= -9999.0:
                sentinel_counts["secondary"] += 1
                secondary = 0.0
            if x_dose <= -9999.0:
                sentinel_counts["x"] += 1
                x_dose = 0.0
            if y_dose <= -9999.0:
                sentinel_counts["y"] += 1
                y_dose = 0.0
            primary_positive = max(primary, 0.0)
            if primary < 0:
                negative_primary_rows += 1

            # [raw primary, raw secondary, raw x, raw y, positive primary,
            #  mapped positive primary, row count, invalid-spatial positive primary]
            bucket = second_values[epoch_second]
            bucket[0] += primary
            bucket[1] += secondary
            bucket[2] += x_dose
            bucket[3] += y_dose
            bucket[4] += primary_positive
            bucket[6] += 1.0
            totals[:5] += [primary, secondary, x_dose, y_dose, primary_positive]

            valid_spatial = (
                x_pos > -9999.0
                and y_pos > -9999.0
                and x_width > 0.0
                and y_width > 0.0
                and grid_min <= x_pos < grid_max
                and grid_min <= y_pos < grid_max
            )
            if primary_positive > 0 and valid_spatial:
                x_idx = int((x_pos - grid_min) // pixel_mm)
                y_idx = int((y_pos - grid_min) // pixel_mm)
                xq = round(x_width / quantization) * quantization
                yq = round(y_width / quantization) * quantization
                key = (max(xq, quantization), max(yq, quantization))
                histogram = width_histograms.get(key)
                if histogram is None:
                    histogram = np.zeros((n_pixels, n_pixels), dtype=np.float64)
                    width_histograms[key] = histogram
                histogram[y_idx, x_idx] += primary_positive
                bucket[5] += primary_positive
                totals[5] += primary_positive
            elif primary_positive > 0:
                invalid_spatial_rows += 1
                bucket[7] += primary_positive

    heatmap = finalize_heatmap(
        width_histograms, pixel_mm, str(grid["width_interpretation"])
    )
    seconds = np.asarray(sorted(second_values), dtype=np.int64)
    series = (
        np.vstack([second_values[int(second)] for second in seconds])
        if len(seconds)
        else np.empty((0, 8), dtype=np.float64)
    )
    heatmap_sum = float(heatmap.sum())
    mapped_sum = float(totals[5])
    conservation_error = (
        (heatmap_sum - mapped_sum) / mapped_sum if mapped_sum > 0 else 0.0
    )
    summary = {
        **entry,
        "parsed_at": datetime.now(timezone.utc).isoformat(),
        "acquisition_period_us": acquisition_period_us,
        "stream_section_count": stream_section_count,
        "ignored_repeated_stream_sections": max(stream_section_count - 1, 0),
        "acquisition_rows": rows,
        "parse_errors": parse_errors,
        "sentinel_charge_counts": sentinel_counts,
        "negative_primary_rows": negative_primary_rows,
        "invalid_spatial_rows": invalid_spatial_rows,
        "first_time": first_time.isoformat() if first_time else None,
        "last_time": last_time.isoformat() if last_time else None,
        "charge_primary_raw_c": float(totals[0]),
        "charge_secondary_raw_c": float(totals[1]),
        "charge_x_raw_c": float(totals[2]),
        "charge_y_raw_c": float(totals[3]),
        "charge_primary_positive_c": float(totals[4]),
        "charge_spatialized_positive_c": mapped_sum,
        "charge_heatmap_c": heatmap_sum,
        "heatmap_conservation_relative_error": conservation_error,
        "width_bins": len(width_histograms),
    }
    arrays = {
        "seconds": seconds,
        "series": series,
        "heatmap_charge_c": heatmap.astype(np.float32),
    }
    return summary, arrays


def selected_entries(manifest: dict[str, Any], args: argparse.Namespace) -> list[dict[str, Any]]:
    entries = list(manifest["entries"])
    if args.session_filter:
        entries = [item for item in entries if args.session_filter in item["session"]]
    if args.record_filter:
        entries = [item for item in entries if args.record_filter in item["id"]]
    if args.task_id is not None:
        entries = [item for item in entries if item["task"] == args.task_id]
    if args.limit is not None:
        entries = entries[: args.limit]
    return entries


def parse_files(args: argparse.Namespace) -> None:
    manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    entries = selected_entries(manifest, args)
    if not entries:
        raise ValueError("No manifest entries match the requested parse selection")
    print(f"Parsing {len(entries)} record file(s)", flush=True)
    failures: list[dict[str, str]] = []
    for number, entry in enumerate(entries, 1):
        npz_path = output_dir / f"{entry['id']}.npz"
        json_path = output_dir / f"{entry['id']}.json"
        if args.resume and npz_path.exists() and json_path.exists():
            print(f"[{number}/{len(entries)}] resume skip {entry['id']}", flush=True)
            continue
        print(
            f"[{number}/{len(entries)}] {entry['id']} "
            f"({entry['size_bytes'] / 1024**2:.1f} MiB)",
            flush=True,
        )
        try:
            summary, arrays = parse_record(entry, manifest["grid"])
            tmp_npz = npz_path.with_suffix(".npz.tmp")
            with tmp_npz.open("wb") as handle:
                np.savez_compressed(handle, **arrays)
            tmp_npz.replace(npz_path)
            json_dump(json_path, summary)
            print(
                f"  rows={summary['acquisition_rows']} "
                f"Qprim={summary['charge_primary_raw_c']:.8g} C "
                f"mapped={summary['charge_spatialized_positive_c']:.8g} C "
                f"conservation={summary['heatmap_conservation_relative_error']:.3g}",
                flush=True,
            )
        except Exception as exc:  # continue array work but preserve a machine-readable failure
            failures.append({"id": entry["id"], "error": repr(exc)})
            print(f"ERROR {entry['id']}: {exc!r}", file=sys.stderr, flush=True)
    status_path = output_dir / f"_task_{args.task_id if args.task_id is not None else 'selection'}.json"
    json_dump(
        status_path,
        {
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "requested": len(entries),
            "failures": failures,
        },
    )
    if failures:
        raise RuntimeError(f"{len(failures)} file(s) failed; see {status_path}")


def load_parsed(output_dir: Path, entry: dict[str, Any]) -> tuple[dict[str, Any], dict[str, np.ndarray]] | None:
    json_path = output_dir / f"{entry['id']}.json"
    npz_path = output_dir / f"{entry['id']}.npz"
    if not json_path.exists() or not npz_path.exists():
        return None
    summary = json.loads(json_path.read_text(encoding="utf-8"))
    with np.load(npz_path) as data:
        arrays = {name: data[name].copy() for name in data.files}
    return summary, arrays


def save_heatmap_png(
    path: Path,
    values: np.ndarray,
    grid: dict[str, Any],
    title: str,
    colorbar_label: str,
    x_label: str = "Machine X [mm]",
    y_label: str = "Machine Y [mm]",
    centered_rectangle_mm: tuple[float, float] | None = None,
    rectangle_label: str | None = None,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    fig, ax = plt.subplots(figsize=(8.4, 7.2), constrained_layout=True)
    x_min, x_max, y_min, y_max, _pixel_x, _pixel_y = grid_geometry(grid)
    extent = [
        x_min,
        x_max,
        y_min,
        y_max,
    ]
    image = ax.imshow(values, origin="lower", extent=extent, cmap="magma")
    x_limits, y_limits = occupied_axis_limits(values, grid)
    if centered_rectangle_mm is not None:
        rectangle_x_mm, rectangle_y_mm = centered_rectangle_mm
        ax.add_patch(
            Rectangle(
                (-rectangle_x_mm / 2.0, -rectangle_y_mm / 2.0),
                rectangle_x_mm,
                rectangle_y_mm,
                fill=False,
                edgecolor="limegreen",
                linewidth=1.25,
                label=rectangle_label,
            )
        )
        outline_margin_mm = 5.0
        x_limits = (
            min(x_limits[0], -rectangle_x_mm / 2.0 - outline_margin_mm),
            max(x_limits[1], rectangle_x_mm / 2.0 + outline_margin_mm),
        )
        y_limits = (
            min(y_limits[0], -rectangle_y_mm / 2.0 - outline_margin_mm),
            max(y_limits[1], rectangle_y_mm / 2.0 + outline_margin_mm),
        )
        if rectangle_label:
            ax.legend(loc="upper right")
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    fig.colorbar(image, ax=ax, label=colorbar_label)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def occupied_axis_limits(
    values: np.ndarray,
    grid: dict[str, Any],
    relative_threshold: float = 1e-6,
    padding_mm: float = 5.0,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Return plot bounds around non-negligible fluence without cropping data.

    The Gaussian reconstruction is truncated at five sigma.  A small relative
    threshold ignores only the faint numerical edge of that kernel, and a fixed
    physical margin keeps the outer beam penumbra visible.
    """
    finite = np.where(np.isfinite(values), values, 0.0)
    peak = float(np.max(finite)) if finite.size else 0.0
    x_min, x_max, y_min, y_max, pixel_x, pixel_y = grid_geometry(grid)
    if peak <= 0:
        return (x_min, x_max), (y_min, y_max)
    rows, columns = np.where(finite > peak * relative_threshold)
    if not len(rows):
        return (x_min, x_max), (y_min, y_max)
    x_low = x_min + float(columns.min()) * pixel_x - padding_mm
    x_high = x_min + float(columns.max() + 1) * pixel_x + padding_mm
    y_low = y_min + float(rows.min()) * pixel_y - padding_mm
    y_high = y_min + float(rows.max() + 1) * pixel_y + padding_mm
    return (
        (max(x_min, x_low), min(x_max, x_high)),
        (max(y_min, y_low), min(y_max, y_high)),
    )


def grid_geometry(
    grid: dict[str, Any],
) -> tuple[float, float, float, float, float, float]:
    minimum = float(grid["minimum_mm"])
    maximum = float(grid["maximum_mm"])
    pixel = float(grid["pixel_mm"])
    return (
        float(grid.get("x_minimum_mm", minimum)),
        float(grid.get("x_maximum_mm", maximum)),
        float(grid.get("y_minimum_mm", minimum)),
        float(grid.get("y_maximum_mm", maximum)),
        float(grid.get("pixel_x_mm", pixel)),
        float(grid.get("pixel_y_mm", pixel)),
    )


def build_surface_time_rows(
    experiment_products,
    group_ids,
    surface_name,
    area_by_group_id,
):
    prefix = surface_name.lower()
    area_field = (
        "radex_instrument_area_cm2"
        if prefix == "radex"
        else f"{prefix}_target_area_cm2"
    )
    result = defaultdict(lambda: defaultdict(float))
    for group_id in group_ids:
        area_cm2 = float(area_by_group_id[group_id])
        if area_cm2 <= 0.0:
            raise ValueError(f"Invalid {surface_name} area for {group_id}: {area_cm2}")
        for second, values in experiment_products[group_id]["time_rows"].items():
            target = result[int(second)]
            for field in ("charge_primary_raw_c", "charge_primary_positive_c", "charge_secondary_raw_c",
                          "protons_delivered_raw", "protons_delivered_positive", "tuning_protons_estimated"):
                target[field] += float(values.get(field, 0.0))
            for source, suffix in (
                ("protons_delivered_raw", "mean_fluence_raw_protons_cm2"),
                ("protons_delivered_positive", "mean_fluence_positive_protons_cm2"),
                ("tuning_protons_estimated", "mean_tuning_fluence_protons_cm2"),
            ):
                target[f"{prefix}_{suffix}"] += (
                    float(values.get(source, 0.0)) / area_cm2
                )
            existing_area = target.get(area_field)
            if existing_area not in (None, 0.0) and not math.isclose(
                existing_area, area_cm2
            ):
                raise ValueError(
                    f"Overlapping {surface_name} target areas at epoch {second}"
                )
            target[area_field] = area_cm2
    return result


def build_radex_time_rows(experiment_products, group_ids):
    return build_surface_time_rows(
        experiment_products,
        group_ids,
        "RadEx",
        {group_id: RADEX_INSTRUMENT_AREA_CM2 for group_id in group_ids},
    )


def build_crds_time_rows(experiment_products, group_ids):
    return build_surface_time_rows(
        experiment_products,
        group_ids,
        "CRDS",
        {
            group_id: CRDS_TARGETS[group_id]["area_cm2"]
            for group_id in group_ids
        },
    )


def write_time_series_csv(
    path: Path,
    time_rows: dict[int, dict[str, float]],
    bms_timezone: timezone,
) -> None:
    cumulative_raw = 0.0
    cumulative_positive = 0.0
    surface_prefixes = [
        prefix
        for prefix in ("radex", "crds")
        if any(
            f"{prefix}_mean_fluence_raw_protons_cm2" in row
            for row in time_rows.values()
        )
    ]
    if len(surface_prefixes) > 1:
        raise ValueError("A time series may contain only one surface-fluence definition")
    surface_prefix = surface_prefixes[0] if surface_prefixes else None
    cumulative_surface_raw = 0.0
    cumulative_surface_positive = 0.0
    fields = [
        "time_utc",
        "time_bms_local",
        "epoch_second",
        "charge_primary_raw_c",
        "charge_primary_positive_c",
        "charge_secondary_raw_c",
        "protons_delivered_raw",
        "cumulative_protons_delivered_raw",
        "protons_delivered_positive",
        "cumulative_protons_delivered_positive",
        "tuning_protons_estimated",
    ]
    if surface_prefix:
        area_field = (
            "radex_instrument_area_cm2"
            if surface_prefix == "radex"
            else f"{surface_prefix}_target_area_cm2"
        )
        fields.extend(
            [
                area_field,
                f"{surface_prefix}_mean_fluence_raw_protons_cm2",
                f"cumulative_{surface_prefix}_mean_fluence_raw_protons_cm2",
                f"{surface_prefix}_mean_fluence_positive_protons_cm2",
                f"cumulative_{surface_prefix}_mean_fluence_positive_protons_cm2",
                f"{surface_prefix}_mean_tuning_fluence_protons_cm2",
            ]
        )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for second in sorted(time_rows):
            row = time_rows[second]
            cumulative_raw += row.get("protons_delivered_raw", 0.0)
            cumulative_positive += row.get("protons_delivered_positive", 0.0)
            if surface_prefix:
                cumulative_surface_raw += row.get(
                    f"{surface_prefix}_mean_fluence_raw_protons_cm2", 0.0
                )
                cumulative_surface_positive += row.get(
                    f"{surface_prefix}_mean_fluence_positive_protons_cm2", 0.0
                )
            utc_time = datetime.fromtimestamp(second, timezone.utc)
            output = {
                "time_utc": utc_time.isoformat(),
                "time_bms_local": utc_time.astimezone(bms_timezone).isoformat(),
                "epoch_second": second,
                "charge_primary_raw_c": row.get("charge_primary_raw_c", 0.0),
                "charge_primary_positive_c": row.get(
                    "charge_primary_positive_c", 0.0
                ),
                "charge_secondary_raw_c": row.get(
                    "charge_secondary_raw_c", 0.0
                ),
                "protons_delivered_raw": row.get("protons_delivered_raw", 0.0),
                "cumulative_protons_delivered_raw": cumulative_raw,
                "protons_delivered_positive": row.get(
                    "protons_delivered_positive", 0.0
                ),
                "cumulative_protons_delivered_positive": cumulative_positive,
                "tuning_protons_estimated": row.get(
                    "tuning_protons_estimated", 0.0
                ),
            }
            if surface_prefix:
                output.update(
                    {
                        area_field: row.get(area_field, 0.0),
                        f"{surface_prefix}_mean_fluence_raw_protons_cm2": row.get(
                            f"{surface_prefix}_mean_fluence_raw_protons_cm2", 0.0
                        ),
                        f"cumulative_{surface_prefix}_mean_fluence_raw_protons_cm2": cumulative_surface_raw,
                        f"{surface_prefix}_mean_fluence_positive_protons_cm2": row.get(
                            f"{surface_prefix}_mean_fluence_positive_protons_cm2", 0.0
                        ),
                        f"cumulative_{surface_prefix}_mean_fluence_positive_protons_cm2": cumulative_surface_positive,
                        f"{surface_prefix}_mean_tuning_fluence_protons_cm2": row.get(
                            f"{surface_prefix}_mean_tuning_fluence_protons_cm2", 0.0
                        ),
                    }
                )
            writer.writerow(output)


def save_time_series_png(
    path: Path,
    time_rows: dict[int, dict[str, float]],
    bms_timezone: timezone,
    title: str,
    surface_name: str | None = None,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    if not time_rows:
        raise ValueError(f"No time-series data available for {title}")
    first_second = min(time_rows)
    last_second = max(time_rows)
    seconds = np.arange(first_second, last_second + 1, dtype=np.int64)
    delivered = np.zeros(len(seconds), dtype=np.float64)
    tuning = np.zeros(len(seconds), dtype=np.float64)
    for second, values in time_rows.items():
        index = int(second - first_second)
        if surface_name:
            prefix = surface_name.lower()
            delivered[index] = values.get(
                f"{prefix}_mean_fluence_raw_protons_cm2", 0.0
            )
            tuning[index] = values.get(
                f"{prefix}_mean_tuning_fluence_protons_cm2", 0.0
            )
        else:
            delivered[index] = values.get("protons_delivered_raw", 0.0)
            tuning[index] = values.get("tuning_protons_estimated", 0.0)
    cumulative = np.cumsum(delivered)
    times = [
        datetime.fromtimestamp(int(second), timezone.utc).astimezone(bms_timezone)
        for second in seconds
    ]

    fig, (rate_ax, cumulative_ax) = plt.subplots(
        2, 1, figsize=(11.0, 7.5), sharex=True, constrained_layout=True
    )
    rate_ax.step(times, delivered, where="mid", linewidth=0.9, color="#3366cc")
    if np.max(tuning) >= np.max(delivered) * 1e-3:
        rate_ax.step(
            times,
            tuning,
            where="mid",
            linewidth=0.8,
            color="#dd7711",
            alpha=0.8,
            label="Estimated tuning beam",
        )
        rate_ax.legend(loc="upper right")
    rate_ax.set_ylabel(
        f"Mean {surface_name} fluence rate [protons cm$^{{-2}}$ s$^{{-1}}$]"
        if surface_name
        else "Delivered protons s$^{-1}$"
    )
    rate_ax.set_title(title)
    rate_ax.grid(alpha=0.25)
    rate_ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    cumulative_ax.plot(times, cumulative, linewidth=1.5, color="#228833")
    cumulative_ax.set_ylabel(
        f"Cumulative mean {surface_name} fluence [protons cm$^{{-2}}$]"
        if surface_name
        else "Cumulative delivered protons"
    )
    cumulative_ax.set_xlabel("BMS local time (+02:00)")
    cumulative_ax.grid(alpha=0.25)
    cumulative_ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    cumulative_ax.xaxis.set_major_formatter(
        mdates.DateFormatter("%H:%M", tz=bms_timezone)
    )
    cumulative_ax.text(
        0.99,
        0.05,
        (
            f"Total: {cumulative[-1]:.4e} protons cm$^{{-2}}$"
            if surface_name
            else f"Total: {cumulative[-1]:.4e} protons"
        ),
        transform=cumulative_ax.transAxes,
        ha="right",
        va="bottom",
    )
    fig.savefig(path, dpi=160)
    plt.close(fig)


def monitor_comparison_arrays(
    time_rows: dict[int, dict[str, float]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if not time_rows:
        raise ValueError("No monitor time-series data available")
    first_second = min(time_rows)
    last_second = max(time_rows)
    seconds = np.arange(first_second, last_second + 1, dtype=np.int64)
    primary = np.zeros(len(seconds), dtype=np.float64)
    secondary = np.zeros(len(seconds), dtype=np.float64)
    for second, values in time_rows.items():
        index = int(second - first_second)
        primary[index] = float(values.get("charge_primary_raw_c", 0.0))
        secondary[index] = float(values.get("charge_secondary_raw_c", 0.0))
    difference = secondary - primary
    relative = np.full(len(seconds), np.nan, dtype=np.float64)
    threshold = max(float(np.max(np.abs(primary))) * 1e-3, 1e-15)
    active = primary > threshold
    relative[active] = 100.0 * difference[active] / primary[active]
    return seconds, primary, secondary, relative


def write_monitor_comparison_csv(
    path: Path,
    time_rows: dict[int, dict[str, float]],
    bms_timezone: timezone,
) -> None:
    seconds, primary, secondary, relative = monitor_comparison_arrays(time_rows)
    fields = [
        "time_utc",
        "time_bms_local",
        "epoch_second",
        "dose_prim_c_per_second",
        "dose_sec_c_per_second",
        "dose_sec_minus_prim_c_per_second",
        "dose_sec_minus_prim_percent_of_prim",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for index, second in enumerate(seconds):
            utc_time = datetime.fromtimestamp(int(second), timezone.utc)
            writer.writerow(
                {
                    "time_utc": utc_time.isoformat(),
                    "time_bms_local": utc_time.astimezone(bms_timezone).isoformat(),
                    "epoch_second": int(second),
                    "dose_prim_c_per_second": primary[index],
                    "dose_sec_c_per_second": secondary[index],
                    "dose_sec_minus_prim_c_per_second": (
                        secondary[index] - primary[index]
                    ),
                    "dose_sec_minus_prim_percent_of_prim": (
                        relative[index] if np.isfinite(relative[index]) else ""
                    ),
                }
            )


def save_monitor_comparison_png(
    path: Path,
    time_rows: dict[int, dict[str, float]],
    bms_timezone: timezone,
    title: str,
) -> dict[str, float]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    seconds, primary, secondary, relative = monitor_comparison_arrays(time_rows)
    times = [
        datetime.fromtimestamp(int(second), timezone.utc).astimezone(bms_timezone)
        for second in seconds
    ]
    finite = np.isfinite(relative)
    total_primary = float(primary.sum())
    total_secondary = float(secondary.sum())
    integrated_difference_percent = (
        100.0 * (total_secondary - total_primary) / total_primary
        if total_primary != 0.0
        else float("nan")
    )
    median_difference_percent = (
        float(np.median(relative[finite])) if np.any(finite) else float("nan")
    )
    percentile_05 = (
        float(np.percentile(relative[finite], 5)) if np.any(finite) else float("nan")
    )
    percentile_95 = (
        float(np.percentile(relative[finite], 95)) if np.any(finite) else float("nan")
    )

    fig, (charge_ax, difference_ax) = plt.subplots(
        2, 1, figsize=(11.0, 7.5), sharex=True, constrained_layout=True
    )
    charge_ax.step(
        times, primary, where="mid", linewidth=0.9, label="DOSE_PRIM"
    )
    charge_ax.step(
        times, secondary, where="mid", linewidth=0.9, label="DOSE_SEC"
    )
    charge_ax.set_ylabel("Monitor charge per second [C]")
    charge_ax.set_title(title)
    charge_ax.grid(alpha=0.25)
    charge_ax.legend(loc="upper right")
    charge_ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    difference_ax.plot(
        times, relative, linewidth=0.8, color="#aa3377"
    )
    difference_ax.axhline(0.0, color="black", linewidth=0.7)
    difference_ax.set_ylabel("(SEC - PRIM) / PRIM [%]")
    difference_ax.set_xlabel("BMS local time (+02:00)")
    difference_ax.grid(alpha=0.25)
    difference_ax.xaxis.set_major_formatter(
        mdates.DateFormatter("%H:%M", tz=bms_timezone)
    )
    difference_ax.text(
        0.99,
        0.05,
        (
            f"Integrated difference: {integrated_difference_percent:+.3f}%\n"
            f"Active-second median: {median_difference_percent:+.3f}%\n"
            f"5th-95th percentile: {percentile_05:+.3f}% to {percentile_95:+.3f}%"
        ),
        transform=difference_ax.transAxes,
        ha="right",
        va="bottom",
    )
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return {
        "dose_prim_integrated_charge_c": total_primary,
        "dose_sec_integrated_charge_c": total_secondary,
        "dose_sec_minus_prim_integrated_percent": integrated_difference_percent,
        "active_second_median_difference_percent": median_difference_percent,
        "active_second_difference_percentile_05": percentile_05,
        "active_second_difference_percentile_95": percentile_95,
    }


def reduce_results(args: argparse.Namespace) -> None:
    manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    parsed_dir = Path(args.parsed_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    experiment_lookup = {item["id"]: item for item in manifest["experiments"]}
    parsed: dict[str, tuple[dict[str, Any], dict[str, np.ndarray]]] = {}
    missing: list[str] = []
    for entry in manifest["entries"]:
        item = load_parsed(parsed_dir, entry)
        if item is None:
            missing.append(entry["id"])
        else:
            parsed[entry["id"]] = item
    if missing and not args.allow_partial:
        raise RuntimeError(f"{len(missing)} parsed files are missing")

    grid = manifest["grid"]
    pixel_area_cm2 = (float(grid["pixel_mm"]) / 10.0) ** 2
    shape = next((arrays["heatmap_charge_c"].shape for _, arrays in parsed.values()), None)
    if shape is None:
        raise RuntimeError("No parsed products found")
    charge_map_all = np.zeros(shape, dtype=np.float64)
    charge_map_documented_part = np.zeros(shape, dtype=np.float64)
    fluence_map = np.zeros(shape, dtype=np.float64)
    time_rows: dict[int, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    calibration_rows: list[dict[str, Any]] = []
    experiment_products: dict[str, dict[str, Any]] = {}

    grouped: dict[
        str, list[tuple[dict[str, Any], dict[str, Any], dict[str, np.ndarray]]]
    ] = defaultdict(list)
    for entry in manifest["entries"]:
        item = parsed.get(entry["id"])
        if item is None:
            continue
        summary, arrays = item
        charge_map_all += arrays["heatmap_charge_c"]
        group_id = entry["experiment_id"] or f"undocumented:{entry['session']}"
        grouped[group_id].append((entry, summary, arrays))

    for group_id, items in grouped.items():
        experiment = experiment_lookup.get(group_id)
        part_items = [item for item in items if item[0]["kind"] == "part"]
        part_charge_raw = sum(item[1]["charge_primary_raw_c"] for item in part_items)
        part_charge_positive = sum(
            item[1]["charge_primary_positive_c"] for item in part_items
        )
        planned_charge = sum(
            float(item[0]["planned_charge_c"] or 0.0) for item in part_items
        )
        valid_k_factors = all(
            item[0].get("k_factor") is not None
            and float(item[0]["k_factor"]) > 0.0
            for item in part_items
        )
        part_charge_map = sum(
            (item[2]["heatmap_charge_c"].astype(np.float64) for item in part_items),
            start=np.zeros(shape, dtype=np.float64),
        )
        expected_part_ids = {
            entry["id"]
            for entry in manifest["entries"]
            if entry["experiment_id"] == group_id and entry["kind"] == "part"
        }
        complete = expected_part_ids.issubset(parsed)
        if experiment is None or not part_items or not complete or not valid_k_factors:
            calibration_rows.append(
                {
                    "group_id": group_id,
                    "calibrated": False,
                    "complete": complete,
                    "parsed_files": len(items),
                    "planned_charge_c": planned_charge,
                    "part_raw_charge_c": part_charge_raw,
                    "part_positive_charge_c": part_charge_positive,
                    "calibration_error": (
                        "Missing or non-positive K_FACTOR"
                        if part_items and not valid_k_factors
                        else None
                    ),
                }
            )
            continue

        planned_protons = float(experiment["calibration_value"])
        prescribed_mu = float(experiment["prescribed_mu"])
        protons_per_mu = float(experiment["protons_per_mu"])
        delivered_mu_raw = 0.0
        delivered_protons_raw = 0.0
        delivered_protons_positive = 0.0
        spatialized_protons = 0.0
        group_fluence = np.zeros(shape, dtype=np.float64)
        k_factors: list[float] = []
        for entry, summary, arrays in part_items:
            k_factor = float(entry["k_factor"])
            k_factors.append(k_factor)
            protons_per_c = monitor_protons_per_c(protons_per_mu, k_factor)
            delivered_mu_raw += monitor_mu_from_charge(
                float(summary["charge_primary_raw_c"]), k_factor
            )
            delivered_protons_raw += (
                float(summary["charge_primary_raw_c"]) * protons_per_c
            )
            delivered_protons_positive += (
                float(summary["charge_primary_positive_c"]) * protons_per_c
            )
            entry_heatmap = arrays["heatmap_charge_c"].astype(np.float64)
            entry_spatialized_protons = float(entry_heatmap.sum()) * protons_per_c
            spatialized_protons += entry_spatialized_protons
            group_fluence += entry_heatmap * protons_per_c / pixel_area_cm2
        delivered_fraction = delivered_mu_raw / prescribed_mu
        spatialized_charge = float(part_charge_map.sum())
        calibration_summary = {
            "group_id": group_id,
            "calibrated": True,
            "complete": True,
            "section": experiment["section"],
            "label": experiment["label"],
            "energy_mev": experiment["energy_mev"],
            "calibration_energy_mev": experiment["calibration_energy_mev"],
            "prescribed_mu": prescribed_mu,
            "protons_per_mu": protons_per_mu,
            "planned_protons": planned_protons,
            "planned_charge_c": planned_charge,
            "part_raw_charge_c": part_charge_raw,
            "part_positive_charge_c": part_charge_positive,
            "spatialized_positive_charge_c": spatialized_charge,
            "delivered_to_prescribed_mu_ratio": delivered_fraction,
            "delivered_mu_raw": delivered_mu_raw,
            "delivered_protons_raw": delivered_protons_raw,
            "delivered_protons_positive": delivered_protons_positive,
            "spatialized_protons": spatialized_protons,
            "monitor_charge_c_per_mu": MONITOR_CHARGE_C_PER_MU,
            "k_factor_min": min(k_factors),
            "k_factor_max": max(k_factors),
            "parsed_files": len(items),
        }
        calibration_rows.append(calibration_summary)
        charge_map_documented_part += part_charge_map
        fluence_map += group_fluence
        group_time_rows: dict[int, dict[str, float]] = defaultdict(
            lambda: defaultdict(float)
        )
        experiment_products[group_id] = {
            "experiment": experiment,
            "charge_part_c": part_charge_map,
            "fluence": group_fluence,
            "calibration": calibration_summary,
            "time_rows": group_time_rows,
        }
        for entry, _summary, arrays in items:
            k_factor = float(entry["k_factor"] or 0.0)
            if k_factor <= 0.0:
                continue
            protons_per_c = monitor_protons_per_c(protons_per_mu, k_factor)
            series = arrays["series"]
            seconds = arrays["seconds"]
            for idx, second in enumerate(seconds):
                q_raw = float(series[idx, 0])
                q_positive = float(series[idx, 4])
                row = time_rows[int(second)]
                row["charge_primary_raw_c"] += q_raw
                row["charge_primary_positive_c"] += q_positive
                row["charge_secondary_raw_c"] += float(series[idx, 1])
                group_row = group_time_rows[int(second)]
                group_row["charge_primary_raw_c"] += q_raw
                group_row["charge_primary_positive_c"] += q_positive
                group_row["charge_secondary_raw_c"] += float(series[idx, 1])
                if entry["kind"] == "part":
                    row["protons_delivered_raw"] += q_raw * protons_per_c
                    row["protons_delivered_positive"] += q_positive * protons_per_c
                    group_row["protons_delivered_raw"] += q_raw * protons_per_c
                    group_row["protons_delivered_positive"] += (
                        q_positive * protons_per_c
                    )
                else:
                    row["tuning_protons_estimated"] += q_positive * protons_per_c
                    group_row["tuning_protons_estimated"] += (
                        q_positive * protons_per_c
                    )

    radex_group_ids = sorted(
        group_id
        for group_id, product in experiment_products.items()
        if product["experiment"]["section"] == "Experiment 1"
    )
    radex_campaign_time_rows = build_radex_time_rows(
        experiment_products, radex_group_ids
    )
    crds_group_ids = sorted(
        group_id for group_id in CRDS_TARGETS if group_id in experiment_products
    )
    if set(crds_group_ids) != set(CRDS_TARGETS):
        missing_crds = sorted(set(CRDS_TARGETS) - set(crds_group_ids))
        raise RuntimeError(f"Missing documented CRDS experiments: {missing_crds}")
    crds_campaign_time_rows = build_crds_time_rows(
        experiment_products, crds_group_ids
    )
    radex_area_report = {"width_x_mm": RADEX_INSTRUMENT_WIDTH_X_MM,
        "width_y_mm": RADEX_INSTRUMENT_WIDTH_Y_MM, "area_cm2": RADEX_INSTRUMENT_AREA_CM2,
        "uniformity_assurance_percent": 5.0,
        "method": ("Monitor-derived delivered protons divided by the explicitly requested "
        "154.0 mm by 66.8 mm RadEx surface. The facility assured fluence uniformity "
        "within 5%. No X/Y magnification or exp_12/exp_13 fit is applied.")}
    crds_area_report = {
        "method": (
            "Monitor-derived delivered protons divided by the requested CRDS target "
            "area for each irradiation: 5.2 cm by 5.4 cm for exp_10 through exp_12, "
            "and 7.0 cm by 7.4 cm for exp_13."
        ),
        "experiments": {
            group_id: dict(CRDS_TARGETS[group_id]) for group_id in crds_group_ids
        },
    }
    for group_id in radex_group_ids:
        product = experiment_products[group_id]
        c = product["calibration"]
        c.update({"radex_instrument_width_x_mm": RADEX_INSTRUMENT_WIDTH_X_MM,
        "radex_instrument_width_y_mm": RADEX_INSTRUMENT_WIDTH_Y_MM,
        "radex_instrument_area_cm2": RADEX_INSTRUMENT_AREA_CM2,
        "radex_mean_fluence_raw_protons_cm2": c["delivered_protons_raw"]/RADEX_INSTRUMENT_AREA_CM2,
        "radex_mean_fluence_positive_protons_cm2": c["delivered_protons_positive"]/RADEX_INSTRUMENT_AREA_CM2})
    for group_id in crds_group_ids:
        product = experiment_products[group_id]
        c = product["calibration"]
        target = CRDS_TARGETS[group_id]
        c.update({
            "crds_target_width_x_cm": target["width_x_cm"],
            "crds_target_width_y_cm": target["width_y_cm"],
            "crds_target_area_cm2": target["area_cm2"],
            "crds_mean_fluence_raw_protons_cm2": (
                c["delivered_protons_raw"] / target["area_cm2"]
            ),
            "crds_mean_fluence_positive_protons_cm2": (
                c["delivered_protons_positive"] / target["area_cm2"]
            ),
        })

    np.savez_compressed(
        output_dir / "campaign_heatmaps.npz",
        charge_all_c=charge_map_all.astype(np.float32),
        charge_documented_part_c=charge_map_documented_part.astype(np.float32),
        fluence_protons_cm2=fluence_map.astype(np.float32),
        x_edges_mm=np.linspace(
            grid["minimum_mm"], grid["maximum_mm"], shape[1] + 1, dtype=np.float64
        ),
        y_edges_mm=np.linspace(
            grid["minimum_mm"], grid["maximum_mm"], shape[0] + 1, dtype=np.float64
        ),
    )
    save_heatmap_png(
        output_dir / "fluence_protons_cm2.png",
        fluence_map,
        grid,
        "Campaign delivered fluence",
        "Protons cm$^{-2}$",
    )
    experiment_dir = output_dir / "experiments"
    experiment_dir.mkdir(parents=True, exist_ok=True)
    for group_id, product in experiment_products.items():
        experiment = product["experiment"]
        np.savez_compressed(
            experiment_dir / f"{group_id}.npz",
            charge_part_c=product["charge_part_c"].astype(np.float32),
            fluence_protons_cm2=product["fluence"].astype(np.float32),
        )
        description = (
            f"{experiment['section']} — {experiment['label']}, "
            f"{experiment['energy_mev']:g} MeV, {experiment['prescribed_mu']:g} MU"
        )
        save_heatmap_png(
            experiment_dir / f"{group_id}_fluence.png",
            product["fluence"],
            grid,
            description + "\nDelivered fluence from primary monitor charge",
            "Protons cm$^{-2}$",
        )
    named_aggregates = {
        "RadEx_64MeV": [
            group_id
            for group_id, product in experiment_products.items()
            if product["experiment"]["section"] == "Experiment 1"
            and float(product["experiment"]["energy_mev"]) == 64.0
        ],
        "Radex_85MeV": [
            group_id
            for group_id, product in experiment_products.items()
            if product["experiment"]["section"] == "Experiment 1"
            and float(product["experiment"]["energy_mev"]) == 85.0
        ],
    }
    named_aggregates["RadEx_Total"] = (
        named_aggregates["RadEx_64MeV"] + named_aggregates["Radex_85MeV"]
    )
    named_aggregates.update(
        {
            "CRDS_5.2x5.4cm": ["exp_10", "exp_11", "exp_12"],
            "CRDS_7.0x7.4cm": ["exp_13"],
            "CRDS_Total": crds_group_ids,
        }
    )
    aggregate_root = output_dir / "aggregates"
    aggregate_root.mkdir(parents=True, exist_ok=True)
    aggregate_summaries: list[dict[str, Any]] = []
    bms_timezone = timezone(
        timedelta(hours=float(manifest["timezone_offset_hours"]))
    )
    for name, group_ids in named_aggregates.items():
        if not group_ids:
            raise RuntimeError(f"Named aggregate {name} has no experiment groups")
        aggregate_dir = aggregate_root / name
        aggregate_dir.mkdir(parents=True, exist_ok=True)
        aggregate_charge = sum(
            (
                experiment_products[group_id]["charge_part_c"]
                for group_id in group_ids
            ),
            start=np.zeros(shape, dtype=np.float64),
        )
        aggregate_fluence = sum(
            (experiment_products[group_id]["fluence"] for group_id in group_ids),
            start=np.zeros(shape, dtype=np.float64),
        )
        aggregate_time: dict[int, dict[str, float]] = defaultdict(
            lambda: defaultdict(float)
        )
        for group_id in group_ids:
            for second, values in experiment_products[group_id]["time_rows"].items():
                target = aggregate_time[second]
                for field, value in values.items():
                    target[field] += value
        surface_name = "CRDS" if name.startswith("CRDS") else "RadEx"
        aggregate_time_surface = (
            build_crds_time_rows(experiment_products, group_ids)
            if surface_name == "CRDS"
            else build_radex_time_rows(experiment_products, group_ids)
        )
        summaries = [
            experiment_products[group_id]["calibration"] for group_id in group_ids
        ]
        summary = {
            "name": name,
            "experiment_ids": group_ids,
            "irradiation_count": len(group_ids),
            "energies_mev": sorted(
                {float(item["energy_mev"]) for item in summaries}
            ),
            "prescribed_mu": sum(float(item["prescribed_mu"]) for item in summaries),
            "delivered_mu_raw": sum(
                float(item["delivered_mu_raw"]) for item in summaries
            ),
            "planned_protons": sum(
                float(item["planned_protons"]) for item in summaries
            ),
            "delivered_protons_raw": sum(
                float(item["delivered_protons_raw"]) for item in summaries
            ),
            "delivered_protons_positive": sum(
                float(item["delivered_protons_positive"]) for item in summaries
            ),
            "spatialized_protons": sum(
                float(item["spatialized_protons"]) for item in summaries
            ),
        }
        if surface_name == "RadEx":
            summary.update({
                "radex_instrument_width_x_mm": RADEX_INSTRUMENT_WIDTH_X_MM,
                "radex_instrument_width_y_mm": RADEX_INSTRUMENT_WIDTH_Y_MM,
                "radex_instrument_area_cm2": RADEX_INSTRUMENT_AREA_CM2,
            })
        else:
            dimensions = sorted(
                {
                    (
                        CRDS_TARGETS[group_id]["width_x_cm"],
                        CRDS_TARGETS[group_id]["width_y_cm"],
                    )
                    for group_id in group_ids
                }
            )
            summary.update({
                "crds_target_dimensions_cm": [list(item) for item in dimensions],
                "crds_target_areas_cm2": sorted(
                    {CRDS_TARGETS[group_id]["area_cm2"] for group_id in group_ids}
                ),
            })
        summary["nominal_comparison_valid"] = "exp_09" not in group_ids
        summary["nominal_comparison_note"] = (
            "Includes early-canceled exp_09; its 10,000 MU summary value was an "
            "estimate, so delivered/planned comparison is not meaningful."
            if not summary["nominal_comparison_valid"]
            else "All included irradiations have prescribed MU targets."
        )
        summary["delivered_to_planned_protons_ratio"] = (
            summary["delivered_protons_raw"] / summary["planned_protons"]
            if summary["nominal_comparison_valid"]
            else None
        )
        summary["spatialized_to_delivered_protons_ratio"] = (
            summary["spatialized_protons"] / summary["delivered_protons_raw"]
        )
        if surface_name == "RadEx":
            summary["radex_mean_fluence_raw_protons_cm2"] = (
                summary["delivered_protons_raw"] / RADEX_INSTRUMENT_AREA_CM2
            )
            summary["radex_mean_fluence_positive_protons_cm2"] = (
                summary["delivered_protons_positive"] / RADEX_INSTRUMENT_AREA_CM2
            )
        else:
            summary["crds_mean_fluence_raw_protons_cm2"] = sum(
                float(item["crds_mean_fluence_raw_protons_cm2"])
                for item in summaries
            )
            summary["crds_mean_fluence_positive_protons_cm2"] = sum(
                float(item["crds_mean_fluence_positive_protons_cm2"])
                for item in summaries
            )
        aggregate_summaries.append(summary)
        json_dump(aggregate_dir / "summary.json", summary)
        np.savez_compressed(
            aggregate_dir / "heatmap.npz",
            charge_part_c=aggregate_charge.astype(np.float32),
            fluence_protons_cm2=aggregate_fluence.astype(np.float32),
            x_edges_mm=np.linspace(
                grid["minimum_mm"],
                grid["maximum_mm"],
                shape[1] + 1,
                dtype=np.float64,
            ),
            y_edges_mm=np.linspace(
                grid["minimum_mm"],
                grid["maximum_mm"],
                shape[0] + 1,
                dtype=np.float64,
            ),
        )
        save_heatmap_png(
            aggregate_dir / "fluence.png",
            aggregate_fluence,
            grid,
            f"{name} — aggregated delivered fluence",
            "Protons cm$^{-2}$",
        )
        if name == "RadEx_Total":
            instrument_size_mm = (154.0, 66.8)
            instrument_label = "RadEx instrument (154 x 66.8 mm)"
            save_heatmap_png(
                aggregate_dir / "fluence_with_instrument_outline.png",
                aggregate_fluence,
                grid,
                f"{name} — aggregated delivered fluence with instrument outline",
                "Protons cm$^{-2}$",
                centered_rectangle_mm=instrument_size_mm,
                rectangle_label=instrument_label,
            )
        write_time_series_csv(
            aggregate_dir / "fluence_vs_time_1s.csv",
            aggregate_time,
            bms_timezone,
        )
        save_time_series_png(
            aggregate_dir / "fluence_vs_time.png",
            aggregate_time,
            bms_timezone,
            f"{name} — delivered protons versus time",
        )
        write_monitor_comparison_csv(
            aggregate_dir / "dose_prim_vs_sec_1s.csv",
            aggregate_time,
            bms_timezone,
        )
        monitor_summary = save_monitor_comparison_png(
            aggregate_dir / "dose_prim_vs_sec.png",
            aggregate_time,
            bms_timezone,
            f"{name} — primary versus secondary dose monitor",
        )
        summary.update(monitor_summary)
        json_dump(aggregate_dir / "summary.json", summary)
        write_time_series_csv(
            aggregate_dir / f"fluence_vs_time_1s_{surface_name}_area.csv",
            aggregate_time_surface,
            bms_timezone,
        )
        save_time_series_png(
            aggregate_dir / f"fluence_vs_time_{surface_name}_area.png",
            aggregate_time_surface,
            bms_timezone,
            f"{name} — mean {surface_name}-surface fluence versus time",
            surface_name=surface_name,
        )
    aggregate_fields = sorted({key for row in aggregate_summaries for key in row})
    with (aggregate_root / "aggregate_summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=aggregate_fields)
        writer.writeheader()
        writer.writerows(aggregate_summaries)

    calibration_fields = sorted({key for row in calibration_rows for key in row})
    with (output_dir / "calibration_qc.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=calibration_fields)
        writer.writeheader()
        writer.writerows(calibration_rows)

    write_time_series_csv(
        output_dir / "fluence_vs_time_1s.csv", time_rows, bms_timezone
    )
    write_monitor_comparison_csv(
        output_dir / "dose_prim_vs_sec_1s.csv", time_rows, bms_timezone
    )
    campaign_monitor_summary = save_monitor_comparison_png(
        output_dir / "dose_prim_vs_sec.png",
        time_rows,
        bms_timezone,
        "Campaign — primary versus secondary dose monitor",
    )
    write_time_series_csv(
        output_dir / "fluence_vs_time_1s_RadEx_area.csv",
        radex_campaign_time_rows,
        bms_timezone,
    )
    save_time_series_png(
        output_dir / "fluence_vs_time_RadEx_area.png",
        radex_campaign_time_rows,
        bms_timezone,
        "Campaign — mean RadEx-surface fluence versus time",
        surface_name="RadEx",
    )
    write_time_series_csv(
        output_dir / "fluence_vs_time_1s_CRDS_area.csv",
        crds_campaign_time_rows,
        bms_timezone,
    )
    save_time_series_png(
        output_dir / "fluence_vs_time_CRDS_area.png",
        crds_campaign_time_rows,
        bms_timezone,
        "Campaign — mean CRDS-target fluence versus time",
        surface_name="CRDS",
    )

    report = {
        "schema_version": SCHEMA_VERSION,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "parsed_files": len(parsed),
        "missing_files": missing,
        "calibrated_experiments": sum(
            row.get("calibrated", False) for row in calibration_rows
        ),
        "grid": grid,
        "pixel_area_cm2": pixel_area_cm2,
        "spatial_heatmap_note": ("Recorder-plane diagnostic heatmaps retain the configured spot-width model, "
        "but spot widths are not used for the RadEx or CRDS mean-fluence estimates."),
        "calibration_method": (
            "Each DOSE_PRIM(C) sample is summed as incremental primary-monitor charge. "
            "For every part record, MU = charge / (3.0e-9 C/MU * K_FACTOR), and the "
            "facility protons/MU value converts that MU to protons. Nominal 64 MeV "
            "irradiations use the facility's 65 MeV conversion row. Map-specification "
            "TOTAL_CHARGE is retained only for QC because canceled/restarted attempts "
            "repeat the full nominal value."
        ),
        "dose_monitor_comparison": campaign_monitor_summary,
        "radex_surface_fluence": radex_area_report,
        "radex_time_series_interpretation": ("RadEx-area time-series fluence is delivered protons "
        "per second divided by the fixed 102.872 cm2 surface."),
        "crds_surface_fluence": crds_area_report,
        "crds_time_series_interpretation": (
            "CRDS-area time-series fluence is delivered protons per second divided "
            "by 28.08 cm2 for exp_10 through exp_12 and 51.8 cm2 for exp_13. "
            "Cumulative CRDS fluence is the sum of the per-irradiation fluences."
        ),
    }
    json_dump(output_dir / "reduction_report.json", report)
    print(
        f"Reduced {len(parsed)} files; {len(missing)} missing; "
        f"{report['calibrated_experiments']} fully calibrated experiment groups"
    )


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    inventory = subparsers.add_parser("inventory", help="build a file manifest")
    inventory.add_argument("--input-dir", required=True)
    inventory.add_argument("--experiment-csv", required=True)
    inventory.add_argument("--conversion-csv", required=True)
    inventory.add_argument("--output", required=True)
    inventory.add_argument("--experiment-date", default="2025-06-28")
    inventory.add_argument("--utc-offset-hours", type=int, default=2)
    inventory.add_argument("--match-tolerance-minutes", type=float, default=5.0)
    inventory.add_argument("--tasks", type=int, default=100)
    inventory.add_argument("--grid-min-mm", type=float, default=-250.0)
    inventory.add_argument("--grid-max-mm", type=float, default=250.0)
    inventory.add_argument("--pixel-mm", type=float, default=1.0)
    inventory.add_argument("--width-quantization-mm", type=float, default=0.25)
    inventory.add_argument(
        "--width-interpretation", choices=("sigma", "fwhm"), default="sigma"
    )
    inventory.set_defaults(func=build_inventory)

    parse = subparsers.add_parser("parse", help="parse assigned map_record files")
    parse.add_argument("--manifest", required=True)
    parse.add_argument("--output-dir", required=True)
    parse.add_argument("--task-id", type=int)
    parse.add_argument("--session-filter")
    parse.add_argument("--record-filter")
    parse.add_argument("--limit", type=int)
    parse.add_argument("--resume", action="store_true")
    parse.set_defaults(func=parse_files)

    reduce = subparsers.add_parser("reduce", help="combine parsed charge products")
    reduce.add_argument("--manifest", required=True)
    reduce.add_argument("--parsed-dir", required=True)
    reduce.add_argument("--output-dir", required=True)
    reduce.add_argument("--allow-partial", action="store_true")
    reduce.set_defaults(func=reduce_results)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = make_parser()
    args = parser.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
