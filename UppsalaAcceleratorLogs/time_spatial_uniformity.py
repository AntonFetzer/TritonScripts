#!/usr/bin/env python3
"""Evaluate short-window spatial uniformity directly from PBS recorder rows."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
from scipy.ndimage import gaussian_filter

from UppsalaAcceleratorLogs import uppsala_analysis as ua


def epoch_second(timestamp: str, cache: dict[str, int]) -> int:
    # Milliseconds vary on every row, but this key changes only once per second.
    key = timestamp[:19] + timestamp[-6:]
    value = cache.get(key)
    if value is None:
        value = int(datetime.fromisoformat(key).timestamp())
        cache[key] = value
    return value


def field_geometry(final_map: np.ndarray, grid: dict[str, Any]) -> dict[str, float]:
    values = np.where(np.isfinite(final_map), final_map, 0.0).astype(np.float64)
    total = float(values.sum())
    rows, columns = np.indices(values.shape)
    center_column = float((columns * values).sum() / total)
    center_row = float((rows * values).sum() / total)
    width_x, width_y = ua.central_field_widths(values, grid)
    return {
        "center_x_mm": float(grid["minimum_mm"])
        + (center_column + 0.5) * float(grid["pixel_mm"]),
        "center_y_mm": float(grid["minimum_mm"])
        + (center_row + 0.5) * float(grid["pixel_mm"]),
        "width_x_50_mm": width_x,
        "width_y_50_mm": width_y,
    }


def parse_experiment_histogram(
    entries: list[dict[str, Any]],
    start_second: int,
    stop_second: int,
    minimum_mm: float,
    maximum_mm: float,
    pixel_mm: float,
) -> tuple[np.ndarray, dict[str, float]]:
    n_seconds = stop_second - start_second + 1
    n_pixels = int(round((maximum_mm - minimum_mm) / pixel_mm))
    histogram = np.zeros((n_seconds, n_pixels, n_pixels), dtype=np.float32)
    time_cache: dict[str, int] = {}
    weighted_x_width = 0.0
    weighted_y_width = 0.0
    weighted_x_width_squared = 0.0
    weighted_y_width_squared = 0.0
    mapped_charge = 0.0
    rows = 0

    for entry in entries:
        with Path(entry["record_path"]).open(
            "r", encoding="utf-8", errors="replace", newline=""
        ) as handle:
            reader = csv.reader(handle)
            in_stream = False
            for row in reader:
                if not row:
                    continue
                if row[0].strip() == "SUBMAP_NUMBER":
                    in_stream = True
                    continue
                if not in_stream or len(row) != 41:
                    continue
                try:
                    second = epoch_second(row[1].strip(), time_cache)
                    primary = float(row[26])
                    x_pos = float(row[4])
                    y_pos = float(row[5])
                    x_width = float(row[2])
                    y_width = float(row[3])
                except (ValueError, IndexError):
                    continue
                rows += 1
                primary = max(primary, 0.0)
                if (
                    primary <= 0.0
                    or primary >= 9999.0
                    or not start_second <= second <= stop_second
                    or not minimum_mm <= x_pos < maximum_mm
                    or not minimum_mm <= y_pos < maximum_mm
                    or x_width <= 0.0
                    or y_width <= 0.0
                ):
                    continue
                ti = second - start_second
                xi = int((x_pos - minimum_mm) // pixel_mm)
                yi = int((y_pos - minimum_mm) // pixel_mm)
                histogram[ti, yi, xi] += primary
                mapped_charge += primary
                weighted_x_width += primary * x_width
                weighted_y_width += primary * y_width
                weighted_x_width_squared += primary * x_width * x_width
                weighted_y_width_squared += primary * y_width * y_width

    if mapped_charge <= 0:
        raise RuntimeError("No valid positive spatial charge was found")
    sigma_x_mm = weighted_x_width / mapped_charge
    sigma_y_mm = weighted_y_width / mapped_charge
    sigma_x_std_mm = math.sqrt(
        max(weighted_x_width_squared / mapped_charge - sigma_x_mm**2, 0.0)
    )
    sigma_y_std_mm = math.sqrt(
        max(weighted_y_width_squared / mapped_charge - sigma_y_mm**2, 0.0)
    )
    gaussian_filter(
        histogram,
        sigma=(0.0, sigma_y_mm / pixel_mm, sigma_x_mm / pixel_mm),
        mode="constant",
        truncate=5.0,
        output=histogram,
    )
    return histogram, {
        "parsed_rows": rows,
        "mapped_positive_charge_c": mapped_charge,
        "charge_weighted_sigma_x_mm": sigma_x_mm,
        "charge_weighted_sigma_y_mm": sigma_y_mm,
        "charge_weighted_sigma_x_std_mm": sigma_x_std_mm,
        "charge_weighted_sigma_y_std_mm": sigma_y_std_mm,
    }


def analyze_windows(
    maps: np.ndarray,
    start_second: int,
    geometry: dict[str, float],
    minimum_mm: float,
    pixel_mm: float,
    windows: tuple[int, ...],
    roi_fractions: tuple[float, ...],
) -> list[dict[str, float]]:
    n_pixels = maps.shape[1]
    centers = minimum_mm + (np.arange(n_pixels) + 0.5) * pixel_mm
    xx, yy = np.meshgrid(centers, centers)
    masks = {}
    for fraction in roi_fractions:
        masks[fraction] = (
            np.abs(xx - geometry["center_x_mm"])
            <= geometry["width_x_50_mm"] * fraction / 2.0
        ) & (
            np.abs(yy - geometry["center_y_mm"])
            <= geometry["width_y_50_mm"] * fraction / 2.0
        )

    results: list[dict[str, float]] = []
    second_totals = maps.sum(axis=(1, 2), dtype=np.float64)
    total_cumulative = np.concatenate(([0.0], np.cumsum(second_totals)))
    for fraction, mask in masks.items():
        roi_seconds = maps[:, mask].astype(np.float64)
        roi_cumulative = np.vstack(
            (
                np.zeros((1, roi_seconds.shape[1]), dtype=np.float64),
                np.cumsum(roi_seconds, axis=0),
            )
        )
        for window in windows:
            if window > maps.shape[0]:
                continue
            rolling_values = roi_cumulative[window:] - roi_cumulative[:-window]
            rolling_totals = (
                total_cumulative[window:] - total_cumulative[:-window]
            )
            means = rolling_values.mean(axis=1)
            minimums = rolling_values.min(axis=1)
            maximums = rolling_values.max(axis=1)
            standard_deviations = rolling_values.std(axis=1)
            valid = (rolling_totals > 0.0) & (means > 0.0)
            for offset in np.flatnonzero(valid):
                mean = float(means[offset])
                minimum = float(minimums[offset])
                maximum = float(maximums[offset])
                max_deviation = max(
                    maximum / mean - 1.0, 1.0 - minimum / mean
                )
                results.append(
                    {
                        "window_seconds": float(window),
                        "window_end_epoch_second": float(
                            start_second + offset + window - 1
                        ),
                        "roi_fraction_of_halfmax_width": fraction,
                        "total_spatial_charge_c": float(rolling_totals[offset]),
                        "roi_mean_charge_density": mean,
                        "roi_min_over_mean": minimum / mean,
                        "roi_max_over_mean": maximum / mean,
                        "roi_max_abs_deviation_fraction": max_deviation,
                        "roi_coefficient_of_variation": float(
                            standard_deviations[offset] / mean
                        ),
                    }
                )
    return results


def summarize(rows: list[dict[str, float]]) -> list[dict[str, float]]:
    grouped: dict[tuple[int, float], list[dict[str, float]]] = defaultdict(list)
    for row in rows:
        grouped[
            (
                int(row["window_seconds"]),
                float(row["roi_fraction_of_halfmax_width"]),
            )
        ].append(row)
    result = []
    for (window, fraction), values in sorted(grouped.items()):
        charges = np.asarray([row["total_spatial_charge_c"] for row in values])
        deviations = np.asarray(
            [row["roi_max_abs_deviation_fraction"] for row in values]
        )
        # "Substantial" excludes only very low-charge leading/trailing windows:
        # at least 10% of the median among the upper half of window charges.
        upper = charges[charges >= np.median(charges)]
        threshold = 0.1 * float(np.median(upper))
        selected = deviations[charges >= threshold]
        result.append(
            {
                "window_seconds": window,
                "roi_fraction_of_halfmax_width": fraction,
                "all_active_window_count": int(len(deviations)),
                "substantial_window_charge_threshold_c": threshold,
                "substantial_window_count": int(len(selected)),
                "substantial_fraction_within_1pct": float(np.mean(selected <= 0.01)),
                "substantial_median_max_deviation_fraction": float(
                    np.median(selected)
                ),
                "substantial_p95_max_deviation_fraction": float(
                    np.percentile(selected, 95)
                ),
                "substantial_worst_max_deviation_fraction": float(selected.max()),
            }
        )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--experiment-index", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--domain-min-mm", type=float, default=-160.0)
    parser.add_argument("--domain-max-mm", type=float, default=160.0)
    parser.add_argument("--pixel-mm", type=float, default=1.0)
    args = parser.parse_args()

    manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    experiments = manifest["experiments"]
    experiment = experiments[args.experiment_index]
    experiment_id = experiment["id"]
    entries = [
        entry
        for entry in manifest["entries"]
        if entry.get("experiment_id") == experiment_id
        and entry.get("kind") == "part"
    ]
    if not entries:
        raise RuntimeError(f"No part records found for {experiment_id}")

    start_second = int(datetime.fromisoformat(experiment["start"]).timestamp())
    stop_second = int(datetime.fromisoformat(experiment["stop"]).timestamp())
    final_path = Path(args.results_dir) / "experiments" / f"{experiment_id}.npz"
    with np.load(final_path, allow_pickle=False) as data:
        final_map = data["fluence_protons_cm2"]
    geometry = field_geometry(final_map, manifest["grid"])

    maps, parse_summary = parse_experiment_histogram(
        entries,
        start_second,
        stop_second,
        args.domain_min_mm,
        args.domain_max_mm,
        args.pixel_mm,
    )
    rows = analyze_windows(
        maps,
        start_second,
        geometry,
        args.domain_min_mm,
        args.pixel_mm,
        windows=(10, 30, 60, 90, 120, 180, 240, 300, 420, 540),
        roi_fractions=(0.25, 0.50, 0.75, 1.00),
    )
    summaries = summarize(rows)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f"{experiment_id}_window_uniformity.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "experiment": experiment,
        "entries": len(entries),
        "geometry": geometry,
        "parse": parse_summary,
        "window_definition": "Trailing window evaluated at every integer second",
        "uniformity_metric": (
            "Maximum absolute pixel deviation from the spatial mean within "
            "the centered ROI rectangle"
        ),
        "roi_definition": "Centered fractions of measured half-maximum X/Y widths",
        "summaries": summaries,
    }
    (output_dir / f"{experiment_id}_summary.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
