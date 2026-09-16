#!/usr/bin/env python3
"""Measure complete-field repaint cadence from parsed accelerator log parts."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from datetime import datetime
from pathlib import Path

import numpy as np

from UppsalaAcceleratorLogs import uppsala_analysis as ua


def correlation(left: np.ndarray, right: np.ndarray) -> float:
    a = left.astype(np.float64, copy=False).ravel()
    b = right.astype(np.float64, copy=False).ravel()
    a -= a.mean()
    b -= b.mean()
    denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / denominator) if denominator else float("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--parsed-dir", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--complete-duration-seconds", type=float, default=30.0)
    args = parser.parse_args()

    manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    parsed_dir = Path(args.parsed_dir)
    results_dir = Path(args.results_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    summaries = []
    for path in parsed_dir.glob("*.json"):
        row = json.loads(path.read_text(encoding="utf-8"))
        if (
            row.get("experiment_id")
            and row.get("kind") == "part"
            and row.get("first_time")
            and row.get("last_time")
        ):
            row["_start"] = datetime.fromisoformat(row["first_time"])
            row["_stop"] = datetime.fromisoformat(row["last_time"])
            row["_duration"] = (row["_stop"] - row["_start"]).total_seconds()
            summaries.append(row)

    detail_rows = []
    experiment_rows = []
    for experiment in manifest["experiments"]:
        experiment_id = experiment["id"]
        rows = sorted(
            [row for row in summaries if row["experiment_id"] == experiment_id],
            key=lambda row: row["_start"],
        )
        complete = [
            row
            for row in rows
            if row["_duration"] >= args.complete_duration_seconds
        ]
        with np.load(
            results_dir / "experiments" / f"{experiment_id}.npz",
            allow_pickle=False,
        ) as data:
            reference = data["fluence_protons_cm2"].astype(np.float64)

        correlations = []
        normalized_maps: dict[int, np.ndarray] = {}
        grid = manifest["grid"]
        minimum_mm = float(grid["minimum_mm"])
        pixel_mm = float(grid["pixel_mm"])
        for row in complete:
            with np.load(
                parsed_dir / f"{row['id']}.npz", allow_pickle=False
            ) as data:
                heatmap = data["heatmap_charge_c"]
            value = correlation(heatmap, reference)
            correlations.append(value)
            heatmap64 = heatmap.astype(np.float64)
            map_sum = float(heatmap64.sum())
            yy, xx = np.indices(heatmap64.shape)
            center_x_mm = (
                minimum_mm
                + (float((xx * heatmap64).sum() / map_sum) + 0.5) * pixel_mm
            )
            center_y_mm = (
                minimum_mm
                + (float((yy * heatmap64).sum() / map_sum) + 0.5) * pixel_mm
            )
            width_x, width_y = ua.central_field_widths(heatmap64, grid)
            normalized_maps[int(row["sequence"])] = heatmap64 / map_sum
            detail_rows.append(
                {
                    "experiment_id": experiment_id,
                    "sequence": row["sequence"],
                    "first_time": row["first_time"],
                    "last_time": row["last_time"],
                    "duration_seconds": row["_duration"],
                    "spatial_map_correlation_to_experiment": value,
                    "map_center_x_mm": center_x_mm,
                    "map_center_y_mm": center_y_mm,
                    "map_width_x_50_mm": width_x,
                    "map_width_y_50_mm": width_y,
                    "element_count": row.get("element_count"),
                    "acquisition_rows": row.get("acquisition_rows"),
                    "spatialized_positive_charge_c": row.get(
                        "charge_spatialized_positive_c"
                    ),
                }
            )

        lag_summaries = []
        for lag in range(1, min(31, len(rows))):
            lag_correlations = []
            lag_periods = []
            by_sequence = {int(row["sequence"]): row for row in complete}
            for sequence, left_map in normalized_maps.items():
                right_map = normalized_maps.get(sequence + lag)
                if right_map is None:
                    continue
                lag_correlations.append(correlation(left_map, right_map))
                lag_periods.append(
                    (
                        by_sequence[sequence + lag]["_start"]
                        - by_sequence[sequence]["_start"]
                    ).total_seconds()
                )
            if lag_correlations:
                lag_summaries.append(
                    {
                        "lag_parts": lag,
                        "pair_count": len(lag_correlations),
                        "median_map_correlation": statistics.median(
                            lag_correlations
                        ),
                        "minimum_map_correlation": min(lag_correlations),
                        "median_elapsed_seconds": statistics.median(lag_periods),
                    }
                )
        repeating_candidates = [
            item for item in lag_summaries if item["pair_count"] >= 2
        ]
        high_correlation_candidates = [
            item
            for item in repeating_candidates
            if item["median_map_correlation"] >= 0.98
        ]
        fundamental_repeating = (
            min(high_correlation_candidates, key=lambda item: item["lag_parts"])
            if high_correlation_candidates
            else None
        )
        best_repeating = (
            max(
                repeating_candidates,
                key=lambda item: item["median_map_correlation"],
            )
            if repeating_candidates
            else None
        )

        consecutive_complete_pairs = [
            (left, right)
            for left, right in zip(rows, rows[1:])
            if left["_duration"] >= args.complete_duration_seconds
            and right["_duration"] >= args.complete_duration_seconds
        ]
        periods = [
            (right["_start"] - left["_start"]).total_seconds()
            for left, right in consecutive_complete_pairs
        ]
        gaps = [
            (right["_start"] - left["_stop"]).total_seconds()
            for left, right in consecutive_complete_pairs
        ]
        durations = [row["_duration"] for row in complete]
        period = statistics.median(periods) if periods else float("nan")
        experiment_rows.append(
            {
                "experiment_id": experiment_id,
                "energy_mev": experiment["energy_mev"],
                "part_count": len(rows),
                "complete_repaint_count": len(complete),
                "partial_or_aborted_part_count": len(rows) - len(complete),
                "median_complete_repaint_duration_seconds": statistics.median(
                    durations
                ),
                "median_complete_repaint_period_seconds": period,
                "repaint_frequency_hz": 1.0 / period if periods else float("nan"),
                "repaints_per_minute": 60.0 / period if periods else float("nan"),
                "median_inter_repaint_gap_seconds": (
                    statistics.median(gaps) if gaps else float("nan")
                ),
                "median_complete_map_correlation": statistics.median(
                    correlations
                ),
                "minimum_complete_map_correlation": min(correlations),
                "recorder_sample_rate_hz": statistics.median(
                    row["acquisition_rows"] / row["_duration"] for row in complete
                ),
                "median_map_elements": statistics.median(
                    row["element_count"] for row in complete
                ),
                "best_repeating_map_lag_parts": (
                    best_repeating["lag_parts"] if best_repeating else None
                ),
                "best_repeating_map_period_seconds": (
                    best_repeating["median_elapsed_seconds"]
                    if best_repeating
                    else None
                ),
                "best_repeating_map_correlation": (
                    best_repeating["median_map_correlation"]
                    if best_repeating
                    else None
                ),
                "fundamental_high_correlation_lag_parts": (
                    fundamental_repeating["lag_parts"]
                    if fundamental_repeating
                    else None
                ),
                "fundamental_high_correlation_period_seconds": (
                    fundamental_repeating["median_elapsed_seconds"]
                    if fundamental_repeating
                    else None
                ),
                "fundamental_high_correlation_value": (
                    fundamental_repeating["median_map_correlation"]
                    if fundamental_repeating
                    else None
                ),
                "lag_summaries_json": json.dumps(lag_summaries),
            }
        )

    for filename, rows in (
        ("scan_cadence_summary.csv", experiment_rows),
        ("scan_cadence_parts.csv", detail_rows),
    ):
        with (output_dir / filename).open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    (output_dir / "scan_cadence_report.json").write_text(
        json.dumps(
            {
                "complete_repaint_definition": (
                    f"part duration >= {args.complete_duration_seconds:g} s"
                ),
                "interpretation": (
                    "A complete part is treated as one full-field repaint when "
                    "its normalized spatial map correlates with the experiment's "
                    "integrated map."
                ),
                "experiments": experiment_rows,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps(experiment_rows, indent=2))


if __name__ == "__main__":
    main()
