#!/usr/bin/env python3
"""Create final RadEx isocentre fluence maps and location tables."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from UppsalaAcceleratorLogs.register_fluence_maps import (
    load_heatmap,
    registered_full_grid,
)
from UppsalaAcceleratorLogs.uppsala_analysis import (
    RADEX_INSTRUMENT_WIDTH_X_MM,
    RADEX_INSTRUMENT_WIDTH_Y_MM,
)


# Best PLD-to-aggregate registration.  Definition:
# x_iso = scale_x * x_chamber + offset_x, likewise for y.
RADEX_CHAMBER_TO_ISOCENTRE = {
    "scale_x_chamber_to_isocentre": 1.480648361103907,
    "scale_y_chamber_to_isocentre": 1.3288011066256862,
    "offset_x_isocentre_mm": 0.6111479120787414,
    "offset_y_isocentre_mm": 2.3558219638048454,
    "area_scale_chamber_to_isocentre": 1.9674871807583802,
    "fluence_density_jacobian_chamber_to_isocentre": 0.5082625237814986,
    "pearson_correlation": 0.99234553966844,
}


RADEX_RADFETS = (
    (0, "VT01, 100% Pb", 60.0, 16.7),
    (1, "VT01, 75% Pb / 25% PE", 30.0, 16.7),
    (2, "VT01, 50% Pb / 50% PE", 0.0, 16.7),
    (3, "VT01, 25% Pb / 75% PE", -30.0, 16.7),
    (4, "VT01, 100% PE", -60.0, 16.7),
    (5, "VT01, exposed", 60.0, -16.7),
    (6, "VT01, 1 mm Al", 30.0, -16.7),
    (7, "VT05, 2 mm Al", 0.0, -16.7),
    (8, "VT05, 4 mm Al", -30.0, -16.7),
    (9, "VT05, 6 mm Al", -60.0, -16.7),
    (10, "VT01, readout PCB back", 28.0, -20.0),
    (11, "VT05, readout PCB back", 0.0, 0.0),
)


def bilinear_sample(
    values: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    x_mm: float,
    y_mm: float,
) -> float:
    """Sample cell-centred values with bilinear interpolation."""
    pixel_x = float(np.diff(x_edges).mean())
    pixel_y = float(np.diff(y_edges).mean())
    column = (x_mm - x_edges[0]) / pixel_x - 0.5
    row = (y_mm - y_edges[0]) / pixel_y - 0.5
    column0 = int(math.floor(column))
    row0 = int(math.floor(row))
    if not (0 <= column0 < values.shape[1] - 1 and 0 <= row0 < values.shape[0] - 1):
        raise ValueError(f"Sample point ({x_mm}, {y_mm}) mm is outside the heatmap")
    dx = column - column0
    dy = row - row0
    return float(
        (1.0 - dx) * (1.0 - dy) * values[row0, column0]
        + dx * (1.0 - dy) * values[row0, column0 + 1]
        + (1.0 - dx) * dy * values[row0 + 1, column0]
        + dx * dy * values[row0 + 1, column0 + 1]
    )


def location_rows(
    fluence: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
) -> list[dict[str, Any]]:
    rows = []
    for channel, shielding, x_mm, y_mm in RADEX_RADFETS:
        rows.append(
            {
                "channel": channel,
                "shielding": shielding,
                "gdml_x_mm": x_mm,
                "gdml_y_mm": y_mm,
                "delivered_fluence_protons_cm2": bilinear_sample(
                    fluence, x_edges, y_edges, x_mm, y_mm
                ),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(
    path: Path,
    title: str,
    rows: list[dict[str, Any]],
) -> None:
    lines = [
        f"# {title}: delivered fluence at RadFET locations",
        "",
        (
            "Values are bilinearly sampled from the ionisation-chamber heatmap "
            "after transformation to the isocentre. GDML coordinates are assumed "
            "to be nominally centred and aligned with the machine axes; no "
            "hypothetical DUT translation or rotation is applied."
        ),
        "",
        "| Ch. | Shielding | GDML X [mm] | GDML Y [mm] | Delivered fluence [protons/cm2] |",
        "| ---: | --- | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row['channel']} | {row['shielding']} | "
            f"{row['gdml_x_mm']:+g} | {row['gdml_y_mm']:+g} | "
            f"{row['delivered_fluence_protons_cm2']:.6e} |"
        )
    lines.extend(
        [
            "",
            (
                "The chamber-to-isocentre coordinate scale is 1.480648 in X "
                "and 1.328801 in Y. Fluence density includes the inverse area "
                "Jacobian 0.508263, which preserves the integrated proton count."
            ),
            "",
            (
                "These are incident proton fluences, not absorbed-dose values. "
                "Absorbed dose requires the channel- and energy-dependent Geant4 "
                "dose coefficient."
            ),
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8", newline="\n")


def save_heatmap(
    path: Path,
    title: str,
    fluence: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    rows: list[dict[str, Any]],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    fig, ax = plt.subplots(figsize=(9.0, 7.2), constrained_layout=True)
    image = ax.imshow(
        fluence,
        origin="lower",
        extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]],
        cmap="magma",
    )
    ax.add_patch(
        Rectangle(
            (-RADEX_INSTRUMENT_WIDTH_X_MM / 2.0, -RADEX_INSTRUMENT_WIDTH_Y_MM / 2.0),
            RADEX_INSTRUMENT_WIDTH_X_MM,
            RADEX_INSTRUMENT_WIDTH_Y_MM,
            fill=False,
            edgecolor="limegreen",
            linewidth=1.4,
            label="RadEx outline (154 x 66.8 mm)",
        )
    )
    for row in rows:
        x_mm = float(row["gdml_x_mm"])
        y_mm = float(row["gdml_y_mm"])
        ax.plot(x_mm, y_mm, marker="o", markersize=4.5, color="cyan", markeredgecolor="black")
        ax.annotate(
            str(row["channel"]),
            (x_mm, y_mm),
            xytext=(4, 4),
            textcoords="offset points",
            color="cyan",
            fontsize=8,
            weight="bold",
        )
    ax.set_xlim(-100, 100)
    ax.set_ylim(-60, 60)
    ax.set_xlabel("Isocentre X [mm]")
    ax.set_ylabel("Isocentre Y [mm]")
    ax.set_title(title)
    ax.legend(loc="upper right")
    fig.colorbar(image, ax=ax, label="Delivered protons cm$^{-2}$")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def process_map(
    input_path: Path,
    output_dir: Path,
    title: str,
    grid_min_mm: float,
    grid_max_mm: float,
) -> dict[str, Any]:
    values, x_edges, y_edges = load_heatmap(
        input_path,
        ("fluence_protons_cm2",),
        grid_min_mm,
        grid_max_mm,
    )
    registered = registered_full_grid(
        values,
        x_edges,
        y_edges,
        x_edges,
        y_edges,
        RADEX_CHAMBER_TO_ISOCENTRE,
        normalize_input=False,
    ) * RADEX_CHAMBER_TO_ISOCENTRE[
        "fluence_density_jacobian_chamber_to_isocentre"
    ]
    pixel_area_cm2 = float(np.diff(x_edges).mean() * np.diff(y_edges).mean() / 100.0)
    input_integral = float(values.sum() * pixel_area_cm2)
    output_integral = float(registered.sum() * pixel_area_cm2)
    rows = location_rows(registered, x_edges, y_edges)
    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_dir / "heatmap_isocentre.npz",
        fluence_protons_cm2=registered.astype(np.float32),
        x_edges_mm=x_edges,
        y_edges_mm=y_edges,
        scale_x_chamber_to_isocentre=RADEX_CHAMBER_TO_ISOCENTRE[
            "scale_x_chamber_to_isocentre"
        ],
        scale_y_chamber_to_isocentre=RADEX_CHAMBER_TO_ISOCENTRE[
            "scale_y_chamber_to_isocentre"
        ],
        offset_x_isocentre_mm=RADEX_CHAMBER_TO_ISOCENTRE["offset_x_isocentre_mm"],
        offset_y_isocentre_mm=RADEX_CHAMBER_TO_ISOCENTRE["offset_y_isocentre_mm"],
        fluence_density_jacobian=RADEX_CHAMBER_TO_ISOCENTRE[
            "fluence_density_jacobian_chamber_to_isocentre"
        ],
    )
    write_csv(output_dir / "radfet_location_fluence.csv", rows)
    write_markdown(output_dir / "radfet_location_fluence.md", title, rows)
    save_heatmap(
        output_dir / "fluence_isocentre.png",
        f"{title} — delivered fluence at isocentre",
        registered,
        x_edges,
        y_edges,
        rows,
    )
    summary = {
        "title": title,
        "source_heatmap": str(input_path.resolve()),
        "transform": RADEX_CHAMBER_TO_ISOCENTRE,
        "input_integrated_protons": input_integral,
        "isocentre_integrated_protons": output_integral,
        "integral_ratio": output_integral / input_integral,
        "radfet_coordinate_assumption": (
            "Nominal GDML coordinates centred and aligned with machine X/Y; "
            "no additional DUT pose transform"
        ),
    }
    with (output_dir / "isocentre_summary.json").open(
        "w", encoding="utf-8", newline="\n"
    ) as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return summary


def parse_map(text: str) -> tuple[str, Path]:
    if "=" not in text:
        raise argparse.ArgumentTypeError("map must be TITLE=PATH")
    title, path = text.split("=", 1)
    return title, Path(path)


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--map", action="append", type=parse_map, required=True)
    parser.add_argument("--grid-min-mm", type=float, default=-250.0)
    parser.add_argument("--grid-max-mm", type=float, default=250.0)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = make_parser().parse_args(argv)
    for title, input_path in args.map:
        summary = process_map(
            input_path,
            input_path.parent,
            title,
            args.grid_min_mm,
            args.grid_max_mm,
        )
        print(
            f"{title}: wrote isocentre products; "
            f"integral ratio={summary['integral_ratio']:.8f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
