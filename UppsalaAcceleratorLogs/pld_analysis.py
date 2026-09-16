#!/usr/bin/env python3
"""Parse an IBA Scanalgo PLD plan and reconstruct planned isocentre fluence."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from UppsalaAcceleratorLogs.uppsala_analysis import (
    RADEX_INSTRUMENT_WIDTH_X_MM,
    RADEX_INSTRUMENT_WIDTH_Y_MM,
    facility_ppmu,
    finalize_heatmap,
    json_dump,
    parse_conversion_csv,
    save_heatmap_png,
)


SPOT_SIGMA_MM = {
    60.0: (7.0, 7.0),
    65.0: (6.6, 6.6),
    70.0: (6.3, 6.3),
    75.0: (6.0, 6.0),
    80.0: (5.7, 5.7),
    85.0: (5.4, 5.4),
}


@dataclass(frozen=True)
class PldSpot:
    layer_index: int
    spot_id: str
    energy_mev: float
    paintings: int
    x_iec_mm: float
    y_iec_mm: float
    meterset_weight: float
    meterset_rate: float

    @property
    def machine_x_mm(self) -> float:
        # Scanalgo writes IEC transverse coordinates.  The campaign recorder
        # heatmaps use the IBA machine axes, whose transverse X/Y axes are
        # interchanged for this beam orientation.
        return self.y_iec_mm

    @property
    def machine_y_mm(self) -> float:
        return self.x_iec_mm


@dataclass(frozen=True)
class PldPlan:
    source: str
    total_mu: float
    declared_cumulative_weight: float
    declared_layer_count: int
    spots: tuple[PldSpot, ...]
    layers: tuple[dict[str, float | int | str], ...]

    @property
    def parsed_weight_sum(self) -> float:
        return math.fsum(spot.meterset_weight for spot in self.spots)


def _finite_float(value: str, field: str, line_number: int) -> float:
    try:
        number = float(value)
    except ValueError as exc:
        raise ValueError(f"Line {line_number}: invalid {field}: {value!r}") from exc
    if not math.isfinite(number):
        raise ValueError(f"Line {line_number}: non-finite {field}: {value!r}")
    return number


def parse_pld(path: Path) -> PldPlan:
    """Parse and validate the Beam/Layer/Element PLD format.

    Each planned spot is represented by two consecutive Element records: a
    zero-weight positioning element followed by the delivered-weight element.
    The returned plan contains one collapsed :class:`PldSpot` per pair.
    """
    path = Path(path)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = [row for row in csv.reader(handle) if row]
    if not rows or rows[0][0].strip() != "Beam" or len(rows[0]) < 10:
        raise ValueError(f"{path} does not start with a valid PLD Beam record")

    header = rows[0]
    total_mu = _finite_float(header[-3], "total MU", 1)
    declared_cumulative = _finite_float(header[-2], "cumulative meterset weight", 1)
    try:
        declared_layers = int(header[-1])
    except ValueError as exc:
        raise ValueError(f"Line 1: invalid layer count: {header[-1]!r}") from exc
    if total_mu <= 0 or declared_cumulative <= 0 or declared_layers <= 0:
        raise ValueError("PLD Beam totals and layer count must be positive")

    layers: list[dict[str, float | int | str]] = []
    spots: list[PldSpot] = []
    row_index = 1
    while row_index < len(rows):
        line_number = row_index + 1
        layer_row = rows[row_index]
        if layer_row[0].strip() != "Layer" or len(layer_row) < 6:
            raise ValueError(f"Line {line_number}: expected a six-field Layer record")
        spot_id = layer_row[1].strip()
        energy = _finite_float(layer_row[2], "energy", line_number)
        cumulative = _finite_float(layer_row[3], "layer cumulative weight", line_number)
        try:
            element_count = int(layer_row[4])
            paintings = int(layer_row[5])
        except ValueError as exc:
            raise ValueError(f"Line {line_number}: invalid element count or paintings") from exc
        if element_count <= 0 or element_count % 2 or paintings <= 0:
            raise ValueError(
                f"Line {line_number}: element count must be positive and even, "
                "and paintings must be positive"
            )
        layer_index = len(layers) + 1
        element_rows = rows[row_index + 1 : row_index + 1 + element_count]
        if len(element_rows) != element_count:
            raise ValueError(f"Layer {layer_index}: truncated Element record list")
        layer_weight = 0.0
        for pair_offset in range(0, element_count, 2):
            pair = element_rows[pair_offset : pair_offset + 2]
            first_line = row_index + pair_offset + 2
            parsed: list[tuple[float, float, float, float]] = []
            for element_offset, element in enumerate(pair):
                element_line = first_line + element_offset
                if element[0].strip() != "Element" or len(element) < 5:
                    raise ValueError(f"Line {element_line}: invalid Element record")
                parsed.append(
                    (
                        _finite_float(element[1], "IEC X", element_line),
                        _finite_float(element[2], "IEC Y", element_line),
                        _finite_float(element[3], "meterset weight", element_line),
                        _finite_float(element[4], "meterset rate", element_line),
                    )
                )
            move, delivery = parsed
            if move[:2] != delivery[:2]:
                raise ValueError(f"Lines {first_line}-{first_line + 1}: spot coordinates differ")
            if not math.isclose(move[2], 0.0, abs_tol=1e-12):
                raise ValueError(f"Line {first_line}: positioning element has nonzero weight")
            if delivery[2] < 0:
                raise ValueError(f"Line {first_line + 1}: delivered weight is negative")
            layer_weight += delivery[2]
            spots.append(
                PldSpot(
                    layer_index=layer_index,
                    spot_id=spot_id,
                    energy_mev=energy,
                    paintings=paintings,
                    x_iec_mm=delivery[0],
                    y_iec_mm=delivery[1],
                    meterset_weight=delivery[2],
                    meterset_rate=delivery[3],
                )
            )
        layers.append(
            {
                "layer_index": layer_index,
                "spot_id": spot_id,
                "energy_mev": energy,
                "declared_cumulative_weight": cumulative,
                "parsed_weight_sum": layer_weight,
                "element_count": element_count,
                "spot_count": element_count // 2,
                "paintings": paintings,
            }
        )
        row_index += 1 + element_count

    if len(layers) != declared_layers:
        raise ValueError(
            f"Declared {declared_layers} layer(s), parsed {len(layers)}"
        )
    if not spots or sum(int(layer["element_count"]) for layer in layers) != 2 * len(spots):
        raise ValueError("PLD contains no complete spot pairs")
    return PldPlan(
        source=str(path.resolve()),
        total_mu=total_mu,
        declared_cumulative_weight=declared_cumulative,
        declared_layer_count=declared_layers,
        spots=tuple(spots),
        layers=tuple(layers),
    )


def spot_sigma(energy_mev: float) -> tuple[float, float]:
    for table_energy, sigma in SPOT_SIGMA_MM.items():
        if math.isclose(energy_mev, table_energy, abs_tol=1e-9):
            return sigma
    raise ValueError(
        f"No supplied isocentre spot sigma for {energy_mev:g} MeV; "
        f"available energies are {sorted(SPOT_SIGMA_MM)}"
    )


def reconstruct_plan(
    plan: PldPlan,
    conversion: list[dict[str, float]],
    grid: dict[str, float | str],
) -> tuple[dict[str, object], dict[str, np.ndarray], list[dict[str, object]]]:
    """Convert parsed spots into an absolute planned-fluence heatmap."""
    grid_min = float(grid["minimum_mm"])
    grid_max = float(grid["maximum_mm"])
    pixel_mm = float(grid["pixel_mm"])
    n_pixels = int(round((grid_max - grid_min) / pixel_mm))
    if n_pixels <= 0 or not math.isclose(grid_min + n_pixels * pixel_mm, grid_max):
        raise ValueError("Grid extent must contain a positive whole number of pixels")
    if plan.parsed_weight_sum <= 0:
        raise ValueError("PLD has no positive delivered meterset weight")

    width_histograms: dict[tuple[float, float], np.ndarray] = {}
    spot_rows: list[dict[str, object]] = []
    total_protons = 0.0
    for spot in plan.spots:
        sigma_x_iec, sigma_y_iec = spot_sigma(spot.energy_mev)
        # The axes swap along with the spot centre coordinates.
        sigma_x_machine, sigma_y_machine = sigma_y_iec, sigma_x_iec
        histogram = width_histograms.setdefault(
            (sigma_x_machine, sigma_y_machine),
            np.zeros((n_pixels, n_pixels), dtype=np.float64),
        )
        protons_per_mu, calibration_energy = facility_ppmu(spot.energy_mev, conversion)
        # Normalize by the parsed weights so the heatmap integrates to the Beam
        # record's prescribed total MU.  This also handles rounded PLD weights;
        # the declared-vs-parsed difference remains visible in summary.json.
        planned_mu = plan.total_mu * spot.meterset_weight / plan.parsed_weight_sum
        planned_protons = planned_mu * protons_per_mu
        column = int(math.floor((spot.machine_x_mm - grid_min) / pixel_mm))
        row = int(math.floor((spot.machine_y_mm - grid_min) / pixel_mm))
        if not (0 <= row < n_pixels and 0 <= column < n_pixels):
            raise ValueError(
                f"Spot at machine ({spot.machine_x_mm}, {spot.machine_y_mm}) mm "
                "falls outside the reconstruction grid"
            )
        histogram[row, column] += planned_protons
        total_protons += planned_protons
        spot_rows.append(
            {
                "layer_index": spot.layer_index,
                "spot_id": spot.spot_id,
                "energy_mev": spot.energy_mev,
                "paintings": spot.paintings,
                "x_iec_mm": spot.x_iec_mm,
                "y_iec_mm": spot.y_iec_mm,
                "machine_x_mm": spot.machine_x_mm,
                "machine_y_mm": spot.machine_y_mm,
                "meterset_weight": spot.meterset_weight,
                "planned_mu": planned_mu,
                "planned_protons": planned_protons,
                "calibration_energy_mev": calibration_energy,
                "protons_per_mu": protons_per_mu,
                "sigma_x_machine_mm": sigma_x_machine,
                "sigma_y_machine_mm": sigma_y_machine,
            }
        )

    protons_per_pixel = finalize_heatmap(width_histograms, pixel_mm, "sigma")
    pixel_area_cm2 = pixel_mm * pixel_mm / 100.0
    fluence = protons_per_pixel / pixel_area_cm2
    energies = sorted({spot.energy_mev for spot in plan.spots})
    summary: dict[str, object] = {
        "source_pld": plan.source,
        "coordinate_transform": "machine_x = pld_iec_y; machine_y = pld_iec_x",
        "device_plane": "isocentre",
        "geometric_magnification": 1.0,
        "declared_layer_count": plan.declared_layer_count,
        "parsed_layer_count": len(plan.layers),
        "spot_count": len(plan.spots),
        "energies_mev": energies,
        "total_mu": plan.total_mu,
        "declared_cumulative_meterset_weight": plan.declared_cumulative_weight,
        "parsed_spot_weight_sum": plan.parsed_weight_sum,
        "weight_sum_minus_declared": plan.parsed_weight_sum - plan.declared_cumulative_weight,
        "weight_normalization": "parsed spot-weight sum, scaled to Beam total MU",
        "planned_protons": total_protons,
        "heatmap_integral_protons": float(protons_per_pixel.sum()),
        "heatmap_conservation_relative_error": (
            float(protons_per_pixel.sum()) - total_protons
        ) / total_protons,
        "pixel_area_cm2": pixel_area_cm2,
        "grid": grid,
        "layers": list(plan.layers),
        "spot_sigma_table_mm": {
            f"{energy:g}": {"sigma_x": sigma[0], "sigma_y": sigma[1]}
            for energy, sigma in SPOT_SIGMA_MM.items()
        },
        "radex_outline_mm": {
            "machine_x": RADEX_INSTRUMENT_WIDTH_X_MM,
            "machine_y": RADEX_INSTRUMENT_WIDTH_Y_MM,
        },
    }
    x_edges = np.linspace(grid_min, grid_max, n_pixels + 1, dtype=np.float64)
    arrays = {
        "planned_protons_per_pixel": protons_per_pixel.astype(np.float32),
        "planned_fluence_protons_cm2": fluence.astype(np.float32),
        "x_edges_mm": x_edges,
        "y_edges_mm": x_edges.copy(),
    }
    return summary, arrays, spot_rows


def write_spots_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("No parsed spots to write")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> None:
    plan = parse_pld(Path(args.pld))
    conversion = parse_conversion_csv(Path(args.conversion_csv))
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    grid: dict[str, float | str] = {
        "minimum_mm": args.grid_min_mm,
        "maximum_mm": args.grid_max_mm,
        "pixel_mm": args.pixel_mm,
        "width_interpretation": "sigma",
    }
    summary, arrays, spots = reconstruct_plan(plan, conversion, grid)
    json_dump(output_dir / "summary.json", summary)
    write_spots_csv(output_dir / "planned_spots.csv", spots)
    np.savez_compressed(output_dir / "planned_fluence.npz", **arrays)
    title = (
        f"{Path(args.pld).stem} — planned {summary['energies_mev'][0]:g} MeV "
        "fluence at isocentre"
    )
    save_heatmap_png(
        output_dir / "planned_fluence_with_instrument_outline.png",
        arrays["planned_fluence_protons_cm2"],
        grid,
        title,
        "Planned protons cm$^{-2}$",
        centered_rectangle_mm=(
            RADEX_INSTRUMENT_WIDTH_X_MM,
            RADEX_INSTRUMENT_WIDTH_Y_MM,
        ),
        rectangle_label="RadEx instrument (154 x 66.8 mm)",
    )
    print(
        f"Parsed {len(plan.spots)} spots in {len(plan.layers)} layer(s); "
        f"wrote planned isocentre products to {output_dir}"
    )


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pld", required=True)
    parser.add_argument("--conversion-csv", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--grid-min-mm", type=float, default=-250.0)
    parser.add_argument("--grid-max-mm", type=float, default=250.0)
    parser.add_argument("--pixel-mm", type=float, default=1.0)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = make_parser().parse_args(argv)
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
