#!/usr/bin/env python3
"""Create final RadEx isocentre fluence maps and location tables.

Two spot models are available:

* ``facility`` (default since 2026-10-23): the unblurred, charge-weighted
  spot-position map of each aggregate (``heatmap_positions.npz``) is carried
  to the isocentre with the registered chamber-to-isocentre transform and only
  then blurred with the facility's isocentre spot sigma, interpolated linearly
  in energy from ``pld_analysis.SPOT_SIGMA_MM`` (Skandion email 2026-08-24).
  ``RadEx_Total`` is the sum of the 64 MeV and 85 MeV maps, each blurred with
  its own sigma.
* ``recorder`` (the method before 2026-10-23, kept for comparison): the
  recorder-plane ``heatmap.npz``, in which every sample is already blurred with
  the recorder X_WIDTH/Y_WIDTH as a chamber-plane sigma, is magnified as a
  whole.  This magnifies the blur by 1.48 / 1.33 (effective isocentre sigma
  11.7 x 9.0 mm at 64 MeV, 11.0 x 7.6 mm at 85 MeV), which has no physical
  basis: the point-source projection governs spot positions, not spot width.

Two figures that the channel-scatter audit of 2026-10-16 drew by hand are
produced here since 2026-10-26: ``row_profiles_spot_models.png`` in every
RadEx aggregate (the delivered field along the two RadFET rows under both
spot models, against the shared-area fluence) and
``spot_lattice_64_vs_85.png`` with ``spot_lattice_summary.json`` in
``RadEx_Total`` (the delivered spot lattice of the two energies in the
recorder frame, from the fine projections in ``heatmap_positions.npz``;
the equality of the two lattices is what justifies reusing the 85 MeV
registration at 64 MeV).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from UppsalaAcceleratorLogs.pld_analysis import SPOT_SIGMA_MM
from UppsalaAcceleratorLogs.register_fluence_maps import (
    load_heatmap,
    registered_full_grid,
)
from UppsalaAcceleratorLogs.uppsala_analysis import (
    RADEX_INSTRUMENT_AREA_CM2,
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
    (6, "VT01, 2 mm Al", 30.0, -16.7),
    (7, "VT05, 4 mm Al", 0.0, -16.7),
    (8, "VT05, 5 mm Al", -30.0, -16.7),
    (9, "VT05, 6 mm Al", -60.0, -16.7),
    (10, "VT01, readout PCB back", 28.0, -20.0),
    (11, "VT05, readout PCB back", 0.0, 0.0),
)


# Nominal energy of each RadEx aggregate's irradiations.  RadEx_Total is
# built from its two single-energy components so that each keeps its own sigma.
RADEX_AGGREGATE_COMPONENTS = {
    "RadEx_64MeV": (("RadEx_64MeV", 64.0),),
    "Radex_85MeV": (("Radex_85MeV", 85.0),),
    "RadEx_Total": (("RadEx_64MeV", 64.0), ("Radex_85MeV", 85.0)),
}

SPOT_SIGMA_SOURCE = (
    "Facility (Skandion) isocentre spot sigma (rms) per nominal energy, sent by "
    "email on 2026-08-24 (transcribed in the campaign's "
    "metadata/Skandion-email-2026-08-24-spot-sizes.md); table "
    "UppsalaAcceleratorLogs.pld_analysis.SPOT_SIGMA_MM, interpolated linearly "
    "in nominal energy. The same spot tune (3.0) was used at 64 and 85 MeV."
)

GAUSSIAN_TRUNCATE = 5.0

#: Isocentre y of the two RadFET rows (sensor board) and the half extent of
#: the planned spot-centre lattice (78 x 34 spots at 2 mm pitch, PLD).
RADFET_ROW_Y_MM = (16.7, -16.7)
PLD_SPOT_PITCH_MM = 2.0
PLANNED_HALF_EXTENT_MM = (77.0, 33.0)


def facility_spot_sigma(energy_mev: float) -> tuple[float, float]:
    """Facility isocentre sigma (x, y) in mm, linear in energy between table rows.

    No extrapolation: np.interp would silently clamp outside the table.
    """
    energies = np.asarray(sorted(SPOT_SIGMA_MM), dtype=float)
    if not (energies[0] - 1e-9 <= energy_mev <= energies[-1] + 1e-9):
        raise ValueError(
            f"{energy_mev:g} MeV is outside the facility spot-sigma table "
            f"({energies[0]:g}-{energies[-1]:g} MeV)"
        )
    sigma_x = np.asarray([SPOT_SIGMA_MM[e][0] for e in energies], dtype=float)
    sigma_y = np.asarray([SPOT_SIGMA_MM[e][1] for e in energies], dtype=float)
    return (
        float(np.interp(energy_mev, energies, sigma_x)),
        float(np.interp(energy_mev, energies, sigma_y)),
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


def sample_radfets(
    fluence: np.ndarray, x_edges: np.ndarray, y_edges: np.ndarray
) -> np.ndarray:
    return np.asarray(
        [
            bilinear_sample(fluence, x_edges, y_edges, x_mm, y_mm)
            for _channel, _label, x_mm, y_mm in RADEX_RADFETS
        ],
        dtype=float,
    )


def transform_positions_to_isocentre(
    values: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    out_x_edges: np.ndarray,
    out_y_edges: np.ndarray,
    transform: dict[str, float],
) -> tuple[np.ndarray, float]:
    """Carry an unblurred chamber-plane fluence map to the isocentre.

    Each chamber pixel's protons (fluence x chamber pixel area) are placed at
    the transformed pixel centre, x_iso = scale_x * x_chamber + offset_x, and
    shared over the four nearest isocentre pixels (cloud-in-cell); dividing by
    the isocentre pixel area gives fluence.  This is the same affine transform
    and inverse-Jacobian density scaling as ``registered_full_grid(...) *
    jacobian`` in ``process_map``, written as a forward deposit: resampling an
    unblurred point histogram by pull-back interpolation does not conserve
    charge (about -0.2 % globally and up to 1 % locally on a synthetic 2 mm
    spot lattice), whereas the forward deposit conserves it exactly.

    Returns the isocentre fluence map and the protons that fell outside the
    output grid.
    """
    pixel_in_x = float(np.diff(x_edges).mean())
    pixel_in_y = float(np.diff(y_edges).mean())
    pixel_out_x = float(np.diff(out_x_edges).mean())
    pixel_out_y = float(np.diff(out_y_edges).mean())
    protons = np.asarray(values, dtype=np.float64) * (pixel_in_x * pixel_in_y / 100.0)
    rows, columns = np.nonzero(protons)
    weights = protons[rows, columns]
    x_chamber = x_edges[0] + (columns + 0.5) * pixel_in_x
    y_chamber = y_edges[0] + (rows + 0.5) * pixel_in_y
    x_iso = (
        transform["scale_x_chamber_to_isocentre"] * x_chamber
        + transform["offset_x_isocentre_mm"]
    )
    y_iso = (
        transform["scale_y_chamber_to_isocentre"] * y_chamber
        + transform["offset_y_isocentre_mm"]
    )
    column_f = (x_iso - out_x_edges[0]) / pixel_out_x - 0.5
    row_f = (y_iso - out_y_edges[0]) / pixel_out_y - 0.5
    column0 = np.floor(column_f).astype(np.int64)
    row0 = np.floor(row_f).astype(np.int64)
    dx = column_f - column0
    dy = row_f - row0
    n_rows = len(out_y_edges) - 1
    n_columns = len(out_x_edges) - 1
    deposited = np.zeros((n_rows, n_columns), dtype=np.float64)
    outside = 0.0
    for d_column, d_row, fraction in (
        (0, 0, (1.0 - dx) * (1.0 - dy)),
        (1, 0, dx * (1.0 - dy)),
        (0, 1, (1.0 - dx) * dy),
        (1, 1, dx * dy),
    ):
        target_row = row0 + d_row
        target_column = column0 + d_column
        inside = (
            (target_row >= 0)
            & (target_row < n_rows)
            & (target_column >= 0)
            & (target_column < n_columns)
        )
        share = weights * fraction
        np.add.at(deposited, (target_row[inside], target_column[inside]), share[inside])
        outside += float(share[~inside].sum())
    return deposited / (pixel_out_x * pixel_out_y / 100.0), outside


def gaussian_blur_mm(
    values: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    sigma_x_mm: float,
    sigma_y_mm: float,
) -> np.ndarray:
    from scipy.ndimage import gaussian_filter

    pixel_x = float(np.diff(x_edges).mean())
    pixel_y = float(np.diff(y_edges).mean())
    return gaussian_filter(
        np.asarray(values, dtype=np.float64),
        sigma=(sigma_y_mm / pixel_y, sigma_x_mm / pixel_x),
        mode="constant",
        truncate=GAUSSIAN_TRUNCATE,
    )


def recorder_blur_isocentre_map(
    input_path: Path,
    grid_min_mm: float,
    grid_max_mm: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Superseded recorder-blur method: magnify the already-blurred chamber map."""
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
    return values, registered, x_edges, y_edges


def facility_sigma_isocentre_map(
    aggregates_root: Path,
    name: str,
    grid_min_mm: float,
    grid_max_mm: float,
) -> dict[str, Any]:
    """Facility-sigma isocentre map of a RadEx aggregate, with its checks."""
    if name not in RADEX_AGGREGATE_COMPONENTS:
        raise ValueError(
            f"No energy assignment for aggregate {name}; "
            f"known: {sorted(RADEX_AGGREGATE_COMPONENTS)}"
        )
    total = None
    cross_check = None
    components = []
    x_edges = y_edges = None
    for component, energy_mev in RADEX_AGGREGATE_COMPONENTS[name]:
        path = aggregates_root / component / "heatmap_positions.npz"
        if not path.exists():
            raise FileNotFoundError(
                f"{path} is missing: re-parse the RadEx part records and reduce "
                "(uppsala_analysis parse --experiments exp_01,...,exp_09 --kind part "
                "--no-resume; then reduce)"
            )
        with np.load(path) as data:
            if "energies_mev" in data.files:
                energies = [float(value) for value in data["energies_mev"]]
                if energies != [energy_mev]:
                    raise ValueError(
                        f"{path} holds energies {energies}, expected [{energy_mev}]"
                    )
        values, x_edges, y_edges = load_heatmap(
            path, ("fluence_protons_cm2",), grid_min_mm, grid_max_mm
        )
        pixel_area_cm2 = float(np.diff(x_edges).mean() * np.diff(y_edges).mean() / 100.0)
        transformed, outside = transform_positions_to_isocentre(
            values, x_edges, y_edges, x_edges, y_edges, RADEX_CHAMBER_TO_ISOCENTRE
        )
        sigma_x, sigma_y = facility_spot_sigma(energy_mev)
        blurred = gaussian_blur_mm(transformed, x_edges, y_edges, sigma_x, sigma_y)
        # Independent numerical check of the ordering: for an affine scaling,
        # blurring at the isocentre with sigma equals blurring in the chamber
        # plane with sigma / scale and then resampling as process_map does.
        chamber_blurred = gaussian_blur_mm(
            values,
            x_edges,
            y_edges,
            sigma_x / RADEX_CHAMBER_TO_ISOCENTRE["scale_x_chamber_to_isocentre"],
            sigma_y / RADEX_CHAMBER_TO_ISOCENTRE["scale_y_chamber_to_isocentre"],
        )
        pulled_back = registered_full_grid(
            chamber_blurred,
            x_edges,
            y_edges,
            x_edges,
            y_edges,
            RADEX_CHAMBER_TO_ISOCENTRE,
            normalize_input=False,
        ) * RADEX_CHAMBER_TO_ISOCENTRE["fluence_density_jacobian_chamber_to_isocentre"]
        components.append(
            {
                "aggregate": component,
                "source_heatmap": str(path.resolve()),
                "nominal_energy_mev": energy_mev,
                "spot_sigma_x_mm": sigma_x,
                "spot_sigma_y_mm": sigma_y,
                "chamber_integrated_protons": float(values.sum() * pixel_area_cm2),
                "transformed_integrated_protons": float(
                    transformed.sum() * pixel_area_cm2
                ),
                "protons_outside_grid_after_transform": outside,
                "blurred_integrated_protons": float(blurred.sum() * pixel_area_cm2),
            }
        )
        total = blurred if total is None else total + blurred
        cross_check = pulled_back if cross_check is None else cross_check + pulled_back
    assert total is not None and cross_check is not None
    assert x_edges is not None and y_edges is not None
    return {
        "fluence": total,
        "cross_check_fluence": cross_check,
        "x_edges": x_edges,
        "y_edges": y_edges,
        "components": components,
    }


def row_profile(
    fluence: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    y_mm: float,
    x_mm: np.ndarray,
) -> np.ndarray:
    """Bilinear samples of a map along one row of constant isocentre y."""
    return np.asarray(
        [bilinear_sample(fluence, x_edges, y_edges, float(x), y_mm) for x in x_mm],
        dtype=float,
    )


def save_row_profiles(
    path: Path,
    name: str,
    maps: list[tuple[str, np.ndarray, np.ndarray, np.ndarray]],
    shared_fluence: float,
) -> None:
    """The delivered field along the two RadFET rows under each spot model.

    One panel per RadFET row (isocentre y = +16.7 and -16.7 mm), the local
    fluence divided by the shared-area value, the channels of that row marked,
    the planned spot-centre extent shaded.  This is the figure that shows what
    the spot model does at the edge channels: with the recorder-width blur the
    field fell by 11-13 % at Ch 0, 4, 5 and 9; with the facility sigma it is
    flat to 2-3 % over Ch 0-9.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x_mm = np.linspace(-100.0, 100.0, 401)
    fig, axes = plt.subplots(
        len(RADFET_ROW_Y_MM), 1, figsize=(9.0, 7.6), sharex=True, constrained_layout=True
    )
    for ax, y_mm in zip(np.atleast_1d(axes), RADFET_ROW_Y_MM):
        for index, (label, fluence, x_edges, y_edges) in enumerate(maps):
            ax.plot(
                x_mm,
                row_profile(fluence, x_edges, y_edges, y_mm, x_mm) / shared_fluence,
                color=f"C{index}",
                label=label,
            )
        ax.axvspan(
            -PLANNED_HALF_EXTENT_MM[0], PLANNED_HALF_EXTENT_MM[0], color="0.9", zorder=0
        )
        for channel, _shielding, channel_x, channel_y in RADEX_RADFETS:
            if abs(channel_y - y_mm) < 0.1:
                ax.axvline(channel_x, color="0.6", linewidth=0.6)
                ax.text(channel_x, 0.55, f"Ch {channel}", ha="center", fontsize=8)
        ax.set_ylim(0.5, 1.1)
        ax.set_ylabel("local / shared-area fluence")
        ax.set_title(f"isocentre y = {y_mm:+.1f} mm")
        ax.grid(alpha=0.3)
    np.atleast_1d(axes)[0].legend(fontsize=8, loc="lower center")
    np.atleast_1d(axes)[-1].set_xlabel("isocentre x [mm]")
    fig.suptitle(
        f"{name}: delivered field along the RadFET rows\n"
        "grey band = planned spot-centre extent (PLD, 78 x 34 spots at 2 mm)"
    )
    fig.savefig(path, dpi=150)
    plt.close(fig)


def spot_lines(
    projection: np.ndarray,
    edges: np.ndarray,
    min_separation_mm: float = 0.6,
    threshold: float = 0.05,
    smoothing_mm: float = 0.15,
) -> tuple[np.ndarray, np.ndarray]:
    """Centre and charge of every spot line in a fine 1-D position projection.

    The projection is smoothed with a Gaussian of ``smoothing_mm``; local maxima
    above ``threshold`` of the global maximum and at least ``min_separation_mm``
    apart (the higher wins) are the lines; every bin is assigned to the nearest
    line (split at the midpoints), which gives each line's charge-weighted
    centre and total charge.  Returns (positions, weights) in ascending position.
    """
    from scipy.ndimage import gaussian_filter1d

    projection = np.asarray(projection, dtype=np.float64)
    centres = (edges[:-1] + edges[1:]) / 2.0
    pixel = float(np.diff(edges).mean())
    smooth = gaussian_filter1d(projection, smoothing_mm / pixel, mode="constant")
    if smooth.max() <= 0:
        return np.empty(0), np.empty(0)
    is_max = (
        (smooth[1:-1] > smooth[:-2])
        & (smooth[1:-1] >= smooth[2:])
        & (smooth[1:-1] > threshold * smooth.max())
    )
    candidates = np.flatnonzero(is_max) + 1
    kept: list[int] = []
    for index in candidates[np.argsort(smooth[candidates])[::-1]]:
        if all(abs(centres[index] - centres[other]) >= min_separation_mm for other in kept):
            kept.append(int(index))
    peaks = np.asarray(sorted(kept), dtype=int)
    bounds = np.concatenate(
        [[edges[0]], (centres[peaks[:-1]] + centres[peaks[1:]]) / 2.0, [edges[-1]]]
    )
    assignment = np.clip(np.searchsorted(bounds, centres, side="right") - 1, 0, len(peaks) - 1)
    weights = np.bincount(assignment, weights=projection, minlength=len(peaks))
    moments = np.bincount(assignment, weights=projection * centres, minlength=len(peaks))
    positions = moments / np.where(weights > 0, weights, np.nan)
    return positions, weights


def lattice_statistics(
    positions: np.ndarray, weights: np.ndarray, planned_pitch_mm: float
) -> dict[str, Any]:
    """Span, pitch and centre of the lines that carry at least 20 % of the median charge."""
    full = weights > 0.2 * np.median(weights) if len(weights) else np.zeros(0, bool)
    kept = positions[full]
    pitch = np.diff(kept)
    median_weight = float(np.median(weights[full])) if full.any() else float("nan")
    return {
        "lines_detected": int(len(positions)),
        "lines_used": int(full.sum()),
        "first_mm": float(kept[0]) if len(kept) else None,
        "last_mm": float(kept[-1]) if len(kept) else None,
        "span_mm": float(kept[-1] - kept[0]) if len(kept) > 1 else None,
        "centre_mm": float((kept[0] + kept[-1]) / 2.0) if len(kept) else None,
        "pitch_median_mm": float(np.median(pitch)) if len(pitch) else None,
        "pitch_mean_mm": float(pitch.mean()) if len(pitch) else None,
        "pitch_sd_mm": float(pitch.std()) if len(pitch) > 1 else None,
        "implied_scale_chamber_to_isocentre": (
            planned_pitch_mm / float(np.median(pitch)) if len(pitch) else None
        ),
        "edge_line_charge_over_median": (
            [float(weights[full][0] / median_weight), float(weights[full][-1] / median_weight)]
            if full.any()
            else None
        ),
    }


def write_spot_lattice_products(
    aggregates_root: Path, output_dir: Path
) -> dict[str, Any] | None:
    """Delivered spot lattice of the two energies in the recorder frame.

    Reads the fine projections of ``RadEx_64MeV`` and ``Radex_85MeV``
    (``heatmap_positions.npz``, parsed since 2026-10-26), detects the spot
    lines along x and y, and writes ``spot_lattice_64_vs_85.png`` (charge per
    line over its median against line position) and
    ``spot_lattice_summary.json`` (lines, span, pitch, centre, implied scale,
    and the differences between the energies, also projected to the isocentre
    with the registered transform).  Returns None when the projections are not
    there yet (records parsed before 2026-10-26).
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    energies = {"RadEx_64MeV": 64.0, "Radex_85MeV": 85.0}
    projections: dict[str, dict[str, Any]] = {}
    for component, energy_mev in energies.items():
        path = aggregates_root / component / "heatmap_positions.npz"
        if not path.exists():
            return None
        with np.load(path) as data:
            if "projection_x_fine_protons" not in data.files:
                return None
            projections[component] = {
                "energy_mev": energy_mev,
                "x": data["projection_x_fine_protons"].astype(np.float64),
                "y": data["projection_y_fine_protons"].astype(np.float64),
                "edges": data["projection_fine_edges_mm"].astype(np.float64),
                "source": str(path.resolve()),
            }
    summary: dict[str, Any] = {
        "method": (
            "Charge-weighted spot positions of every part record, projected on the "
            "recorder x and y axes in 0.05 mm bins (uppsala_analysis "
            "FINE_PROJECTION_PIXEL_MM), summed per energy in protons. Spot lines: local "
            "maxima of the 0.15 mm-smoothed projection above 5 % of its maximum and at "
            "least 0.6 mm apart; each bin belongs to the nearest line. Span, pitch and "
            "centre use the lines carrying at least 20 % of the median line charge."
        ),
        "planned_pitch_mm": PLD_SPOT_PITCH_MM,
        "registration": RADEX_CHAMBER_TO_ISOCENTRE,
        "energies": {},
    }
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.2), constrained_layout=True)
    for index, (component, item) in enumerate(projections.items()):
        per_axis: dict[str, Any] = {"source": item["source"]}
        for axis_index, axis in enumerate(("x", "y")):
            positions, weights = spot_lines(item[axis], item["edges"])
            statistics = lattice_statistics(positions, weights, PLD_SPOT_PITCH_MM)
            per_axis[axis] = statistics
            full = weights > 0.2 * np.median(weights) if len(weights) else np.zeros(0, bool)
            axes[axis_index].plot(
                positions[full],
                weights[full] / np.median(weights[full]),
                marker="o",
                markersize=3,
                linewidth=0.8,
                color=f"C{index}",
                label=(
                    f"{item['energy_mev']:g} MeV: {statistics['lines_used']} lines, "
                    f"span {statistics['span_mm']:.2f} mm, pitch "
                    f"{statistics['pitch_median_mm']:.3f} mm"
                ),
            )
        summary["energies"][component] = per_axis
    for axis_index, axis in enumerate(("x", "y")):
        axes[axis_index].set_xlabel(f"recorder (IC2/3) {axis} position [mm]")
        axes[axis_index].set_ylabel("charge per spot line / median")
        axes[axis_index].grid(alpha=0.3)
        axes[axis_index].legend(fontsize=8)
    components = list(projections)
    differences: dict[str, Any] = {}
    for axis, scale_key in (("x", "scale_x_chamber_to_isocentre"), ("y", "scale_y_chamber_to_isocentre")):
        first = summary["energies"][components[0]][axis]
        second = summary["energies"][components[1]][axis]
        if first["span_mm"] and second["span_mm"]:
            centre_difference = second["centre_mm"] - first["centre_mm"]
            differences[axis] = {
                "span_relative_difference": second["span_mm"] / first["span_mm"] - 1.0,
                "centre_difference_chamber_mm": centre_difference,
                "centre_difference_isocentre_mm": (
                    centre_difference * RADEX_CHAMBER_TO_ISOCENTRE[scale_key]
                ),
                "implied_scales": [
                    first["implied_scale_chamber_to_isocentre"],
                    second["implied_scale_chamber_to_isocentre"],
                ],
                "registered_scale": RADEX_CHAMBER_TO_ISOCENTRE[scale_key],
            }
    summary["85_minus_64"] = differences
    fig.suptitle(
        "Delivered spot lattice in the recorder frame, all RadEx spills "
        "(from the fine spot-position projections)"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / "spot_lattice_64_vs_85.png", dpi=150)
    plt.close(fig)
    write_summary(output_dir / "spot_lattice_summary.json", summary)
    return summary


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


RECORDER_METHOD_NOTE = (
    "The chamber-to-isocentre coordinate scale is 1.480648 in X "
    "and 1.328801 in Y. Fluence density includes the inverse area "
    "Jacobian 0.508263, which preserves the integrated proton count. "
    "Spot model: recorder X_WIDTH/Y_WIDTH applied as a chamber-plane sigma "
    "and magnified with the map (superseded 2026-10-23; kept for comparison)."
)


def write_markdown(
    path: Path,
    title: str,
    rows: list[dict[str, Any]],
    method_note: str = RECORDER_METHOD_NOTE,
    shared_fluence: float | None = None,
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
    ]
    if shared_fluence:
        lines.extend(
            [
                "| Ch. | Shielding | GDML X [mm] | GDML Y [mm] | Delivered fluence [protons/cm2] | Local / shared |",
                "| ---: | --- | ---: | ---: | ---: | ---: |",
            ]
        )
    else:
        lines.extend(
            [
                "| Ch. | Shielding | GDML X [mm] | GDML Y [mm] | Delivered fluence [protons/cm2] |",
                "| ---: | --- | ---: | ---: | ---: |",
            ]
        )
    for row in rows:
        line = (
            f"| {row['channel']} | {row['shielding']} | "
            f"{row['gdml_x_mm']:+g} | {row['gdml_y_mm']:+g} | "
            f"{row['delivered_fluence_protons_cm2']:.6e} |"
        )
        if shared_fluence:
            line += f" {row['delivered_fluence_protons_cm2'] / shared_fluence:.4f} |"
        lines.append(line)
    lines.append("")
    if shared_fluence:
        lines.extend(
            [
                (
                    f"Shared value: delivered protons / {RADEX_INSTRUMENT_AREA_CM2:.3f} cm2 "
                    f"= {shared_fluence:.6e} protons/cm2."
                ),
                "",
            ]
        )
    lines.extend(
        [
            method_note,
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


def write_summary(path: Path, summary: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")


def save_isocentre_npz(
    path: Path,
    fluence: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    **extra: Any,
) -> None:
    np.savez_compressed(
        path,
        fluence_protons_cm2=fluence.astype(np.float32),
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
        **extra,
    )


def process_map(
    input_path: Path,
    output_dir: Path,
    title: str,
    grid_min_mm: float,
    grid_max_mm: float,
) -> dict[str, Any]:
    """Recorder-blur isocentre products (method before 2026-10-23)."""
    values, registered, x_edges, y_edges = recorder_blur_isocentre_map(
        input_path, grid_min_mm, grid_max_mm
    )
    pixel_area_cm2 = float(np.diff(x_edges).mean() * np.diff(y_edges).mean() / 100.0)
    input_integral = float(values.sum() * pixel_area_cm2)
    output_integral = float(registered.sum() * pixel_area_cm2)
    rows = location_rows(registered, x_edges, y_edges)
    output_dir.mkdir(parents=True, exist_ok=True)
    save_isocentre_npz(output_dir / "heatmap_isocentre.npz", registered, x_edges, y_edges)
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
        "spot_model": "recorder",
        "spot_model_note": (
            "Superseded 2026-10-23: recorder X_WIDTH/Y_WIDTH applied as a "
            "chamber-plane sigma before the chamber-to-isocentre magnification."
        ),
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
    write_summary(output_dir / "isocentre_summary.json", summary)
    return summary


def process_aggregate_facility(
    aggregates_root: Path,
    name: str,
    output_dir: Path,
    grid_min_mm: float,
    grid_max_mm: float,
) -> dict[str, Any]:
    """Facility-sigma isocentre products for one RadEx aggregate."""
    result = facility_sigma_isocentre_map(
        aggregates_root, name, grid_min_mm, grid_max_mm
    )
    fluence = result["fluence"]
    x_edges = result["x_edges"]
    y_edges = result["y_edges"]
    pixel_area_cm2 = float(np.diff(x_edges).mean() * np.diff(y_edges).mean() / 100.0)
    components = result["components"]
    input_integral = sum(item["chamber_integrated_protons"] for item in components)
    output_integral = float(fluence.sum() * pixel_area_cm2)

    aggregate_summary_path = aggregates_root / name / "summary.json"
    aggregate_summary = json.loads(aggregate_summary_path.read_text(encoding="utf-8"))
    shared_fluence = float(aggregate_summary["radex_mean_fluence_raw_protons_cm2"])
    spatialized = float(aggregate_summary["spatialized_protons"])
    delivered = float(aggregate_summary["delivered_protons_raw"])

    own_positions = aggregates_root / name / "heatmap_positions.npz"
    own_integral = None
    if own_positions.exists():
        own_values, own_x, own_y = load_heatmap(
            own_positions, ("fluence_protons_cm2",), grid_min_mm, grid_max_mm
        )
        own_integral = float(
            own_values.sum() * np.diff(own_x).mean() * np.diff(own_y).mean() / 100.0
        )

    rows = location_rows(fluence, x_edges, y_edges)
    local = np.asarray([row["delivered_fluence_protons_cm2"] for row in rows])
    cross_local = sample_radfets(result["cross_check_fluence"], x_edges, y_edges)

    recorder_comparison = None
    recorder_path = aggregates_root / name / "heatmap.npz"
    if recorder_path.exists():
        _values, recorder_map, rx, ry = recorder_blur_isocentre_map(
            recorder_path, grid_min_mm, grid_max_mm
        )
        recorder_local = sample_radfets(recorder_map, rx, ry)
        recorder_comparison = {
            "source_heatmap": str(recorder_path.resolve()),
            "radfet_fluence_protons_cm2": recorder_local.tolist(),
            "radfet_local_to_shared": (recorder_local / shared_fluence).tolist(),
            "facility_over_recorder": (local / recorder_local).tolist(),
        }

    sigma_text = "; ".join(
        f"{item['nominal_energy_mev']:g} MeV: sigma {item['spot_sigma_x_mm']:.3f} x "
        f"{item['spot_sigma_y_mm']:.3f} mm"
        for item in components
    )
    method_note = (
        "Spot model (since 2026-10-23): the unblurred charge-weighted spot-position "
        "map is transformed to the isocentre (scale 1.480648 in X and 1.328801 in Y, "
        "forward deposit, which carries the inverse area Jacobian 0.508263 and "
        "conserves the integrated proton count) and then blurred with the facility "
        f"isocentre spot sigma ({sigma_text}; Skandion email 2026-08-24, "
        "pld_analysis.SPOT_SIGMA_MM interpolated in energy)."
        + (" Each energy is blurred with its own sigma before summation." if len(components) > 1 else "")
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    save_isocentre_npz(
        output_dir / "heatmap_isocentre.npz",
        fluence,
        x_edges,
        y_edges,
        spot_sigma_energies_mev=np.asarray(
            [item["nominal_energy_mev"] for item in components], dtype=float
        ),
        spot_sigma_x_mm=np.asarray(
            [item["spot_sigma_x_mm"] for item in components], dtype=float
        ),
        spot_sigma_y_mm=np.asarray(
            [item["spot_sigma_y_mm"] for item in components], dtype=float
        ),
    )
    write_csv(output_dir / "radfet_location_fluence.csv", rows)
    write_markdown(
        output_dir / "radfet_location_fluence.md",
        name,
        rows,
        method_note=method_note,
        shared_fluence=shared_fluence,
    )
    save_heatmap(
        output_dir / "fluence_isocentre.png",
        f"{name} — delivered fluence at isocentre\nfacility spot sigma: {sigma_text}",
        fluence,
        x_edges,
        y_edges,
        rows,
    )
    profile_maps = [(f"facility spot sigma ({sigma_text})", fluence, x_edges, y_edges)]
    if recorder_comparison is not None:
        profile_maps.append(
            ("recorder-width blur, magnified (superseded 2026-10-23)", recorder_map, rx, ry)
        )
    save_row_profiles(
        output_dir / "row_profiles_spot_models.png", name, profile_maps, shared_fluence
    )
    summary = {
        "title": name,
        "spot_model": "facility",
        "spot_model_note": (
            "Unblurred charge-weighted spot positions (heatmap_positions.npz) "
            "transformed to the isocentre by forward cloud-in-cell deposit with the "
            "registered scales/offsets (equivalent to the inverse-Jacobian density "
            "correction), then blurred with the facility isocentre spot sigma "
            "(scipy.ndimage.gaussian_filter, truncate 5) and sampled bilinearly at "
            "the RadFET GDML coordinates. Replaces the recorder-blur method "
            "(recorder X_WIDTH/Y_WIDTH as chamber-plane sigma, magnified with the map)."
        ),
        "spot_sigma_source": SPOT_SIGMA_SOURCE,
        "spot_sigma": [
            {
                "aggregate": item["aggregate"],
                "nominal_energy_mev": item["nominal_energy_mev"],
                "sigma_x_mm": item["spot_sigma_x_mm"],
                "sigma_y_mm": item["spot_sigma_y_mm"],
            }
            for item in components
        ],
        "components": components,
        "source_heatmaps": [item["source_heatmap"] for item in components],
        "transform": RADEX_CHAMBER_TO_ISOCENTRE,
        "input_integrated_protons": input_integral,
        "isocentre_integrated_protons": output_integral,
        "integral_ratio": output_integral / input_integral,
        "aggregate_positions_integrated_protons": own_integral,
        "aggregate_spatialized_protons": spatialized,
        "aggregate_delivered_protons_raw": delivered,
        "isocentre_to_spatialized_protons_ratio": output_integral / spatialized,
        "shared_area_fluence_protons_cm2": shared_fluence,
        "radfet_local_to_shared": (local / shared_fluence).tolist(),
        "cross_check_chamber_blur_then_resample": {
            "method": (
                "Blur the chamber-plane positions with sigma/scale, then resample "
                "with registered_full_grid x Jacobian (the recorder-path transform)."
            ),
            "radfet_fluence_protons_cm2": cross_local.tolist(),
            "max_abs_relative_difference": float(np.max(np.abs(cross_local / local - 1.0))),
        },
        "recorder_blur_comparison": recorder_comparison,
        "row_profile_figure": "row_profiles_spot_models.png",
        "radfet_coordinate_assumption": (
            "Nominal GDML coordinates centred and aligned with machine X/Y; "
            "no additional DUT pose transform"
        ),
    }
    write_summary(output_dir / "isocentre_summary.json", summary)
    return summary


def parse_map(text: str) -> tuple[str, Path]:
    if "=" not in text:
        raise argparse.ArgumentTypeError("map must be TITLE=PATH")
    title, path = text.split("=", 1)
    return title, Path(path)


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--aggregates",
        type=Path,
        help="results/final/aggregates directory; processes --names there",
    )
    parser.add_argument(
        "--names",
        nargs="+",
        default=list(RADEX_AGGREGATE_COMPONENTS),
        help="RadEx aggregates to process (default: all three)",
    )
    parser.add_argument(
        "--spot-model",
        choices=("facility", "recorder"),
        default="facility",
        help=(
            "facility (default): facility isocentre sigma applied after the "
            "transform; recorder: superseded recorder-width blur, for comparison"
        ),
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        help="write products to OUTPUT_ROOT/<name>/ instead of the aggregate directory",
    )
    parser.add_argument(
        "--map",
        action="append",
        type=parse_map,
        help="TITLE=PATH of a recorder heatmap.npz (legacy; --spot-model recorder only)",
    )
    parser.add_argument("--grid-min-mm", type=float, default=-250.0)
    parser.add_argument("--grid-max-mm", type=float, default=250.0)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = make_parser()
    args = parser.parse_args(argv)
    if args.map:
        if args.spot_model != "recorder":
            parser.error("--map processes recorder heatmaps; add --spot-model recorder")
        for title, input_path in args.map:
            output_dir = (
                args.output_root / title if args.output_root else input_path.parent
            )
            summary = process_map(
                input_path, output_dir, title, args.grid_min_mm, args.grid_max_mm
            )
            print(
                f"{title}: wrote recorder-blur isocentre products; "
                f"integral ratio={summary['integral_ratio']:.8f}"
            )
        return 0
    if not args.aggregates:
        parser.error("give --aggregates (or --map with --spot-model recorder)")
    for name in args.names:
        output_dir = args.output_root / name if args.output_root else args.aggregates / name
        if args.spot_model == "recorder":
            summary = process_map(
                args.aggregates / name / "heatmap.npz",
                output_dir,
                name,
                args.grid_min_mm,
                args.grid_max_mm,
            )
            print(
                f"{name}: wrote recorder-blur isocentre products; "
                f"integral ratio={summary['integral_ratio']:.8f}"
            )
            continue
        summary = process_aggregate_facility(
            args.aggregates, name, output_dir, args.grid_min_mm, args.grid_max_mm
        )
        sigmas = ", ".join(
            f"{item['nominal_energy_mev']:g} MeV {item['sigma_x_mm']:.3f} mm"
            for item in summary["spot_sigma"]
        )
        local = summary["radfet_local_to_shared"]
        print(
            f"{name}: facility sigma ({sigmas}); integral ratio="
            f"{summary['integral_ratio']:.8f}; local/shared Ch0-11 "
            + " ".join(f"{value:.4f}" for value in local)
        )
    if args.spot_model == "facility" and "RadEx_Total" in args.names:
        output_dir = (
            args.output_root / "RadEx_Total" if args.output_root else args.aggregates / "RadEx_Total"
        )
        lattice = write_spot_lattice_products(args.aggregates, output_dir)
        if lattice is None:
            print(
                "RadEx_Total: no spot-lattice products (heatmap_positions.npz has no fine "
                "projections; re-parse the RadEx part records and reduce)"
            )
        else:
            for axis in ("x", "y"):
                item = lattice["85_minus_64"].get(axis)
                if item:
                    print(
                        f"RadEx_Total: spot lattice {axis}: span differs by "
                        f"{100 * item['span_relative_difference']:+.3f} %, centres by "
                        f"{item['centre_difference_chamber_mm']:+.3f} mm (chamber); implied "
                        f"scales {item['implied_scales'][0]:.4f} / {item['implied_scales'][1]:.4f} "
                        f"against the registered {item['registered_scale']:.4f}"
                    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
