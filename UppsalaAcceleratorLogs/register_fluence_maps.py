#!/usr/bin/env python3
"""Register an ionisation-chamber fluence map to a PLD isocentre map.

The transform is

    x_isocentre = scale_x * x_chamber + offset_x
    y_isocentre = scale_y * y_chamber + offset_y

Independent scale factors describe the chamber-to-isocentre geometric
magnification.  Translation and intensity are nuisance parameters.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np


def map_array(data: np.lib.npyio.NpzFile, preferred: tuple[str, ...]) -> np.ndarray:
    for key in preferred:
        if key in data.files:
            values = data[key].astype(np.float64)
            if values.ndim != 2:
                raise ValueError(f"{key} is not a two-dimensional heatmap")
            return values
    raise ValueError(f"None of {preferred} found; available arrays: {data.files}")


def edges_for(
    data: np.lib.npyio.NpzFile,
    shape: tuple[int, int],
    fallback_min_mm: float,
    fallback_max_mm: float,
) -> tuple[np.ndarray, np.ndarray]:
    if "x_edges_mm" in data.files and "y_edges_mm" in data.files:
        x_edges = data["x_edges_mm"].astype(np.float64)
        y_edges = data["y_edges_mm"].astype(np.float64)
    else:
        x_edges = np.linspace(fallback_min_mm, fallback_max_mm, shape[1] + 1)
        y_edges = np.linspace(fallback_min_mm, fallback_max_mm, shape[0] + 1)
    if len(x_edges) != shape[1] + 1 or len(y_edges) != shape[0] + 1:
        raise ValueError("Heatmap dimensions do not match coordinate edges")
    if not (np.all(np.diff(x_edges) > 0) and np.all(np.diff(y_edges) > 0)):
        raise ValueError("Heatmap coordinate edges must increase monotonically")
    return x_edges, y_edges


def normalized(values: np.ndarray) -> np.ndarray:
    clean = np.where(np.isfinite(values) & (values > 0), values, 0.0)
    peak = float(clean.max()) if clean.size else 0.0
    if peak <= 0:
        raise ValueError("Heatmap has no positive finite values")
    return clean / peak


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    aa = a - a.mean()
    bb = b - b.mean()
    denominator = math.sqrt(float(np.dot(aa, aa) * np.dot(bb, bb)))
    return float(np.dot(aa, bb) / denominator) if denominator else float("nan")


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    denominator = math.sqrt(float(np.dot(a, a) * np.dot(b, b)))
    return float(np.dot(a, b) / denominator) if denominator else float("nan")


def query_grid(
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    downsample: int,
    x_half_window_mm: float,
    y_half_window_mm: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = (x_edges[:-1] + x_edges[1:]) / 2.0
    y = (y_edges[:-1] + y_edges[1:]) / 2.0
    query_x, query_y = np.meshgrid(x[::downsample], y[::downsample])
    mask = (
        (np.abs(query_x) <= x_half_window_mm)
        & (np.abs(query_y) <= y_half_window_mm)
    )
    return query_x[mask], query_y[mask], mask


def sample_registered(
    observed: np.ndarray,
    observed_x_edges: np.ndarray,
    observed_y_edges: np.ndarray,
    query_x: np.ndarray,
    query_y: np.ndarray,
    parameters: np.ndarray | tuple[float, float, float, float],
) -> np.ndarray:
    try:
        from scipy.ndimage import map_coordinates
    except ImportError as exc:
        raise RuntimeError(
            "SciPy is required; on Triton load scicomp-python-env/2025.2"
        ) from exc
    scale_x, scale_y, offset_x, offset_y = parameters
    source_x = (query_x - offset_x) / scale_x
    source_y = (query_y - offset_y) / scale_y
    pixel_x = float(np.diff(observed_x_edges).mean())
    pixel_y = float(np.diff(observed_y_edges).mean())
    coordinates = np.asarray(
        [
            (source_y - observed_y_edges[0]) / pixel_y - 0.5,
            (source_x - observed_x_edges[0]) / pixel_x - 0.5,
        ]
    )
    return map_coordinates(
        observed,
        coordinates,
        order=1,
        mode="constant",
        cval=0.0,
    )


def fit_registration(
    planned: np.ndarray,
    observed: np.ndarray,
    planned_x_edges: np.ndarray,
    planned_y_edges: np.ndarray,
    observed_x_edges: np.ndarray,
    observed_y_edges: np.ndarray,
    x_half_window_mm: float = 110.0,
    y_half_window_mm: float = 65.0,
    downsample: int = 2,
    initial: tuple[float, float, float, float] | None = None,
    global_search: bool = True,
    seed: int = 20260916,
) -> dict[str, float]:
    try:
        from scipy.optimize import differential_evolution, minimize
    except ImportError as exc:
        raise RuntimeError(
            "SciPy is required; on Triton load scicomp-python-env/2025.2"
        ) from exc
    planned_norm = normalized(planned)
    observed_norm = normalized(observed)
    query_x, query_y, mask = query_grid(
        planned_x_edges,
        planned_y_edges,
        downsample,
        x_half_window_mm,
        y_half_window_mm,
    )
    planned_sample = planned_norm[::downsample, ::downsample][mask]
    bounds = ((1.0, 2.0), (1.0, 1.8), (-15.0, 15.0), (-15.0, 15.0))

    def objective(parameters: np.ndarray) -> float:
        observed_sample = sample_registered(
            observed_norm,
            observed_x_edges,
            observed_y_edges,
            query_x,
            query_y,
            parameters,
        )
        return -pearson(planned_sample, observed_sample)

    if global_search:
        result = differential_evolution(
            objective,
            bounds=bounds,
            seed=seed,
            tol=1e-7,
            maxiter=80,
            popsize=10,
            polish=False,
            workers=1,
        )
        start = result.x
    elif initial is not None:
        start = np.asarray(initial, dtype=float)
    else:
        start = np.asarray((1.48, 1.33, 0.6, 2.3), dtype=float)
    polished = minimize(
        objective,
        start,
        method="Powell",
        bounds=bounds,
        options={"xtol": 1e-7, "ftol": 1e-10, "maxiter": 2000},
    )
    parameters = polished.x
    observed_sample = sample_registered(
        observed_norm,
        observed_x_edges,
        observed_y_edges,
        query_x,
        query_y,
        parameters,
    )
    gain = float(np.dot(planned_sample, observed_sample) / np.dot(observed_sample, observed_sample))
    area_scale = float(parameters[0] * parameters[1])
    return {
        "scale_x_chamber_to_isocentre": float(parameters[0]),
        "scale_y_chamber_to_isocentre": float(parameters[1]),
        "offset_x_isocentre_mm": float(parameters[2]),
        "offset_y_isocentre_mm": float(parameters[3]),
        "pearson_correlation": pearson(planned_sample, observed_sample),
        "cosine_similarity": cosine(planned_sample, observed_sample),
        "normalized_intensity_gain": gain,
        "normalized_rmse": float(
            np.sqrt(np.mean((planned_sample - gain * observed_sample) ** 2))
        ),
        "area_scale_chamber_to_isocentre": area_scale,
        "fluence_density_jacobian_chamber_to_isocentre": 1.0 / area_scale,
        "x_half_window_mm": x_half_window_mm,
        "y_half_window_mm": y_half_window_mm,
        "sampling_mm": float(np.diff(planned_x_edges).mean()) * downsample,
    }


def fit_translation_only(
    planned: np.ndarray,
    observed: np.ndarray,
    planned_x_edges: np.ndarray,
    planned_y_edges: np.ndarray,
    observed_x_edges: np.ndarray,
    observed_y_edges: np.ndarray,
    x_half_window_mm: float,
    y_half_window_mm: float,
    downsample: int,
) -> dict[str, float]:
    from scipy.optimize import minimize

    p = normalized(planned)
    o = normalized(observed)
    qx, qy, mask = query_grid(
        planned_x_edges,
        planned_y_edges,
        downsample,
        x_half_window_mm,
        y_half_window_mm,
    )
    ps = p[::downsample, ::downsample][mask]

    def objective(offsets: np.ndarray) -> float:
        values = sample_registered(
            o,
            observed_x_edges,
            observed_y_edges,
            qx,
            qy,
            (1.0, 1.0, offsets[0], offsets[1]),
        )
        return -pearson(ps, values)

    result = minimize(
        objective,
        (0.0, 0.0),
        method="Powell",
        bounds=((-15, 15), (-15, 15)),
    )
    return {
        "offset_x_isocentre_mm": float(result.x[0]),
        "offset_y_isocentre_mm": float(result.x[1]),
        "pearson_correlation": float(-result.fun),
    }


def load_heatmap(
    path: Path,
    preferred_keys: tuple[str, ...],
    fallback_min_mm: float,
    fallback_max_mm: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with np.load(path) as data:
        values = map_array(data, preferred_keys)
        x_edges, y_edges = edges_for(
            data,
            values.shape,
            fallback_min_mm,
            fallback_max_mm,
        )
    return values, x_edges, y_edges


def registered_full_grid(
    observed: np.ndarray,
    observed_x_edges: np.ndarray,
    observed_y_edges: np.ndarray,
    planned_x_edges: np.ndarray,
    planned_y_edges: np.ndarray,
    fit: dict[str, float],
    normalize_input: bool = True,
) -> np.ndarray:
    x = (planned_x_edges[:-1] + planned_x_edges[1:]) / 2.0
    y = (planned_y_edges[:-1] + planned_y_edges[1:]) / 2.0
    query_x, query_y = np.meshgrid(x, y)
    parameters = (
        fit["scale_x_chamber_to_isocentre"],
        fit["scale_y_chamber_to_isocentre"],
        fit["offset_x_isocentre_mm"],
        fit["offset_y_isocentre_mm"],
    )
    values = normalized(observed) if normalize_input else observed
    return sample_registered(
        values,
        observed_x_edges,
        observed_y_edges,
        query_x,
        query_y,
        parameters,
    )


def save_diagnostic(
    path: Path,
    planned: np.ndarray,
    registered: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    fit: dict[str, float],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p = normalized(planned)
    gain = fit["normalized_intensity_gain"]
    registered_scaled = gain * registered
    difference = registered_scaled - p
    extent = [x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    for ax, values, title, cmap, limits in (
        (axes[0], p, "PLD planned at isocentre", "magma", (0, 1)),
        (axes[1], registered_scaled, "Chamber map scaled to isocentre", "magma", (0, 1)),
        (axes[2], difference, "Registered chamber - PLD", "coolwarm", (-0.25, 0.25)),
    ):
        image = ax.imshow(
            values,
            origin="lower",
            extent=extent,
            cmap=cmap,
            vmin=limits[0],
            vmax=limits[1],
        )
        ax.set_xlim(-110, 110)
        ax.set_ylim(-65, 65)
        ax.set_xlabel("Isocentre X [mm]")
        ax.set_ylabel("Isocentre Y [mm]")
        ax.set_title(title)
        fig.colorbar(image, ax=ax, shrink=0.82)
    fig.suptitle(
        f"sx={fit['scale_x_chamber_to_isocentre']:.4f}, "
        f"sy={fit['scale_y_chamber_to_isocentre']:.4f}, "
        f"r={fit['pearson_correlation']:.5f}"
    )
    fig.savefig(path, dpi=180)
    plt.close(fig)


def parse_validation(text: str) -> tuple[str, Path]:
    if "=" not in text:
        raise argparse.ArgumentTypeError("validation map must be NAME=PATH")
    name, path = text.split("=", 1)
    return name, Path(path)


def run(args: argparse.Namespace) -> None:
    planned, planned_x, planned_y = load_heatmap(
        Path(args.planned),
        ("planned_fluence_protons_cm2", "fluence_protons_cm2"),
        args.grid_min_mm,
        args.grid_max_mm,
    )
    observed, observed_x, observed_y = load_heatmap(
        Path(args.observed),
        ("fluence_protons_cm2", "planned_fluence_protons_cm2"),
        args.grid_min_mm,
        args.grid_max_mm,
    )
    fit = fit_registration(
        planned,
        observed,
        planned_x,
        planned_y,
        observed_x,
        observed_y,
        args.x_half_window_mm,
        args.y_half_window_mm,
        args.downsample,
    )
    translation_only = fit_translation_only(
        planned,
        observed,
        planned_x,
        planned_y,
        observed_x,
        observed_y,
        args.x_half_window_mm,
        args.y_half_window_mm,
        args.downsample,
    )
    sensitivity = []
    initial = (
        fit["scale_x_chamber_to_isocentre"],
        fit["scale_y_chamber_to_isocentre"],
        fit["offset_x_isocentre_mm"],
        fit["offset_y_isocentre_mm"],
    )
    for x_half, y_half in args.sensitivity_window:
        sensitivity.append(
            fit_registration(
                planned,
                observed,
                planned_x,
                planned_y,
                observed_x,
                observed_y,
                x_half,
                y_half,
                args.downsample,
                initial=initial,
                global_search=False,
            )
        )
    validations: dict[str, Any] = {}
    for name, validation_path in args.validation_map:
        values, x_edges, y_edges = load_heatmap(
            validation_path,
            ("fluence_protons_cm2", "planned_fluence_protons_cm2"),
            args.grid_min_mm,
            args.grid_max_mm,
        )
        validations[name] = fit_registration(
            planned,
            values,
            planned_x,
            planned_y,
            x_edges,
            y_edges,
            args.x_half_window_mm,
            args.y_half_window_mm,
            args.downsample,
            initial=initial,
            global_search=False,
        )
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    result = {
        "transform_definition": (
            "x_isocentre = scale_x * x_chamber + offset_x; "
            "y_isocentre = scale_y * y_chamber + offset_y"
        ),
        "planned_map": str(Path(args.planned).resolve()),
        "observed_map": str(Path(args.observed).resolve()),
        "fit": fit,
        "translation_only_reference": translation_only,
        "sensitivity_windows": sensitivity,
        "validation_maps": validations,
    }
    with (output_dir / "registration_summary.json").open(
        "w", encoding="utf-8", newline="\n"
    ) as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    registered = registered_full_grid(
        observed,
        observed_x,
        observed_y,
        planned_x,
        planned_y,
        fit,
    )
    registered_absolute = registered_full_grid(
        observed,
        observed_x,
        observed_y,
        planned_x,
        planned_y,
        fit,
        normalize_input=False,
    ) / fit["area_scale_chamber_to_isocentre"]
    np.savez_compressed(
        output_dir / "registered_heatmaps.npz",
        planned_normalized=normalized(planned).astype(np.float32),
        chamber_registered_normalized=registered.astype(np.float32),
        chamber_registered_fluence_conserving=registered_absolute.astype(np.float32),
        x_edges_mm=planned_x,
        y_edges_mm=planned_y,
    )
    save_diagnostic(
        output_dir / "registration_diagnostic.png",
        planned,
        registered,
        planned_x,
        planned_y,
        fit,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


def window(text: str) -> tuple[float, float]:
    try:
        x, y = text.split(",", 1)
        return float(x), float(y)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("window must be X_HALF,Y_HALF") from exc


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--planned", required=True)
    parser.add_argument("--observed", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--validation-map", action="append", type=parse_validation, default=[])
    parser.add_argument("--sensitivity-window", action="append", type=window, default=[])
    parser.add_argument("--x-half-window-mm", type=float, default=110.0)
    parser.add_argument("--y-half-window-mm", type=float, default=65.0)
    parser.add_argument("--downsample", type=int, default=2)
    parser.add_argument("--grid-min-mm", type=float, default=-250.0)
    parser.add_argument("--grid-max-mm", type=float, default=250.0)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    run(make_parser().parse_args(argv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
