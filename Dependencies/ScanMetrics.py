"""Derive depth-dose and lateral-profile metrics from measured or simulated scans."""

from __future__ import annotations

import numpy as np

PROFILE_CURVE_TYPES = ("INPLANE_PROFILE", "CROSSPLANE_PROFILE", "DIAGONAL_PROFILE")
DEPTH_CURVE_TYPES = ("PDD", "DEPTH_DOSE")

PROFILE_LEVELS = (90.0, 80.0, 50.0, 20.0)
DEPTH_LEVELS = (90.0, 80.0, 50.0)
FLATNESS_LEVEL = 90.0
SYMMETRY_LEVEL = 80.0


def cellBounds(depths):
    """
    Build scoring-cell boundaries for a depth grid whose first cell starts at the surface.

    Args:
        depths (np.ndarray): Cell centre depths in mm, ascending.

    Returns:
        np.ndarray: ``len(depths) + 1`` boundaries in mm.
    """
    depths = np.asarray(depths, dtype=float)
    if len(depths) < 2:
        raise ValueError("Need at least two depths to derive cell boundaries")
    boundaries = np.empty(len(depths) + 1)
    boundaries[0] = 0.0
    boundaries[1:-1] = (depths[:-1] + depths[1:]) / 2
    boundaries[-1] = depths[-1] + (depths[-1] - depths[-2]) / 2
    return boundaries


def levelCrossing(position, normalized, level, side):
    """
    Interpolate where a normalised lateral profile crosses ``level`` on one side of the axis.

    Args:
        position (np.ndarray): Scan positions in mm, ascending.
        normalized (np.ndarray): Profile normalised to the central-axis value, in percent.
        level (float): Crossing level in percent.
        side (str): 'left' for negative positions, 'right' for positive ones.

    Returns:
        float: Crossing position in mm, or NaN when the scan does not reach the level.
    """
    if side == "left":
        positions, values = position[position <= 0], normalized[position <= 0]
    elif side == "right":
        positions = position[position >= 0][::-1]
        values = normalized[position >= 0][::-1]
    else:
        raise ValueError(f"Unknown side {side!r}")

    for index in range(1, len(positions)):
        first, second = values[index - 1], values[index]
        if (first - level) * (second - level) <= 0 and first != second:
            fraction = (level - first) / (second - first)
            return positions[index - 1] + fraction * (positions[index] - positions[index - 1])
    return np.nan


def fallingCrossing(depths, values, level):
    """
    Interpolate where a depth-dose curve first falls through ``level`` beyond its peak.

    Args:
        depths (np.ndarray): Depths in mm, ascending.
        values (np.ndarray): Dose normalised to its maximum, in percent.
        level (float): Crossing level in percent.

    Returns:
        float: Depth in mm, or NaN when the curve does not fall through the level.
    """
    peak = int(np.argmax(values))
    for right in range(peak + 1, len(values)):
        if values[right] <= level < values[right - 1]:
            left = right - 1
            fraction = (level - values[left]) / (values[right] - values[left])
            return depths[left] + fraction * (depths[right] - depths[left])
    return np.nan


def _ordered(scan):
    """Return the scan positions and values sorted by ascending position."""
    position, value = np.asarray(scan["position"]), np.asarray(scan["value"])
    order = np.argsort(position)
    return position[order], value[order]


def profileMetrics(scan, levels=PROFILE_LEVELS):
    """
    Derive field width, flatness and symmetry metrics from one lateral profile scan.

    The profile is normalised to the interpolated central-axis value, so a scan
    whose maximum is off axis still reports the central axis as 100 percent.
    Levels the scan does not reach are reported as NaN rather than extrapolated.
    Flatness is the half-span inside the 90 percent width; symmetry is the largest
    left-right difference inside the 80 percent width, in percentage points, and is
    evaluated only where the scan covers both signs of the off-axis position.

    Args:
        scan (dict): One scan from ``Read.ReadMCC.readMCC``.
        levels (tuple[float, ...]): Relative levels in percent to report widths at.

    Returns:
        dict: Normalised profile plus width, flatness and symmetry metrics.
    """
    position, value = _ordered(scan)
    if position.min() > 0 or position.max() < 0:
        raise ValueError(f"{scan['curve_type']}: scan does not cross the central axis")

    centre = float(np.interp(0.0, position, value))
    if not centre > 0:
        raise ValueError(f"{scan['curve_type']}: central-axis value is not positive")
    normalized = 100 * value / centre

    crossings = {level: {side: levelCrossing(position, normalized, level, side)
                         for side in ("left", "right")}
                 for level in set(levels) | {FLATNESS_LEVEL, SYMMETRY_LEVEL}}

    flat = crossings[FLATNESS_LEVEL]
    inside = (position >= flat["left"]) & (position <= flat["right"])
    span = float(normalized[inside].max() - normalized[inside].min())

    field = crossings[SYMMETRY_LEVEL]
    reach = min(-position.min(), position.max())
    if not np.isnan(field["left"]) and not np.isnan(field["right"]):
        reach = min(reach, -field["left"], field["right"])
    symmetric = position[(position >= 0) & (position <= reach)]
    asymmetry = float(np.max(np.abs(np.interp(symmetric, position, normalized)
                                    - np.interp(-symmetric, position, normalized))))

    widths = {level: crossings[level]["right"] - crossings[level]["left"] for level in levels}
    offsets = {level: (crossings[level]["right"] + crossings[level]["left"]) / 2
               for level in levels}
    return {
        "kind": "profile",
        "curve_type": scan["curve_type"],
        "depth_mm": scan["depth"],
        "unit": scan["unit"],
        "position": position,
        "value": value,
        "normalized": normalized,
        "reference_value": centre,
        "levels": tuple(levels),
        "crossings": crossings,
        "widths": widths,
        "offsets": offsets,
        "flatness_percent": span / 2,
        "flatness_span_percent": span,
        "asymmetry_percentage_points": asymmetry,
        "symmetry_reach_mm": float(reach),
        "truncated": bool(np.isnan(list(widths.values())).any()),
    }


def depthDoseMetrics(scan, levels=DEPTH_LEVELS):
    """
    Derive electron depth-dose metrics from one depth-dose scan.

    The curve is normalised to its maximum. ``dmax`` is the sampled depth of that
    maximum, not an interpolated one, so it is limited by the scan step. Ranges are
    linear interpolations of the falling edge. The bremsstrahlung tail is taken as
    the last sampled point and is reported, not subtracted.

    Args:
        scan (dict): One scan from ``Read.ReadMCC.readMCC``.
        levels (tuple[float, ...]): Relative levels in percent to report ranges at.

    Returns:
        dict: Normalised curve plus range, surface and tail metrics.
    """
    depth, value = _ordered(scan)
    peak = float(value.max())
    if not peak > 0:
        raise ValueError(f"{scan['curve_type']}: maximum dose is not positive")
    normalized = 100 * value / peak
    ranges = {level: fallingCrossing(depth, normalized, level) for level in levels}
    return {
        "kind": "depth_dose",
        "curve_type": scan["curve_type"],
        "depth_mm": scan["depth"],
        "unit": scan["unit"],
        "position": depth,
        "value": value,
        "normalized": normalized,
        "reference_value": peak,
        "levels": tuple(levels),
        "dmax_sampled_mm": float(depth[int(np.argmax(value))]),
        "surface_percent": float(normalized[0]),
        "ranges": ranges,
        "tail_depth_mm": float(depth[-1]),
        "tail_percent": float(normalized[-1]),
        "truncated": bool(np.isnan(list(ranges.values())).any()),
    }


def scanMetrics(scan):
    """
    Derive the metrics appropriate to a scan's curve type.

    Args:
        scan (dict): One scan from ``Read.ReadMCC.readMCC``.

    Returns:
        dict: ``profileMetrics`` or ``depthDoseMetrics`` output.
    """
    curve_type = scan["curve_type"]
    if curve_type in PROFILE_CURVE_TYPES:
        return profileMetrics(scan)
    if curve_type in DEPTH_CURVE_TYPES:
        return depthDoseMetrics(scan)
    raise ValueError(f"Unsupported scan curve type {curve_type!r}")
