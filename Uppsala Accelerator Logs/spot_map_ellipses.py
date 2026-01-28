#!/usr/bin/env python3
"""
Batch heat-map generator (cell-center-in-ellipse rule).

- Scans DATA_DIR for *.csv.
- Reads the SUBMAP table starting at 'SUBMAP_NUMBER'.
- Builds a 0.1 mm grid; for each spot, increments all grid cells whose
  **centers** lie inside the axis-aligned ellipse with semi-axes a=wx/2, b=wy/2.
- Accumulates bin counts into a matrix and saves a heat map PDF (same base
  filename as the CSV).

Requirements: pandas, numpy, matplotlib
"""

from pathlib import Path
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ----------------- Configuration -----------------
DATA_DIR        = Path("skandion_irrlogs_20250628_20250628_095105_674.PBS.2.8/Combined")   # <-- CHANGE ME
DX              = 0.01  # mm (grid spacing in X)
DY              = 0.01  # mm (grid spacing in Y)
WIDTH_ROUND_MM  = 0.01  # round widths to this step for bucketing (mm)
MAX_SPOTS       = 400_000  # optional cap; downsample if above
VMIN_PERCENTILE = 0      # imshow contrast floor (percentile of nonzero cells)
VMAX_PERCENTILE = 100     # imshow contrast ceiling
# -------------------------------------------------

# ---- Column names for the SC (scanner chamber) log format ----
HEADER_PREFIX = "layerId,elementId,"
COL_X  = "xPositionSC"
COL_Y  = "yPositionSC"
COL_WX = "xWidth"
COL_WY = "yWidth"



def find_submap_header_line(path: Path) -> int:
    """Return 1-based line number where the data header starts."""
    with path.open("r", errors="replace") as f:
        for i, line in enumerate(f, start=1):
            if line.startswith(HEADER_PREFIX):
                return i
    raise ValueError(f"Could not find data header starting with '{HEADER_PREFIX}' in {path.name}.")



def load_submap_df(path: Path) -> pd.DataFrame:
    start = find_submap_header_line(path)
    return pd.read_csv(path, skiprows=start - 1, low_memory=False)


def filter_valid(df: pd.DataFrame) -> pd.DataFrame:
    req = [COL_X, COL_Y, COL_WX, COL_WY]
    missing = [c for c in req if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    for c in req:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    valid = df[
        df[COL_X].between(-10_000, 10_000)
        & df[COL_Y].between(-10_000, 10_000)
        & df[COL_WX].between(0, 1_000)
        & df[COL_WY].between(0, 1_000)
    ].dropna(subset=req)

    return valid



def grid_extent(x, y, a, b, dx=DX, dy=DY, margin=0.0):
    """
    Compute grid extents (inclusive) from min/max of (center ± semi-axis) with margin.
    Returns (x0, x1, y0, y1) aligned to grid.
    """
    xmin = np.min(x - a) - margin
    xmax = np.max(x + a) + margin
    ymin = np.min(y - b) - margin
    ymax = np.max(y + b) + margin

    x0 = np.floor(xmin / dx) * dx
    x1 = np.ceil(xmax / dx) * dx
    y0 = np.floor(ymin / dy) * dy
    y1 = np.ceil(ymax / dy) * dy
    return x0, x1, y0, y1


def ellipse_mask(a, b, dx=DX, dy=DY):
    """
    Boolean mask for an axis-aligned ellipse using **cell centers**.
    Semi-axes: a, b; spacings: dx, dy. The mask is centered on the middle cell.
    A cell is 1 iff its center satisfies (x/a)^2 + (y/b)^2 <= 1.
    """
    # number of cells to cover semi-axes (rounded up)
    rx = int(np.ceil(a / dx))
    ry = int(np.ceil(b / dy))
    if rx < 0: rx = 0
    if ry < 0: ry = 0

    # grid of centers (relative to ellipse center)
    xs = (np.arange(-rx, rx + 1)) * dx
    ys = (np.arange(-ry, ry + 1)) * dy
    X, Y = np.meshgrid(xs, ys, indexing="xy")

    with np.errstate(divide="ignore", invalid="ignore"):
        inside = ((X / a) ** 2 + (Y / b) ** 2) <= 1.0

    # ensure center included if a or b degenerate
    inside[ry, rx] = True
    return inside.astype(np.uint8)  # 0/1


def build_heatmap(spots_df: pd.DataFrame, dx=DX, dy=DY):
    """
    Rasterize spots into a heat map matrix using the **cell-center-in-ellipse** rule.

    Returns: matrix (ny, nx), x0, y0, dx, dy
    """
    x  = spots_df[COL_X].to_numpy()
    y  = spots_df[COL_Y].to_numpy()
    wx = spots_df[COL_WX].to_numpy()
    wy = spots_df[COL_WY].to_numpy()

    # semi-axes (full width / 2)
    a = wx / 2.0
    b = wy / 2.0

    # Grid extent to cover all ellipses
    x0, x1, y0, y1 = grid_extent(x, y, a, b, dx=dx, dy=dy)
    nx = int(round((x1 - x0) / dx)) + 1
    ny = int(round((y1 - y0) / dy)) + 1
    M = np.zeros((ny, nx), dtype=np.uint32)

    # Optional downsample for speed
    if len(x) > MAX_SPOTS:
        step = max(1, len(x) // MAX_SPOTS)
        x, y, a, b = x[::step], y[::step], a[::step], b[::step]

    # Bucket by rounded full widths, so we can reuse masks
    key_wx = np.round(wx / WIDTH_ROUND_MM) * WIDTH_ROUND_MM
    key_wy = np.round(wy / WIDTH_ROUND_MM) * WIDTH_ROUND_MM
    buckets = defaultdict(list)
    for i, (xi, yi, ai, bi, kwx, kwy) in enumerate(zip(x, y, a, b, key_wx, key_wy)):
        buckets[(float(kwx), float(kwy))].append((xi, yi, ai, bi))

    for (_kwx, _kwy), items in buckets.items():
        # Build a mask once for a representative (a, b)
        _, _, a0, b0 = items[0]
        a0 = max(a0, 1e-9)
        b0 = max(b0, 1e-9)
        mask = ellipse_mask(a0, b0, dx=dx, dy=dy)
        ry, rx = (mask.shape[0] // 2), (mask.shape[1] // 2)

        for xi, yi, ai, bi in items:
            # Center index of this spot on the global grid
            cx = int(round((xi - x0) / dx))
            cy = int(round((yi - y0) / dy))

            # Target slice bounds
            x_start = cx - rx
            x_end   = cx + rx + 1
            y_start = cy - ry
            y_end   = cy + ry + 1

            # Clip to matrix bounds
            xs0 = max(0, x_start)
            ys0 = max(0, y_start)
            xs1 = min(nx, x_end)
            ys1 = min(ny, y_end)
            if xs0 >= xs1 or ys0 >= ys1:
                continue

            # Corresponding slice in mask
            mx0 = xs0 - x_start
            my0 = ys0 - y_start
            mx1 = mx0 + (xs1 - xs0)
            my1 = my0 + (ys1 - ys0)

            # Stamp (add 1s) into matrix
            M[ys0:ys1, xs0:xs1] += mask[my0:my1, mx0:mx1]

    return M, x0, y0, dx, dy


def save_heatmap_pdf(matrix, x0, y0, dx, dy, out_pdf: Path):
    """
    Save matrix as a heat map PDF. Axes are labeled in mm.
    """
    ny, nx = matrix.shape
    x1 = x0 + dx * (nx - 1)
    y1 = y0 + dy * (ny - 1)

    fig, ax = plt.subplots(figsize=(10, 7))
    data = matrix.astype(float)
    nonzero = data[data > 0]
    if nonzero.size > 0:
        vmin = np.percentile(nonzero, VMIN_PERCENTILE)
        vmax = np.percentile(nonzero, VMAX_PERCENTILE)
        vmin = max(vmin, 1.0)
        vmax = max(vmax, vmin + 1.0)
        im = ax.imshow(
            data,
            origin="lower",
            extent=[x0, x1, y0, y1],
            aspect="equal",
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
        )
    else:
        im = ax.imshow(
            data,
            origin="lower",
            extent=[x0, x1, y0, y1],
            aspect="equal",
            interpolation="nearest",
        )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Bin count (cell centers inside ellipse)")
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_title(f"Beam spot heat map ({DX} mm bins)")

    fig.tight_layout()
    fig.savefig(out_pdf, dpi=300)
    plt.close(fig)


def save_heatmap_png(matrix, x0, y0, dx, dy, out_png: Path):
    """
    Save matrix as a heat map PNG. Axes are labeled in mm.
    """
    ny, nx = matrix.shape
    x1 = x0 + dx * (nx - 1)
    y1 = y0 + dy * (ny - 1)

    fig, ax = plt.subplots(figsize=(10, 7))
    data = matrix.astype(float)
    nonzero = data[data > 0]
    if nonzero.size > 0:
        vmin = np.percentile(nonzero, VMIN_PERCENTILE)
        vmax = np.percentile(nonzero, VMAX_PERCENTILE)
        vmin = max(vmin, 1.0)
        vmax = max(vmax, vmin + 1.0)
        im = ax.imshow(
            data,
            origin="lower",
            extent=[x0, x1, y0, y1],
            aspect="equal",
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
        )
    else:
        im = ax.imshow(
            data,
            origin="lower",
            extent=[x0, x1, y0, y1],
            aspect="equal",
            interpolation="nearest",
        )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Bin count (cell centers inside ellipse)")
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_title(f"Beam spot heat map ({DX} mm bins)")

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, format="png")
    plt.close(fig)


def process_file(csv_path: Path):
    df = load_submap_df(csv_path)
    valid = filter_valid(df)
    if valid.empty:
        raise ValueError("No valid spot rows after filtering.")
    M, x0, y0, dx, dy = build_heatmap(valid, dx=DX, dy=DY)
    out_png = csv_path.with_suffix(".png")
    save_heatmap_png(M, x0, y0, dx, dy, out_png)
    print(f"Saved heat map: {out_png.name}  (shape={M.shape[0]}×{M.shape[1]})")


def main() -> int:
    if not DATA_DIR.is_dir():
        print(f"DATA_DIR does not exist or is not a directory: {DATA_DIR}", file=sys.stderr)
        return 1

    # csvs = sorted(DATA_DIR.glob("*.csv"))
    csvs = sorted(DATA_DIR.glob("*_record_*.csv"))
    if not csvs:
        print(f"No CSV files found in {DATA_DIR}")
        return 0

    for csv_path in csvs:
        try:
            print(f"Processing {csv_path.name} …")
            process_file(csv_path)
        except Exception as e:
            print(f"  !! Skipped {csv_path.name}: {e}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
