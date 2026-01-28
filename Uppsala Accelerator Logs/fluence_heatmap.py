"""
Fluence heat-map generator (proton-weighted, cell-center-in-ellipse rule).

- Scans DATA_DIR for *_record_*.csv (mixed headers supported).
- Extracts **scanner chamber** blocks that start with:
    layerId,elementId,xDose,yDose,xPositionSC,yPositionSC,xWidth,yWidth,...
- For each spot (layerId, elementId):
    Qx = ΔxDose, Qy = ΔyDose across all repeats in the file.
    N_protons = (Qx + Qy) / e
    x = median xPositionSC; y = median yPositionSC; wx, wy = median widths.
- Builds a DX × DY (mm) grid. For each spot, adds N_protons to all grid cells
  whose **centers** lie inside the axis-aligned ellipse with semi-axes a=wx/2, b=wy/2.
- Saves a PNG heat map per CSV.

Outputs per file:
  - <csv>.fluence.png   : proton fluence map (protons per bin)
  - <csv>.spots.csv     : per-spot summary (positions, widths, proton counts)

Requirements: pandas, numpy, matplotlib
"""

from pathlib import Path
from io import StringIO
from collections import defaultdict
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ----------------- Configuration -----------------
DATA_DIR        = Path("skandion_irrlogs_20250628_20250628_095105_674.PBS.2.8")  # <-- CHANGE ME
DX              = 0.1   # mm (grid spacing in X)
DY              = 0.1   # mm (grid spacing in Y)
WIDTH_ROUND_MM  = 0.05  # round widths to bucket masks (mm)
MAX_SPOTS       = 800_000  # optional cap; downsample if above (weighted by block stride)
VMIN_PERCENTILE = 1      # imshow contrast floor (percentile of nonzero cells)
VMAX_PERCENTILE = 99.9   # imshow contrast ceiling
# -------------------------------------------------

# Headers
SC_HEADER_PREFIX   = "layerId,elementId,"
SUBMAP_HEADER_PREF = "SUBMAP_NUMBER"

# Scanner-chamber column names
COL_LAY   = "layerId"
COL_ELEM  = "elementId"
COL_X     = "xPositionSC"
COL_Y     = "yPositionSC"
COL_WX    = "xWidth"
COL_WY    = "yWidth"
COL_QX    = "xDose"
COL_QY    = "yDose"

E_CHARGE = 1.602_176_634e-19  # Coulomb

# ----------------- Parsing -----------------

def load_scanner_blocks(path: Path) -> pd.DataFrame:
    """Parse and concatenate all 'layerId,elementId,…' blocks; ignore other blocks."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    dataframes = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith(SC_HEADER_PREFIX):
            header = line
            cols = header.split(",")
            ncols = len(cols)
            j = i + 1
            block = [header]
            while j < len(lines):
                s = lines[j].strip()
                if s.startswith(SC_HEADER_PREFIX) or s.startswith(SUBMAP_HEADER_PREF):
                    break
                # keep only rows that match column count
                if s and not s.startswith("#") and len(s.split(",")) == ncols:
                    block.append(s)
                j += 1
            if len(block) > 1:
                df_block = pd.read_csv(StringIO("\n".join(block)))
                dataframes.append(df_block)
            i = j
        else:
            i += 1

    if not dataframes:
        raise ValueError(f"No scanner-chamber blocks found in {path.name}")

    df = pd.concat(dataframes, ignore_index=True)

    # Coerce numeric types
    for c in (COL_LAY, COL_ELEM, COL_X, COL_Y, COL_WX, COL_WY, COL_QX, COL_QY):
        if c not in df.columns:
            raise ValueError(f"Missing expected column '{c}' in scanner block.")
        df[c] = pd.to_numeric(df[c], errors="coerce")

    return df


# ----------------- Fluence rasterization -----------------

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
    rx = int(np.ceil(max(a, 1e-12) / dx))
    ry = int(np.ceil(max(b, 1e-12) / dy))

    xs = (np.arange(-rx, rx + 1)) * dx
    ys = (np.arange(-ry, ry + 1)) * dy
    X, Y = np.meshgrid(xs, ys, indexing="xy")

    with np.errstate(divide="ignore", invalid="ignore"):
        inside = ((X / max(a, 1e-12)) ** 2 + (Y / max(b, 1e-12)) ** 2) <= 1.0

    inside[ry, rx] = True  # ensure center is included
    return inside.astype(np.uint8)  # 0/1


def build_fluence_heatmap(spots_df: pd.DataFrame, dx=DX, dy=DY):
    """
    Rasterize weighted spots into a fluence matrix using the cell-center-in-ellipse rule.

    Returns: matrix (ny, nx), x0, y0, dx, dy
    """
    x  = spots_df[COL_X].to_numpy()
    y  = spots_df[COL_Y].to_numpy()
    wx = spots_df[COL_WX].to_numpy()
    wy = spots_df[COL_WY].to_numpy()
    w  = spots_df["N_protons"].to_numpy()  # weights

    # semi-axes (full width / 2)
    a = wx / 2.0
    b = wy / 2.0

    # Grid extent to cover all ellipses
    x0, x1, y0, y1 = grid_extent(x, y, a, b, dx=dx, dy=dy)
    nx = int(round((x1 - x0) / dx)) + 1
    ny = int(round((y1 - y0) / dy)) + 1
    M = np.zeros((ny, nx), dtype=np.float64)

    # Optional downsample for speed (take every k-th spot, but scale weights accordingly)
    if len(x) > MAX_SPOTS:
        step = max(1, len(x) // MAX_SPOTS)
        x, y, a, b, w = x[::step], y[::step], a[::step], b[::step], w[::step] * step

    # Bucket by rounded widths to reuse masks
    key_wx = np.round(wx / WIDTH_ROUND_MM) * WIDTH_ROUND_MM
    key_wy = np.round(wy / WIDTH_ROUND_MM) * WIDTH_ROUND_MM

    buckets = defaultdict(list)
    for i, (xi, yi, ai, bi, wi, kwx, kwy) in enumerate(zip(x, y, a, b, w, key_wx, key_wy)):
        buckets[(float(kwx), float(kwy))].append((xi, yi, ai, bi, wi))

    for (_kwx, _kwy), items in buckets.items():
        # Representative mask for this width bucket
        _, _, a0, b0, _ = items[0]
        mask = ellipse_mask(max(a0, 1e-12), max(b0, 1e-12), dx=dx, dy=dy)
        ry, rx = (mask.shape[0] // 2), (mask.shape[1] // 2)

        for xi, yi, ai, bi, wi in items:
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

            # Weighted stamp: add wi to cells where mask==1
            M[ys0:ys1, xs0:xs1] += wi * mask[my0:my1, mx0:mx1]

    return M, x0, y0, dx, dy


# ----------------- Spots & plotting -----------------

def compute_spot_table(df_sc: pd.DataFrame) -> pd.DataFrame:
    """Aggregate rows into (layerId, elementId) spots with protons and geometry."""
    df = df_sc[[COL_LAY, COL_ELEM, COL_X, COL_Y, COL_WX, COL_WY, COL_QX, COL_QY]].copy()
    df.sort_values([COL_LAY, COL_ELEM], inplace=True)

    def end_minus_start(s):
        s = pd.to_numeric(s, errors="coerce").dropna()
        if s.empty:
            return np.nan
        return s.iloc[-1] - s.iloc[0]

    aggs = (
        df.groupby([COL_LAY, COL_ELEM], sort=False)
          .apply(lambda g: pd.Series({
              "x_mm": g[COL_X].median(),
              "y_mm": g[COL_Y].median(),
              "wx_mm": g[COL_WX].median(),
              "wy_mm": g[COL_WY].median(),
              "Qx_C": end_minus_start(g[COL_QX]),
              "Qy_C": end_minus_start(g[COL_QY]),
              "n_samples": len(g),
          }))
          .reset_index()
    )
    aggs["Q_C"] = aggs["Qx_C"].fillna(0.0) + aggs["Qy_C"].fillna(0.0)
    aggs["N_protons"] = aggs["Q_C"] / E_CHARGE
    return aggs


def save_heatmap_png(matrix, x0, y0, dx, dy, out_png: Path, title: str = ""):
    """
    Save matrix as a heat map PNG. Axes are labeled in mm; colorbar = protons/bin.
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
        vmin = max(vmin, np.finfo(float).eps)
        vmax = max(vmax, vmin * 1.01)
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
    cbar.set_label("Protons per bin")
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    if not title:
        title = f"Proton fluence heat map (bins: {DX} × {DY} mm²)"
    ax.set_title(title)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, format="png")
    plt.close(fig)


def process_file(csv_path: Path) -> None:
    df_sc = load_scanner_blocks(csv_path)
    spots = compute_spot_table(df_sc)

    if spots.empty:
        raise ValueError("No spots computed from scanner blocks.")

    # Save spot table next to source CSV
    out_spots = csv_path.with_suffix(".spots.csv")
    spots.to_csv(out_spots, index=False)

    # Build weighted fluence matrix
    M, x0, y0, dx, dy = build_fluence_heatmap(spots.rename(columns={"x_mm": COL_X, "y_mm": COL_Y, "wx_mm": COL_WX, "wy_mm": COL_WY}), dx=DX, dy=DY)
    out_png = csv_path.with_suffix(".fluence.png")
    ttl = f"{csv_path.name} – Proton fluence (ΣN={spots['N_protons'].sum():.3e})"
    save_heatmap_png(M, x0, y0, dx, dy, out_png, title=ttl)

    print(f"Saved: {out_spots.name} (spots={len(spots)})")
    print(f"Saved: {out_png.name}   (matrix={M.shape[0]}×{M.shape[1]})")


def main() -> int:
    if not DATA_DIR.is_dir():
        print(f"DATA_DIR does not exist or is not a directory: {DATA_DIR}", file=sys.stderr)
        return 1

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
