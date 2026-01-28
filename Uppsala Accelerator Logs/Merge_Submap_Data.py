#!/usr/bin/env python3
"""
Collect all valid SUBMAP rows from CSVs matching "*_record_*.csv" in DATA_DIR
and write them to a single CSV. Adds a 'source_file' column.
"""
from pathlib import Path
import sys

import pandas as pd
import numpy as np

# Configuration (adjust if needed)
DATA_DIR = Path("skandion_irrlogs_20250628_20250628_095105_674.PBS.2.8")


def find_submap_header_line(path: Path) -> int:
    """Return 1-based line number where 'SUBMAP_NUMBER' starts."""
    with path.open("r", errors="replace") as f:
        for i, line in enumerate(f, start=1):
            if "SUBMAP_NUMBER" in line:
                return i
    raise ValueError(f"Could not find 'SUBMAP_NUMBER' header in {path.name}.")


def load_submap_df(path: Path) -> pd.DataFrame:
    """Load CSV starting from the SUBMAP table header."""
    start = find_submap_header_line(path)
    return pd.read_csv(path, skiprows=start - 1, low_memory=False)


def filter_valid(df: pd.DataFrame) -> pd.DataFrame:
    """Return rows with required numeric columns and sane ranges."""
    req = ["X_POSITION(mm)", "Y_POSITION(mm)", "X_WIDTH(mm)", "Y_WIDTH(mm)"]
    missing = [c for c in req if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    for c in req:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    valid = df[
        df["X_POSITION(mm)"].between(-10_000, 10_000)
        & df["Y_POSITION(mm)"].between(-10_000, 10_000)
        & df["X_WIDTH(mm)"].between(0, 1_000)
        & df["Y_WIDTH(mm)"].between(0, 1_000)
    ].dropna(subset=req)

    return valid


def collect_valid_rows(out_csv: Path):
    """Collect valid rows from all *_record_*.csv files and write to out_csv."""
    csvs = sorted(DATA_DIR.glob("*_record_*.csv"))
    if not csvs:
        print(f"No CSV files found in {DATA_DIR}", file=sys.stderr)
        return 1

    rows = []
    for csv_path in csvs:
        print(f"Processing {csv_path.name}...")
        try:
            df = load_submap_df(csv_path)
            valid = filter_valid(df)
            if not valid.empty:
                v = valid.copy()
                v["source_file"] = csv_path.name
                rows.append(v)
        except Exception as e:
            print(f"  !! Skipped {csv_path.name}: {e}", file=sys.stderr)

    if not rows:
        print("No valid rows collected.", file=sys.stderr)
        return 1

    combined = pd.concat(rows, ignore_index=True)
    combined.to_csv(out_csv, index=False)
    print(f"Wrote {len(combined)} rows to {out_csv}")
    return 0


def main(argv=None) -> int:
    argv = argv or sys.argv
    out = Path(argv[1]) if len(argv) >= 2 else DATA_DIR / "combined_valid_spots.csv"
    return collect_valid_rows(out)


if __name__ == "__main__":
    raise SystemExit(main())