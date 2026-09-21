"""Read PTW MEPHYSTO CC-Export measurement files (.mcc)."""

from __future__ import annotations

import numpy as np


def readMCC(file):
    """
    Read every scan block of a PTW MEPHYSTO CC-Export (.mcc) measurement file.

    A single file may hold several scans, for example an inplane and a
    crossplane profile exported together, or a single depth-dose curve.
    Scans are returned in file order; callers that expect one specific curve
    must select it by ``curve_type`` rather than assuming a position.

    Args:
        file (str | Path): Path to the .mcc file.

    Returns:
        list[dict]: One dictionary per scan, each with:
            - 'metadata': dict of the raw ``KEY=VALUE`` header entries as strings
            - 'curve_type': SCAN_CURVETYPE, or None when absent
            - 'depth': SCAN_DEPTH in mm as float, NaN when absent
            - 'unit': MEAS_UNIT, or None when absent
            - 'position': numpy array of the first data column
            - 'value': numpy array of the field detector column
            - 'reference': numpy array of the reference detector column,
              or None when the file has no third column
    """
    scans = []
    metadata = None
    rows = None

    with open(file, "r", encoding="latin-1") as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped:
                continue
            token = stripped.split()[0]

            if token == "BEGIN_SCAN":
                metadata, rows = {}, None
                continue
            if token == "END_SCAN":
                if metadata is None:
                    raise ValueError(f"{file}: END_SCAN without BEGIN_SCAN")
                scans.append(_buildScan(file, metadata, rows))
                metadata, rows = None, None
                continue
            if metadata is None:
                continue

            if stripped == "BEGIN_DATA":
                rows = []
                continue
            if stripped == "END_DATA":
                continue
            if rows is not None:
                rows.append([float(field) for field in stripped.split()])
            elif "=" in stripped:
                key, _, value = stripped.partition("=")
                metadata[key.strip()] = value.strip()

    if metadata is not None:
        raise ValueError(f"{file}: BEGIN_SCAN without END_SCAN")
    if not scans:
        raise ValueError(f"{file}: no scan blocks found")
    return scans


def _buildScan(file, metadata, rows):
    """Validate one parsed scan block and convert it to the documented dict."""
    if not rows:
        raise ValueError(f"{file}: scan {metadata.get('SCAN_CURVETYPE')} has no data")

    widths = {len(row) for row in rows}
    if len(widths) != 1:
        raise ValueError(f"{file}: scan has inconsistent column counts {sorted(widths)}")
    columns = widths.pop()
    if columns not in (2, 3):
        raise ValueError(f"{file}: expected 2 or 3 data columns, found {columns}")

    data = np.asarray(rows, dtype=float)
    if not np.isfinite(data).all():
        raise ValueError(f"{file}: scan contains non-finite data")

    depth = metadata.get("SCAN_DEPTH")
    return {
        "metadata": metadata,
        "curve_type": metadata.get("SCAN_CURVETYPE"),
        "depth": float(depth) if depth is not None else float("nan"),
        "unit": metadata.get("MEAS_UNIT"),
        "position": data[:, 0],
        "value": data[:, 1],
        "reference": data[:, 2] if columns == 3 else None,
    }


def readScan(file, curve_type=None):
    """
    Read exactly one scan from a .mcc file.

    Args:
        file (str | Path): Path to the .mcc file.
        curve_type (str | None): SCAN_CURVETYPE to select, for example 'PDD' or
            'INPLANE_PROFILE'. Required when the file holds more than one scan.

    Returns:
        dict: The selected scan, in the form documented by ``readMCC``.
    """
    scans = readMCC(file)
    if curve_type is not None:
        scans = [scan for scan in scans if scan["curve_type"] == curve_type]
        if not scans:
            raise ValueError(f"{file}: no {curve_type} scan")
    if len(scans) != 1:
        raise ValueError(
            f"{file}: expected a single scan, found {len(scans)}: "
            f"{[scan['curve_type'] for scan in scans]}"
        )
    return scans[0]
