"""Read a TotalDose CSV into the dictionary format of Dependencies.TotalDose.totalDose."""

import csv

import numpy as np

HEADER = ["Tile", "Dose [kRad]", "Error [kRad]", "Non-Zero Entries"]


def readTotalDoseCSV(file, scale=1):
    """
    Reads a TotalDose_<name>.csv file, such as the Heavy-Ion totals of Plotting/AllSpeciesShieldingCurves.py.

    The values are as stored: each dose is rounded to two significant digits of its error.

    Args:
        file (str): Path to the CSV file.
        scale (float): Factor for the dose and the error, for example the mission duration in seconds
            to merge a total that is per second with mission doses.

    Returns:
        dict: Dictionary with numpy arrays, one value per tile, usable with mergeTotalDose:
            - 'dose': Dose in kRad
            - 'error': Absolute error in kRad
            - 'entries': Zeros, because the CSV does not store the number of entries
            - 'non-zeros': Number of non-zero entries
    """
    with open(file, newline="") as f:
        rows = [row for row in csv.reader(f, skipinitialspace=True) if row]

    if not rows or rows[0] != HEADER:
        raise ValueError(f"{file}: not a TotalDose CSV")
    rows = rows[1:]
    if not rows or [int(row[0]) for row in rows] != list(range(len(rows))):
        raise ValueError(f"{file}: tiles are not 0..n-1 in order")

    return {
        'dose': np.array([float(row[1]) for row in rows]) * scale,
        'error': np.array([float(row[2]) for row in rows]) * scale,
        'entries': np.zeros(len(rows), dtype=int),
        'non-zeros': np.array([int(float(row[3])) for row in rows]),
    }
