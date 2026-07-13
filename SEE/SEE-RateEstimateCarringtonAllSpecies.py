#!/usr/bin/env python3
"""Estimate SEE rates for the Carrington CREME96 and SAPPHIRE datasets.

Adapted from SEE-RateEstimateCarringtonFinal.py.  By default, both all-species
folders next to this script are analysed.  Alternative dataset directories can
be supplied as positional arguments.
"""

import argparse
import csv
import os
import re
import sys
from pathlib import Path
import numpy as np


# totalLETHistos and its imports live in the shared Python analysis directory.
ANALYSIS_ROOT = Path("/home/anton/triton_work/Python")
sys.path.insert(0, str(ANALYSIS_ROOT))
from Dependencies.TotalLETHistos import totalLETHistos  # noqa: E402


CROSS_SECTIONS = (
    {
        "name": "NanoXplore SEU",
        "L0": 0.11,
        "W": 36.0,
        "S": 4.4,
        "A0": 5.2e-9,
    },
    {
        "name": "Cypress CY62167GE30-45ZXI",
        "L0": 0.1,
        "W": 70.0,
        "S": 1.2,
        "A0": 2.6e-7,
    },
)

# The SAPPHIRE SEF spectra are normalised to a 4015-day (11-year) mission
# fluence (cm-2), so they must be divided by the mission duration to obtain a
# rate.  The CREME96 spectra are already normalised to flux (cm-2 s-1) and
# need no time normalisation.
NORMALISATION_SECONDS = 4015 * 24 * 3600

DATA_ROOT = Path("/home/anton/triton_work/GRAS/LET_Histograms")
DEFAULT_DIRECTORIES = (
    DATA_ROOT / "Carrington-CREME96-All-Species",
    DATA_ROOT / "Carrington-SAPPHIRE-All-Species",
)


def natural_key(path):
    """Sort paths naturally (for example, 2mm before 16mm)."""
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", path.name)]


def cross_section(let_values, parameters):
    """Return the Weibull SEE cross section in cm2/bit."""
    let_values = np.asarray(let_values, dtype=float)
    result = np.zeros_like(let_values)
    mask = let_values > parameters["L0"]
    result[mask] = parameters["A0"] * (
        1.0
        - np.exp(
            -(
                (let_values[mask] - parameters["L0"]) / parameters["W"]
            )
            ** parameters["S"]
        )
    )
    return result


def analyse_directory(directory, parameters):
    """Analyse every species/shielding combination in one dataset."""
    directory = Path(directory).expanduser().resolve()
    if not directory.is_dir():
        print(f"Dataset directory does not exist: {directory}", file=sys.stderr)
        return False

    # Only fluence-normalised datasets (SAPPHIRE SEF) are divided by the
    # mission duration; flux-normalised datasets (CREME96) already are rates.
    normalisation = NORMALISATION_SECONDS if "SAPPHIRE" in directory.name else 1.0

    cross_section_name = parameters["name"]
    output_path = directory / f"SEERates_{cross_section_name}.csv"
    species_folders = sorted(
        (path for path in directory.iterdir() if path.is_dir() and path.name != "NoSEEs"),
        key=natural_key,
    )

    with output_path.open("w", newline="") as output_file:
        writer = csv.writer(output_file)
        writer.writerow(
            (
                "Data",
                "Shielding",
                "Crossection",
                "SEE_Rate",
                "SEE_Error",
                "Relative_SEE_Error",
                "Entries_Contributing_To_SEE",
            )
        )

        for species_folder in species_folders:
            shielding_folders = sorted(
                (path for path in species_folder.iterdir() if path.is_dir()),
                key=natural_key,
            )
            for shielding_folder in shielding_folders:
                result_directory = shielding_folder / "Res"
                let_files = list(result_directory.glob("*LET*.csv")) if result_directory.is_dir() else []
                if not let_files:
                    print(f"No LET histogram found in {result_directory} -> skipping")
                    continue

                let_hist, _ = totalLETHistos(os.fspath(result_directory))
                required_keys = ("lower", "upper", "mean", "value", "error", "entries")
                lengths = [len(let_hist[key]) for key in required_keys]
                if not lengths or any(length != lengths[0] for length in lengths):
                    print(f"Inconsistent histogram lengths in {result_directory} -> skipping")
                    continue

                values = np.asarray(let_hist["value"], dtype=float) / normalisation
                errors = np.asarray(let_hist["error"], dtype=float) / normalisation
                weibull = cross_section(let_hist["mean"], parameters)
                see_values = values * weibull
                see_errors = errors * weibull

                see_rate = float(np.sum(see_values))
                see_rate_error = float(np.sqrt(np.sum(np.square(see_errors))))
                relative_error = see_rate_error / see_rate if see_rate else np.nan
                entries = float(np.sum(np.asarray(let_hist["entries"])[see_values > 0]))

                print(
                    f"{species_folder.name}/{shielding_folder.name}: "
                    f"{see_rate:.6e} +/- {see_rate_error:.2e} s-1 bit-1 "
                    f"({int(entries)} contributing entries)"
                )
                writer.writerow(
                    (
                        species_folder.name,
                        shielding_folder.name,
                        cross_section_name,
                        see_rate,
                        see_rate_error,
                        relative_error,
                        entries,
                    )
                )

    print(f"Wrote {output_path}")
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "directories",
        nargs="*",
        type=Path,
        default=DEFAULT_DIRECTORIES,
        help="dataset directories (defaults to both Carrington all-species folders)",
    )
    args = parser.parse_args()
    success = all(
        analyse_directory(directory, parameters)
        for directory in args.directories
        for parameters in CROSS_SECTIONS
    )
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
