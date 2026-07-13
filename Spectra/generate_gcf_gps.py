#!/usr/bin/env python3
"""Convert a SPENVIS GCF (ISO 15390 cosmic-ray) ion table to Geant4 GPS macros.

The ion table is read with readSpenvis_gcf, which returns total kinetic energy
in MeV and flux in cm-2 s-1 (differential per MeV) as required by Geant4 GPS.
One macro is written for every element having a non-zero spectrum, using the
shared writer from generate_sefflare_gps.

By default this converts the Carrington GEO GCR data set used by this project
and writes the macros into a ``GCR-AllSpecies`` directory beside the input.
Pass positional paths to use another GCF file or output directory.
"""

from __future__ import annotations

import argparse
from pathlib import Path

# Importing generate_sefflare_gps also puts the package root on sys.path,
# which the Read import below relies on.
from generate_sefflare_gps import MASS_NUMBERS, write_macros
from Read.ReadSpenvis_gcf import readSpenvis_gcf


DEFAULT_INPUT = Path("/home/anton/triton_work/Spectra/Carrington/GEO/spenvis_gcf.txt")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path, default=DEFAULT_INPUT)
    parser.add_argument(
        "output",
        nargs="?",
        type=Path,
        help="output directory (default: GCR-AllSpecies beside input)",
    )
    parser.add_argument("--prefix", default="GEO-GCR")
    args = parser.parse_args()

    output = args.output or args.input.parent / "GCR-AllSpecies"
    data = readSpenvis_gcf(args.input, masses=MASS_NUMBERS)
    manifest, skipped = write_macros(data, output, args.prefix, args.input.name)
    print(f"Wrote {len(manifest)} GPS macros and manifest to {output}")
    if skipped:
        print("Skipped zero-flux species: " + ", ".join(skipped))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
