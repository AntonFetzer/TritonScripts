#!/usr/bin/env python3
"""Convert a SPENVIS SEFFLARE solar-ion table to Geant4 GPS macros.

The ion table is read with readSpenvis_sefflare_Ions, which returns total
kinetic energy in MeV and flux in cm-2 s-1 (differential per MeV) as required
by Geant4 GPS.  One macro is written for every element having a non-zero
spectrum.  The generated layout and normalization aliases match SPENVIS's
manual GPS export.

By default this converts the Carrington CREME96 data set used by this project
and writes the macros into a ``GPS-AllSpecies`` directory beside the input.
Pass positional paths to use another SEFFLARE file or output directory.

This module also holds the element tables and macro writer shared by
generate_sef_gps.py and generate_gcf_gps.py.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Read.ReadSpenvis_sefflare import readSpenvis_sefflare_Ions


DEFAULT_INPUT = Path(
    "/home/anton/triton_work/Spectra/Carrington/GEO-Extreme/CREME96/"
    "spenvis_sefflare.txt"
)

SYMBOLS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe "
    "Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In "
    "Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf "
    "Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U"
).split()

ELEMENT_NAMES = (
    "Proton Helium Lithium Beryllium Boron Carbon Nitrogen Oxygen Fluorine "
    "Neon Sodium Magnesium Aluminium Silicon Phosphorus Sulfur Chlorine Argon "
    "Potassium Calcium Scandium Titanium Vanadium Chromium Manganese Iron "
    "Cobalt Nickel Copper Zinc Gallium Germanium Arsenic Selenium Bromine "
    "Krypton Rubidium Strontium Yttrium Zirconium Niobium Molybdenum Technetium "
    "Ruthenium Rhodium Palladium Silver Cadmium Indium Tin Antimony Tellurium "
    "Iodine Xenon Caesium Barium Lanthanum Cerium Praseodymium Neodymium "
    "Promethium Samarium Europium Gadolinium Terbium Dysprosium Holmium Erbium "
    "Thulium Ytterbium Lutetium Hafnium Tantalum Tungsten Rhenium Osmium Iridium "
    "Platinum Gold Mercury Thallium Lead Bismuth Polonium Astatine Radon Francium "
    "Radium Actinium Thorium Protactinium Uranium"
).split()

# Representative isotope for each elemental spectrum.  The source data contain
# elements, not isotope-resolved spectra.  Stable elements use their most
# abundant natural isotope; radioactive elements use the conventional isotope
# represented by bracketed standard atomic weights.  Edit this table if a
# different isotope composition is required by a particular simulation.
MASS_NUMBERS = (
    1, 4, 7, 9, 11, 12, 14, 16, 19, 20, 23, 24, 27, 28, 31, 32, 35, 40,
    39, 40, 45, 48, 51, 52, 55, 56, 59, 58, 63, 64, 69, 74, 75, 80, 79,
    84, 85, 88, 89, 90, 93, 98, 98, 102, 103, 106, 107, 114, 115, 120,
    121, 130, 127, 132, 133, 138, 139, 140, 141, 142, 145, 152, 153, 158,
    159, 164, 165, 166, 169, 174, 175, 180, 181, 184, 187, 192, 193, 195,
    197, 202, 205, 208, 209, 209, 210, 222, 223, 226, 227, 232, 231, 238,
)


def render_macro(data, species_index, mass_number, source_name):
    """Render one SPENVIS-compatible GPS macro, or None for a zero spectrum.

    ``data`` is the reader's array: Data[Z-1, point, (Energy MeV, Integral,
    Differential per MeV)], already in cm-2 (s-1) units.
    """
    atomic_number = species_index + 1
    spectrum = data[species_index]
    points = [(energy, flux) for energy, _, flux in spectrum if flux > 0]
    if not points:
        return None

    # SPENVIS uses the tabulated integral flux at the lowest energy as the
    # absolute spectrum normalization.  GPS normalizes the arb histogram itself.
    normalization = spectrum[0, 1]
    if normalization <= 0:
        raise ValueError(
            f"Z={atomic_number}: positive differential spectrum but "
            "non-positive integral normalization"
        )

    lines = ["#", f"# Generated from {source_name}", "#", "#Source definition", "#"]
    if atomic_number == 1 and mass_number == 1:
        lines.append("/gps/particle  proton")
    else:
        lines.extend(
            (
                "/gps/particle  ion",
                f"/gps/ion        {atomic_number:d}        {mass_number:d}        "
                f"{atomic_number:d}   0.000000E+00",
            )
        )
    lines.extend(("/gps/ene/type  Arb", "/gps/hist/type arb"))
    lines.extend(f"/gps/hist/point   {energy:.6E}   {flux:.6E}" for energy, flux in points)
    lines.extend(
        (
            "/gps/hist/inter Log",
            "/gps/ang/type cos",
            "/gps/ang/mintheta    0.000E+00 deg",
            "/gps/ang/maxtheta    9.000E+01 deg",
            "/gps/source/list ",
            "",
            "#",
            "#Normalisation",
            "#",
            f'/control/alias  NORM_FACTOR_SPECTRUM "   {normalization:.6E} "',
            '/control/alias  NORM_FACTOR_ANGULAR "   2.500000E-01 "',
            "",
        )
    )
    return "\n".join(lines)


def write_macros(data, output_dir: Path, prefix: str, source_name: str):
    """Write all non-zero species macros and a CSV manifest."""
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = []
    skipped = []
    for index, (symbol, name, mass_number) in enumerate(
        zip(SYMBOLS, ELEMENT_NAMES, MASS_NUMBERS)
    ):
        macro = render_macro(data, index, mass_number, source_name)
        if macro is None:
            skipped.append(symbol)
            continue
        filename = f"{prefix}-{name}.mac"
        (output_dir / filename).write_text(macro, encoding="ascii")
        manifest.append((index + 1, symbol, mass_number, data[index, 0, 1], filename))

    manifest_path = output_dir / f"{prefix}-manifest.csv"
    with manifest_path.open("w", encoding="ascii", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("Z", "symbol", "A", "integral_flux_cm-2_s-1", "macro"))
        writer.writerows(manifest)
    return manifest, skipped


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path, default=DEFAULT_INPUT)
    parser.add_argument(
        "output",
        nargs="?",
        type=Path,
        help="output directory (default: GPS-AllSpecies beside input)",
    )
    parser.add_argument("--prefix", default="GEO-CREME96-Solar")
    args = parser.parse_args()

    output = args.output or args.input.parent / "GPS-AllSpecies"
    data = readSpenvis_sefflare_Ions(args.input, masses=MASS_NUMBERS)
    manifest, skipped = write_macros(data, output, args.prefix, args.input.name)
    print(f"Wrote {len(manifest)} GPS macros and manifest to {output}")
    if skipped:
        print("Skipped zero-flux species: " + ", ".join(skipped))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
