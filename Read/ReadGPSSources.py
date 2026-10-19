"""Read the mono-energetic components of a multi-source GPS macro."""

import re

import numpy as np

UNIT_TO_MEV = {"eV": 1e-6, "keV": 1e-3, "MeV": 1.0, "GeV": 1e3}


def readGPSSources(path):
    """Energies and relative intensities of the sources in a GPS macro.

    Each source starts with ``/gps/source/add <intensity>`` and needs a mono
    energy, ``/gps/ene/mono`` or ``/gps/energy``, before the next source starts.
    Comments and other commands are ignored. Raises ValueError for a source
    without an energy, a non-mono energy type or an unknown unit.

    Returns:
        dict: 'energy' (MeV), 'intensity' (as written), 'fraction' (intensity
        divided by its sum, the probability GPS samples each source with) and
        'particle' (one name per source, None where the macro sets none).
    """
    sources = []
    with open(path, "r", encoding="utf-8") as handle:
        for number, raw in enumerate(handle, start=1):
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            words = line.split()
            command = words[0]
            if command == "/gps/source/add":
                sources.append({"intensity": float(words[1]), "energy": None,
                                "particle": None, "line": number})
            elif not sources:
                continue
            elif command == "/gps/particle":
                sources[-1]["particle"] = words[1]
            elif command == "/gps/ene/type" and words[1] != "Mono":
                raise ValueError(f"{path}:{number}: energy type {words[1]}, expected Mono")
            elif command in ("/gps/ene/mono", "/gps/energy"):
                unit = words[2] if len(words) > 2 else "MeV"
                if unit not in UNIT_TO_MEV or not re.fullmatch(r"[-+0-9.eE]+", words[1]):
                    raise ValueError(f"{path}:{number}: cannot read energy '{line}'")
                sources[-1]["energy"] = float(words[1]) * UNIT_TO_MEV[unit]

    if not sources:
        raise ValueError(f"{path}: no /gps/source/add")
    missing = [s["line"] for s in sources if s["energy"] is None]
    if missing:
        raise ValueError(f"{path}: sources starting at lines {missing} set no energy")
    intensity = np.array([s["intensity"] for s in sources], dtype=float)
    if np.any(intensity < 0) or intensity.sum() <= 0:
        raise ValueError(f"{path}: intensities must be non-negative with a positive sum")
    return {
        "energy": np.array([s["energy"] for s in sources], dtype=float),
        "intensity": intensity,
        "fraction": intensity / intensity.sum(),
        "particle": [s["particle"] for s in sources],
    }
