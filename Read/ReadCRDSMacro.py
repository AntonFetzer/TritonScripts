"""Read the source-spectrum alias from a CRDS GRAS macro."""

from pathlib import Path
import re


def spectrumName(runPath):
    """Name of the source spectrum a run used, from its A.mac Spectrum alias."""
    macro = Path(runPath) / "A.mac"
    match = re.search(r"^\s*/control/alias\s+Spectrum\s+(\S+)",
                      macro.read_text(encoding="utf-8"), re.MULTILINE)
    if match is None:
        raise RuntimeError(f"no Spectrum alias in {macro}")
    return match.group(1).rstrip("/").split("/")[-1]
