"""Read GRAS non-ionising dose (NID) modules from one result file.

The displacement-damage counterpart of Read.ReadDose, and a sibling of its
readDoseModules. Nothing already in Read/ can do this: every dose reader here
matches GRAS blocks titled TOTAL DOSE, and a NID module writes its own block
titled TOTAL NID with a different unit column. Dependencies.TotalDose reads a
file containing both and returns only the dose tallies.

The shape differs from the dose case and that drives the interface. A dose
module carries several volumes as columns of one block, so readDose returns
arrays indexed by tile. A NID module carries exactly ONE volume, and several
volumes mean several modules, so this returns a dict keyed by module name, the
way readDoseModules does for layered phantoms.

Units are preserved as written. There is no counterpart to the rad-to-kRad
conversion in readDose: displacement damage is quoted in MeV/g and nothing
here rescales it.
"""
import re

# GRAS writes values in plain or exponential notation, sometimes with leading
# whitespace and always comma separated.
_NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"
_DATA_LINE = re.compile(
    rf"^\s*({_NUMBER})\s*,\s*({_NUMBER})\s*,\s*({_NUMBER})\s*,\s*({_NUMBER})\s*$",
    re.MULTILINE,
)
_MODULE_NAME = re.compile(r"'GRAS_MODULE_NAME'\s*,\s*-1\s*,\s*'([^']+)'")
_UNIT_NAME = re.compile(r"'NID'\s*,\s*'([^']+)'")
_TITLE = re.compile(r"'GRAS_DATA_TITLE'\s*,\s*-1\s*,\s*'TOTAL NID'")


def readNID(file, modulePrefix=None):
    """Read every non-ionising dose module in a GRAS CSV output file.

    Args:
        file (str): Path to the GRAS CSV output file.
        modulePrefix (str or None): Optional module-name prefix to retain, for
            files carrying NID modules belonging to different studies.

    Returns:
        dict: Mapping from module name to a dict with ``nid``, ``error``,
        ``entries``, ``non-zeros`` and ``unit``. Empty if the file holds no NID
        module, which is the normal case for a TID-only run.
    """
    with open(file, "r") as stream:
        contents = stream.read()

    modules = {}
    for block in re.split(r"(?=^'\*',)", contents, flags=re.MULTILINE):
        if not _TITLE.search(block):
            continue
        nameMatch = _MODULE_NAME.search(block)
        unitMatch = _UNIT_NAME.search(block)
        values = _DATA_LINE.findall(block)
        if not nameMatch or not unitMatch or not values:
            continue
        name = nameMatch.group(1)
        if modulePrefix is not None and not name.startswith(modulePrefix):
            continue
        nid, error, entries, nonZeros = values[-1]
        modules[name] = {
            "nid": float(nid),
            "error": float(error),
            "entries": int(float(entries)),
            "non-zeros": int(float(nonZeros)),
            "unit": unitMatch.group(1),
        }

    return modules