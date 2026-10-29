"""Read a dose-space campaign file: measured per-channel doses, their errors and the device map.

A dose campaign is a JSON file next to the measured table it names (for HUS,
GRAS/3D/RadEx/RadEx-HUS/MeasuredDoses/campaign.json). It holds the table, its column
names, the phases, the device type and label of every channel, the first tally of each
detector layer in the GRAS runs, the production spectrum and the systematic fraction.
Relative paths resolve against the campaign file. FitScaledSpectrum,
FitSpectrumModels, JointSpectrumFit, CompareVariant and BuildCoefficientTable read
their run-specific settings from it, so pointing them at another campaign is enough.

Every channel's error is sigma^2 = noise^2 + conversion^2 + (systematic x dose)^2.
"""

import csv
import json
from pathlib import Path

import numpy as np


def loadCampaign(path):
    """Read a campaign file; paths in it become absolute, channel keys integers."""
    path = Path(path).resolve()
    with path.open() as handle:
        campaign = json.load(handle)
    base = path.parent

    def resolve(value):
        return value if value is None else str((base / value).resolve())

    campaign["path"] = str(path)
    campaign["doses"] = resolve(campaign["doses"])
    campaign["spectrum"] = resolve(campaign.get("spectrum"))
    campaign["reference_run"] = resolve(campaign.get("reference_run"))
    for key in ("devices", "labels"):
        campaign[key] = {int(c): v for c, v in campaign.get(key, {}).items()}
    campaign["channels"] = sorted(campaign["devices"])
    if campaign["channels"] != list(range(len(campaign["channels"]))):
        raise ValueError(f"{path}: channels must be 0..n-1")
    joint = campaign.get("joint_fit")
    if joint:
        for key in ("response_matrix", "measured_pdd"):
            joint[key] = resolve(joint[key])
    return campaign


def deviceTypes(campaign):
    """Device types in order of first appearance."""
    return list(dict.fromkeys(campaign["devices"][c] for c in campaign["channels"]))


def layerSlice(campaign, layer):
    """Tally columns of one detector layer."""
    first = campaign["layers"][layer]
    return slice(first, first + len(campaign["channels"]))


def scaleGroups(campaign, channels=None):
    """Calibration-scale groups: one shared scale, and one per device type."""
    channels = np.asarray(campaign["channels"] if channels is None else channels)
    devices = campaign["devices"]
    return {"shared": {"all": channels},
            "per-type": {t: np.array([c for c in channels if devices[c] == t]) for t in deviceTypes(campaign)}}


def phaseDoses(campaign, phase, systematic=None):
    """Measured doses of one phase, ordered by channel.

    Returns:
        dict: measured and sigma [kRad] per channel, the phase fluence, the
            systematic fraction used, and the table rows by channel.
    """
    systematic = campaign["systematic_fraction"] if systematic is None else systematic
    with open(campaign["doses"], newline="") as handle:
        table = {int(r["channel"]): r for r in csv.DictReader(handle) if r["phase"] == phase}
    channels = campaign["channels"]
    missing = [c for c in channels if c not in table]
    if missing:
        raise ValueError(f"{campaign['doses']}: phase {phase!r} lacks channels {missing}")
    for c in channels:
        listed = table[c].get("radfet_type")
        if listed and listed != campaign["devices"][c]:
            raise ValueError(f"channel {c}: table says {listed}, campaign says {campaign['devices'][c]}")

    def column(name):
        return np.array([float(table[c][name] or 0) for c in channels])

    measured = column(campaign["dose_column"])
    sigma = np.sqrt(column(campaign["noise_column"]) ** 2 + column(campaign["conversion_column"]) ** 2
                    + (systematic * measured) ** 2)
    return {"measured": measured, "sigma": sigma, "fluence": float(table[channels[0]][campaign["fluence_column"]]),
            "systematic": systematic, "rows": table}
