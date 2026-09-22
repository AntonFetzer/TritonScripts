"""Campaign definitions shared by the CRDS TID, DDD and series drivers.

Every CRDS campaign runs the same instrument with the same detector macro, so
the tally lists are common and only the campaign directory, the beam and the
run list differ. The drivers take ``--campaign`` and read everything
campaign-specific from here; nothing campaign-specific belongs in a driver.

    Plotting/CRDSTID.py     TID in the RadFET, one run
    Plotting/CRDSDDD.py     DDD in the LED, one run, the same result files
    Plotting/CRDSSeries.py  every run of a campaign, both instruments

Uppsala additionally has a ``delivery`` block: its two exposures were not
uniform over the instrument, so the series driver combines the coefficients
with the local fluence at each part from the accelerator-log maps, through
Dependencies/CRDSDelivery.py.

There is no source area anywhere in this analysis. GRAS normalises in
FLUENCE/CURRENT mode and applies fluence * source_surface / n_primaries before
writing, so every coefficient is already per incident particle/cm2.
"""

from pathlib import Path

from Read.ReadCRDSMacro import spectrumName

CRDS_PATH = Path("/scratch/work/fetzera1/GRAS/CRDS")
UPPSALA_LOGS = Path("/scratch/work/fetzera1/Radiation Testing/2025-06 Proton/"
                    "Uppsala Accelerator Logs/results/final/aggregates")

# Tally order is fixed by CRDS/CRDSRadFETLEDDetector.mac, which every campaign
# listed here uses. Change both together.
TILES = [
    (0, "VT01_gox_0_PV", "gate oxide"),
    (1, "VT01_die_0_PV", "silicon die"),
]
MODULES = [
    ("nidActive", "LED_active_0_PV", "active layer"),
    ("nidSubstrate", "LED_die_0_PV", "GaAs substrate"),
]

# Material of each tallied volume, from the VT01 and LED modules. GRAS scores
# dose and DDD per gram of the volume's own material, so this is the material
# the units refer to.
MATERIALS = {
    "VT01_gox_0_PV": "SiO2",
    "VT01_die_0_PV": "Si",
    "LED_active_0_PV": "GaAs",
    "LED_die_0_PV": "GaAs",
}

# Each run: folder, Slurm job and the number of result files it produced. The
# file count is the completeness guard of Dependencies.AggregateRun; set it to
# what the run actually produced, which an autoSeed collision can make fewer
# than the array size.
#
# ``steps`` pairs runs as (numerator, denominator); the series driver reports
# each tally's change between them. ``sweep``, where the runs vary one numeric
# parameter, names that run key and the plot the series driver draws against it.
# ``notCoefficients`` names (run, quantity) pairs whose tally is booked but is
# not a usable coefficient; they are flagged and left out of the steps.
# ``delivery``, where present, turns coefficients into delivered values; see
# Dependencies/CRDSDelivery.py.
CAMPAIGNS = {
    "HUS": {
        "path": CRDS_PATH / "CRDS-HUS",
        "description": "HUS TrueBeam PDD electron field",
        "particle": "electron",
        "defaultFolder": "PDD-Electron-CRDS-20mmField",
        # Campaign reference normalisation, the same assumed fluence the
        # RadEx-HUS results are quoted at. Not an MU calibration.
        "referenceFluence": 2e12,
        "referenceNote": "assumed",
        "runs": {
            "20 x 20 mm": {"folder": "PDD-Electron-CRDS-20mmField",
                           "job": 20392947, "files": 100,
                           "field": "20 x 20"},
            "100 x 100 mm": {"folder": "PDD-Electron-CRDS-100mmField",
                             "job": 20392948, "files": 100,
                             "field": "100 x 100"},
            "200 x 200 mm": {"folder": "PDD-Electron-CRDS-200mmField",
                             "job": 20392949, "files": 100,
                             "field": "200 x 200"},
            "20 x 20 mm, no table": {"folder": "PDD-Electron-CRDS-noTable",
                                     "job": 20392950, "files": 100,
                                     "field": "20 x 20"},
        },
        # Successive widenings, then the table control, each against its
        # reference.
        "steps": [("20 x 20 mm", "100 x 100 mm"),
                  ("100 x 100 mm", "200 x 200 mm"),
                  ("20 x 20 mm, no table", "20 x 20 mm")],
        "seriesOutput": "FieldSizeControl_CRDS-HUS.csv",
        # The runs differ in field size and in the table, not along one axis.
        "sweep": None,
        # RadEx-HUS Ch9, the exposed 0 mm channel, is the nearest thing to a
        # bare device in the same field with the same detector model: job
        # 20315277 gave 61.435 kRad on its gate oxide at the same reference
        # fluence. CRDS sits at the same 200 mm air distance but has no
        # aluminium top plate, no spacer and a 1.6 mm board instead of
        # 6.08 mm, so this is context, not a prediction.
        "oxideComparison": {"label": "RadEx-HUS Ch9 oxide",
                            "kRad": 61.435,
                            "note": "different shielding and field size, "
                                    "context only"},
        "activeSubstrateNote": (
            "These should agree: DDD is intensive and a 6 MeV electron field\n"
            "  is uniform across 0.18 mm of GaAs. A divergence means it no "
            "longer is."),
    },
    "Kumpula": {
        "path": CRDS_PATH / "CRDS-Kumpula",
        "description": "Kumpula 10 MeV proton beam, filed-LED epoxy sweep",
        "particle": "proton",
        "defaultFolder": "10MeVProton-CRDS-500umEpoxy",
        # The delivered fluence is not yet in the model; it needs the
        # accelerator log, as Uppsala's did. Results stay per proton/cm2
        # unless a fluence is passed with --reference-fluence.
        "referenceFluence": None,
        "referenceNote": "supplied on the command line",
        # 100 um holds 99 files: tasks 50 and 100 drew the same autoSeed pair
        # (931749, 144148) and wrote the same TID_931749_144148.csv. Identical
        # seeds are identical histories, so the survivor is a valid sample and
        # the run is only 1% short.
        "runs": {
            f"{epoxy} um epoxy": {"folder": f"10MeVProton-CRDS-{epoxy}umEpoxy",
                                  "job": job, "files": files,
                                  "field": "25 x 25", "epoxy_mm": epoxy / 1000}
            for epoxy, job, files in [(800, 20396958, 100), (700, 20396957, 100),
                                      (600, 20396956, 100),
                                      (500, 20395871, 100), (400, 20395872, 100),
                                      (300, 20395873, 100), (200, 20395874, 100),
                                      (100, 20395875, 99)]
        },
        # One numeric parameter varies across the runs, so the series driver
        # also plots every tally against it.
        "sweep": {"key": "epoxy_mm", "xlabel": "Epoxy above the GaAs die [mm]",
                  "output": "EpoxySweep_CRDS-Kumpula.pdf"},
        # Each thickness against 500 um, the original upper bracket. The
        # RadFET package was not modified, so its TID must not move across the
        # sweep; the LED is what the sweep is for.
        "steps": [(f"{epoxy} um epoxy", "500 um epoxy")
                  for epoxy in (800, 700, 600, 400, 300, 200, 100)],
        "seriesOutput": "EpoxySweep_CRDS-Kumpula.csv",
        "oxideComparison": None,
        "activeSubstrateNote": (
            "These need not agree: a 10 MeV proton stops in or near the die,\n"
            "  so the active layer and the substrate sample different parts "
            "of the\n  Bragg curve. Track this ratio across the epoxy sweep."),
    },
    "Uppsala": {
        "path": CRDS_PATH / "CRDS-Uppsala",
        "description": "Uppsala Skandion 62 MeV proton beam, two scanned fields",
        "particle": "proton",
        "defaultFolder": "64MeVProton-70x74mmField",
        # No single reference fluence: the delivered fluence differs by part
        # and by exposure, and the ``delivery`` block below accounts for it.
        "referenceFluence": None,
        "referenceNote": "supplied on the command line",
        "runs": {
            "52 x 54 mm": {"folder": "64MeVProton-52x54mmField",
                           "job": 20399188, "files": 100, "field": "52 x 54"},
            "70 x 74 mm": {"folder": "64MeVProton-70x74mmField",
                           "job": 20399189, "files": 200, "field": "70 x 74"},
        },
        # Field size at fixed geometry. Only the RadFET is a coefficient in
        # both runs, so this is the test of the assumption the LED relies on.
        "steps": [("52 x 54 mm", "70 x 74 mm")],
        "seriesOutput": "FieldSizeControl_CRDS-Uppsala.csv",
        "sweep": None,
        # The LED sits at y = +27.85 mm, outside the sharp simulated edge of
        # the 52 x 54 mm field at +27, so its tally there is mostly scatter.
        "notCoefficients": {("52 x 54 mm", "DDD")},
        "oxideComparison": None,
        "activeSubstrateNote": (
            "These should agree: a 62 MeV proton loses under 2 MeV crossing\n"
            "  0.18 mm of GaAs, so the field is uniform through the die."),
        "delivery": {
            "maps": UPPSALA_LOGS,
            # Design-file part centres in GDML coordinates, the same numbers
            # as RadFET_pos_* and LED_pos_* in CRDS-Uppsala-define.xml.
            # CHANGE THEM TOGETHER.
            "parts": {"RadFET": (-37 + 21.57, 35 - 18.40),
                      "LED": (-37 + 20.78, 35 - 7.15)},
            "partOfQuantity": {"TID": "RadFET", "DDD": "LED"},
            # Board-on-field alignment was not recorded; report the spread
            # over a board shift of this size in x and in y.
            "alignmentShift_mm": 1.0,
            # The LED takes the 70 x 74 mm coefficient for both exposures,
            # because it lies in the soft edge of the narrow field.
            "exposures": [
                {"label": "exp_10-12, 5.2 x 5.4 cm", "map": "CRDS_5.2x5.4cm",
                 "coefficientRun": {"RadFET": "52 x 54 mm", "LED": "70 x 74 mm"}},
                {"label": "exp_13, 7.0 x 7.4 cm", "map": "CRDS_7.0x7.4cm",
                 "coefficientRun": {"RadFET": "70 x 74 mm", "LED": "70 x 74 mm"}},
            ],
            "output": "DeliveredDose_CRDS-Uppsala.csv",
        },
    },
}


def campaign(name):
    """Look up a campaign by name, naming the valid ones if it is unknown."""
    if name not in CAMPAIGNS:
        raise KeyError(f"unknown CRDS campaign {name!r}; "
                       f"known: {', '.join(CAMPAIGNS)}")
    return CAMPAIGNS[name]


def runForFolder(spec, folder):
    """Return (label, run) for a run folder, or (None, None) if it is unlisted.

    Pilots and one-off controls need not be listed; the single-run drivers
    then require the job ID on the command line.
    """
    for label, run in spec["runs"].items():
        if run["folder"] == folder:
            return label, run
    return None, None

