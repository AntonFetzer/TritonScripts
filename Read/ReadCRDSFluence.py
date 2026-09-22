"""Read CRDS accelerator-log fluence maps and summary values."""

import json

import numpy as np

from UppsalaAcceleratorLogs.radex_isocentre_outputs import bilinear_sample


def sampleFluence(mapPath, x, y):
    """Fluence in particles/cm2 at (x, y) mm on one isocentre heatmap."""
    with np.load(mapPath) as data:
        return bilinear_sample(data["fluence_protons_cm2"], data["x_edges_mm"],
                               data["y_edges_mm"], x, y)


def fieldAverageFluence(summaryPath):
    """Delivered protons over the requested target area, from summary.json."""
    with open(summaryPath, encoding="utf-8") as stream:
        return float(json.load(stream)["crds_mean_fluence_raw_protons_cm2"])
