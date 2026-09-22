"""Estimate CRDS isocentre maps using the RadEx PLD geometric calibration."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from UppsalaAcceleratorLogs.register_fluence_maps import load_heatmap, registered_full_grid
from UppsalaAcceleratorLogs.radex_isocentre_outputs import RADEX_CHAMBER_TO_ISOCENTRE

GROUPS = ("CRDS_5.2x5.4cm", "CRDS_7.0x7.4cm")
ASSUMPTION = (
    "The 85 MeV RadEx PLD geometric scales are assumed transferable to 64 MeV CRDS. "
    "Each target-size aggregate is recentered by its fluence-weighted chamber centroid. "
    "The sum assumes both target-size fields centred on the same isocentre; actual DUT "
    "translation and rotation are unknown. No CRDS PLD or independent spatial validation "
    "is available. Scale method sensitivity is about +/-0.02 in X and +/-0.03 in Y, "
    "not a statistical confidence interval. Input chamber spot-width modelling is retained."
)

def _load(path):
    values, x, y = load_heatmap(path, ("fluence_protons_cm2",), -250., 250.)
    if not np.all(np.isfinite(values)) or np.any(values < 0) or values.sum() <= 0:
        raise ValueError(f"Invalid fluence in {path}")
    if not np.allclose(np.diff(x), np.diff(x)[0]) or not np.allclose(np.diff(y), np.diff(y)[0]):
        raise ValueError("Registration requires uniform bins")
    return values, x, y

def _integral(values, x, y):
    return float(np.sum(values * np.diff(y)[:, None] * np.diff(x)[None, :]) / 100.)

def _save(folder, values, x, y, summary, dimensions):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("Invalid transformed fluence")
    integral = _integral(values, x, y)
    summary.update(isocentre_integrated_protons=integral,
                   integral_ratio=integral / summary["input_integrated_protons"],
                   coordinate_assumption=ASSUMPTION)
    if abs(summary["integral_ratio"] - 1.) > 1e-4:
        raise ValueError(f"Proton conservation failed: {summary}")
    np.savez_compressed(folder / "heatmap_isocentre.npz",
                        fluence_protons_cm2=values,
                        x_edges_mm=x, y_edges_mm=y,
                        metadata_json=json.dumps(summary))
    (folder / "isocentre_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)
    im = ax.imshow(values, origin="lower", extent=(x[0], x[-1], y[0], y[-1]), cmap="magma", vmin=0)
    for i, (width, height) in enumerate(dimensions):
        ax.add_patch(Rectangle((-width * 5, -height * 5), width * 10, height * 10,
                     fill=False, edgecolor=("limegreen", "cyan")[i], linestyle=("--", ":")[i],
                     label=f"Requested target {width:g} x {height:g} cm"))
    ax.set(xlim=(-65, 65), ylim=(-65, 65), xlabel="Estimated isocentre X [mm]",
           ylabel="Estimated isocentre Y [mm]", title=folder.name + " — estimated isocentre fluence")
    ax.legend(loc="upper right", fontsize=8)
    fig.colorbar(im, ax=ax, label="Delivered protons cm$^{-2}$")
    fig.text(.5, .005, "RadEx-derived scaling; fields assumed centred; DUT pose unknown",
             ha="center", fontsize=8)
    fig.savefig(folder / "fluence_isocentre.png", dpi=180)
    plt.close(fig)
    print(f'{folder.name}: integral ratio={summary["integral_ratio"]:.10f}; protons={integral:.8e}')

def main():
    """Generate both CRDS field-size products and their aligned cumulative sum."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aggregates", type=Path, required=True)
    args = parser.parse_args()
    transformed, originals, summaries, dimensions = [], [], [], []
    common_x = common_y = None
    for name in GROUPS:
        folder = args.aggregates / name
        values, x, y = _load(folder / "heatmap.npz")
        if common_x is not None and (not np.array_equal(x, common_x) or not np.array_equal(y, common_y)):
            raise ValueError("CRDS aggregate grids differ")
        common_x, common_y = x, y
        cx = float(np.dot(values.sum(axis=0), (x[:-1] + x[1:]) / 2) / values.sum())
        cy = float(np.dot(values.sum(axis=1), (y[:-1] + y[1:]) / 2) / values.sum())
        transform = {key: val for key, val in RADEX_CHAMBER_TO_ISOCENTRE.items()
                     if key != "pearson_correlation"}
        transform["offset_x_isocentre_mm"] = -transform["scale_x_chamber_to_isocentre"] * cx
        transform["offset_y_isocentre_mm"] = -transform["scale_y_chamber_to_isocentre"] * cy
        result = registered_full_grid(values, x, y, x, y, transform, normalize_input=False)
        result *= transform["fluence_density_jacobian_chamber_to_isocentre"]
        source_summary = json.loads((folder / "summary.json").read_text())
        dims = source_summary["crds_target_dimensions_cm"]
        summary = dict(source_heatmap=str(folder / "heatmap.npz"), transform=transform,
                       chamber_centroid_mm=[cx, cy], input_integrated_protons=_integral(values, x, y),
                       experiment_ids=source_summary["experiment_ids"],
                       calibration_source="metadata/RA_RadEx_setup.pld; Radex_85MeV registration",
                       spatialized_to_delivered_protons_ratio=source_summary["spatialized_to_delivered_protons_ratio"])
        _save(folder, result, x, y, summary, dims)
        transformed.append(result)
        originals.append(values)
        summaries.append(summary)
        dimensions.extend(dims)
    folder = args.aggregates / "CRDS_Total"
    total, tx, ty = _load(folder / "heatmap.npz")
    if not np.array_equal(tx, common_x) or not np.array_equal(ty, common_y):
        raise ValueError("Total grid differs")
    if not np.allclose(total, sum(originals), rtol=2e-6, atol=1e-6):
        raise ValueError("Original CRDS total does not equal its components")
    summary = dict(source_heatmap=str(folder / "heatmap.npz"),
                   construction="Sum of independently recentered target-size isocentre maps",
                   components=summaries, input_integrated_protons=_integral(total, tx, ty))
    _save(folder, sum(transformed), tx, ty, summary, dimensions)

if __name__ == "__main__":
    main()
