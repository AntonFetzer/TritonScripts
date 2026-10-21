"""Figures for kernel fits of measured per-channel doses: doses with model predictions, pulls, chi2 scans."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plotChannelFit(channels, measured, sigma, predictions, path, groupOf=None, title=None):
    """Measured doses with error bars and every model's prediction, plus the best model's pulls.

    Args:
        channels (sequence): channel labels, in order.
        measured, sigma (arrays): measured dose and its uncertainty.
        predictions (list): (values, label, pulls or None) per model; the first
            model's pulls fill the lower panel when given.
        path (str): output PDF.
        groupOf (sequence): a group name per channel; VT05-style groups get a
            square marker, others a circle.
        title (str): upper panel title.
    """
    x = np.arange(len(channels))
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(10, 8.5), gridspec_kw={"height_ratios": [3, 1]}, sharex=True)
    groups = list(groupOf) if groupOf is not None else [""] * len(channels)
    for group in dict.fromkeys(groups):
        sel = np.array([g == group for g in groups])
        y, s = measured[sel], sigma[sel]
        top.errorbar(x[sel], y, yerr=[np.minimum(s, 0.95 * y), s], fmt="s" if group == "VT05" else "o",
                     color="black", capsize=3, label=f"measured {group}".strip())
    for i, (values, label, _) in enumerate(predictions):
        top.plot(x, np.where(values > 0, values, np.nan), "_", color=f"C{i}", ms=16, mew=2, label=label)
    top.set_yscale("log")
    top.set_ylabel("dose (krad)")
    top.set_title(title or "")
    top.grid(alpha=0.25, which="both")
    top.legend(fontsize=8)
    for i, (_, label, pulls) in enumerate(predictions):
        if pulls is not None:
            bottom.bar(x + (i - 0.5 * (len(predictions) - 1)) * 0.8 / len(predictions), pulls,
                       0.8 / len(predictions), color=f"C{i}", label=label.split(":")[0])
    for level in (-2, 0, 2):
        bottom.axhline(level, color="grey" if level else "black", ls="--" if level else "-", lw=0.8)
    bottom.set_xticks(x, [str(c) for c in channels])
    bottom.set_xlabel("channel")
    bottom.set_ylabel("pull (model - measured) / sigma")
    bottom.grid(alpha=0.25, axis="y")
    bottom.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plotChi2Scans(scans, path, ylim=(0, 25), title=None):
    """Delta-chi2 curves, one panel per scan.

    Args:
        scans (list): (x, curves, xlabel, marks) per panel, where curves is a
            list of (chi2, label) drawn relative to each curve's own minimum and
            marks a list of (x value, label) drawn as vertical lines.
        path (str): output PDF.
    """
    fig, axes = plt.subplots(1, len(scans), figsize=(6.2 * len(scans), 4.5), squeeze=False)
    for ax, (x, curves, xlabel, marks) in zip(axes[0], scans):
        for i, (chi2, label) in enumerate(curves):
            ax.plot(x, chi2 - np.nanmin(chi2), "-", color=f"C{i}", label=label)
        for value, label in marks:
            ax.axvline(value, color="C3", lw=0.8, ls=":", label=label)
        ax.axhline(1.0, color="grey", lw=0.8)
        ax.set_ylim(*ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("delta chi2 to each curve's minimum")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
