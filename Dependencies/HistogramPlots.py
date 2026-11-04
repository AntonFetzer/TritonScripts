"""Draw the repository's standard histogram dictionaries.

Every histogram reader in ``Read/`` and ``Dependencies/`` returns the same
dictionary contract, documented in AGENTS.md::

    ("lower", "upper", "mean", "value", "error", "entries")

Because the contract is shared, so was the drawing code: six scripts in
``Plotting/`` each carried their own copy of the same two idioms, a filled bar
plus a step outline, or a translucent bar plus error bars at the bin means,
followed by the same log-log grid and the same PDF save. They differed only in
which reader produced the histogram and what the axis labels said. Those idioms
live here now, so a change to how histograms look happens once.

Nothing here reads a GRAS file or combines datasets; it only draws what a
reader returned. Bin values are taken as given, including their errors.
"""

import matplotlib.pyplot as plt
import numpy as np

# The dictionary contract every histogram reader in this repository returns.
REQUIRED_FIELDS = ("lower", "upper", "mean", "value", "error", "entries")

# Default colour cycle, carried over from the scripts this module replaces so
# that regenerated figures keep their established appearance.
COLOURS = ["C1", "C0", "C2", "C8", "C3", "C7", "C9", "C4", "C5", "C6"]


def checkHistogram(histogram, field="value"):
    """Validate a histogram dictionary before drawing it.

    Args:
        histogram (dict): a reader's output.
        field (str): the per-bin quantity about to be drawn, usually 'value'
            or 'entries'.

    Returns:
        dict: the same histogram, unchanged.

    Raises:
        TypeError: the histogram is not a dictionary.
        KeyError: a required field or the requested field is missing.
        ValueError: the arrays are not all the same length, or a bin edge is
            not finite.
    """
    if not isinstance(histogram, dict):
        raise TypeError(f"histogram must be a dict, got {type(histogram).__name__}")

    missing = [key for key in REQUIRED_FIELDS if key not in histogram]
    if missing:
        raise KeyError(f"histogram is missing {missing}; expected {list(REQUIRED_FIELDS)}")
    if field not in histogram:
        raise KeyError(f"histogram has no field {field!r}")

    lengths = {key: len(histogram[key]) for key in REQUIRED_FIELDS}
    if len(set(lengths.values())) != 1:
        raise ValueError(f"histogram arrays have mismatched lengths: {lengths}")

    edges = np.concatenate([np.asarray(histogram["lower"], dtype=float),
                            np.asarray(histogram["upper"], dtype=float)])
    if not np.isfinite(edges).all():
        raise ValueError("histogram has non-finite bin edges")
    return histogram


def nonZeroRange(histograms, field="value"):
    """Energy range spanned by the bins with a non-zero value, over several histograms.

    Matplotlib's autoscaling covers every bin, because a bin holding zero is
    valid data, so a spectrum that fills two decades of an eight-decade
    histogram ends up squeezed into a corner. This returns the limits that fit
    the populated bins instead.

    Args:
        histograms (iterable of dict): readers' outputs.
        field (str): the per-bin quantity that decides whether a bin counts,
            usually the one being drawn.

    Returns:
        tuple(float, float) or None: the lower edge of the first and the upper
        edge of the last bin whose ``field`` is non-zero in any histogram, or
        None when every bin of every histogram is zero.
    """
    lowest, highest = np.inf, -np.inf
    for histogram in histograms:
        checkHistogram(histogram, field)
        filled = np.nonzero(np.asarray(histogram[field], dtype=float))[0]
        if filled.size:
            lowest = min(lowest, float(np.asarray(histogram["lower"], dtype=float)[filled].min()))
            highest = max(highest, float(np.asarray(histogram["upper"], dtype=float)[filled].max()))
    if not np.isfinite(lowest):
        return None
    return lowest, highest


def colour(index):
    """Colour for series ``index``, cycling when there are more series than colours."""
    return COLOURS[index % len(COLOURS)]


def stepHistogram(histogram, field="value", label=None, color=None, alpha=0.5,
                  axes=None, fill=True):
    """Draw a histogram as a filled bar with a step outline.

    This is the comparison idiom: several runs overlaid, each a translucent
    fill so overlaps stay readable, with a hard step outline carrying the
    legend entry.

    Args:
        histogram (dict): a reader's output.
        field (str): per-bin quantity to draw, 'value' or 'entries'.
        label (str or None): legend entry; None omits this series from the legend.
        color (str or None): matplotlib colour; None lets matplotlib choose.
        alpha (float): opacity of the fill. The outline is always opaque.
        axes (matplotlib.axes.Axes or None): target axes; None uses the current axes.
        fill (bool): draw the translucent bar under the outline; False draws
            the outline alone, which keeps many overlaid runs readable.

    Returns:
        matplotlib.axes.Axes: the axes drawn on.
    """
    checkHistogram(histogram, field)
    axes = axes or plt.gca()
    lower = np.asarray(histogram["lower"], dtype=float)
    upper = np.asarray(histogram["upper"], dtype=float)
    values = np.asarray(histogram[field], dtype=float)

    if fill:
        axes.bar(lower, values, width=upper - lower, align="edge", alpha=alpha,
                 color=color)
    axes.step(lower, values, where="post", label=label, color=color)
    return axes


def errorbarHistogram(histogram, field="value", label=None, color=None,
                      alpha=0.3, axes=None, fill=True, capsize=5):
    """Draw a histogram as a translucent bar with error bars at the bin means.

    This is the idiom for showing a single series whose per-bin uncertainty
    matters, rather than comparing shapes between runs.

    Args:
        histogram (dict): a reader's output.
        field (str): per-bin quantity to draw, 'value' or 'entries'.
        label (str or None): legend entry.
        color (str or None): matplotlib colour.
        alpha (float): opacity of the fill.
        axes (matplotlib.axes.Axes or None): target axes.
        fill (bool): draw the underlying bar as well as the error bars.
        capsize (float): size of the error bar caps.

    Returns:
        matplotlib.axes.Axes: the axes drawn on.
    """
    checkHistogram(histogram, field)
    axes = axes or plt.gca()
    lower = np.asarray(histogram["lower"], dtype=float)
    upper = np.asarray(histogram["upper"], dtype=float)
    mean = np.asarray(histogram["mean"], dtype=float)
    values = np.asarray(histogram[field], dtype=float)
    errors = np.asarray(histogram["error"], dtype=float)

    if fill:
        axes.bar(lower, values, width=upper - lower, align="edge", alpha=alpha,
                 color=color)
    axes.errorbar(mean, values, errors, fmt=" ", capsize=capsize, elinewidth=1,
                  capthick=1, label=label, color=color)
    return axes


def logLogGrid(axes=None, title=None, xlabel=None, ylabel=None, legend=True,
               xscale="log", yscale="log"):
    """Apply the shared axis styling used by every histogram figure here.

    Args:
        axes (matplotlib.axes.Axes or None): target axes.
        title, xlabel, ylabel (str or None): applied when not None.
        legend (bool): draw a legend, but only if something carries a label.
        xscale, yscale (str): matplotlib scale names. Pass 'linear' to override
            one of them; relative-error axes in this repository are linear from
            zero by convention.

    Returns:
        matplotlib.axes.Axes: the axes styled.
    """
    axes = axes or plt.gca()
    axes.set_xscale(xscale)
    axes.set_yscale(yscale)
    axes.grid(which="both")
    if title is not None:
        axes.set_title(title)
    if xlabel is not None:
        axes.set_xlabel(xlabel)
    if ylabel is not None:
        axes.set_ylabel(ylabel)
    # Calling legend() with nothing labelled emits a matplotlib warning and an
    # empty box, which happens whenever a caller passes label=None throughout.
    handles, labels = axes.get_legend_handles_labels()
    if legend and labels:
        axes.legend()
    return axes


def savePdf(path, figure=None):
    """Save a figure as PDF with the tight bounding box used throughout.

    Args:
        path (str or Path): destination file.
        figure (matplotlib.figure.Figure or None): figure to save; None uses
            the current figure.

    Returns:
        str: the path written, as a string.
    """
    figure = figure or plt.gcf()
    figure.savefig(str(path), format="pdf", bbox_inches="tight")
    return str(path)


def histogramWithCrossSection(histogram, sigma, path, title, xlabel, ylabel, label,
                              sigmaLabel, sigmaAxisLabel, legendLoc=None, capsize=5):
    """Draw a histogram with an SEU cross section on a second axis and save it as PDF.

    The figure of the SEU rate scripts in SEE/: the histogram on the left axis in
    C0 and the cross section on the right axis in C1, both logarithmic.

    Args:
        histogram (dict): a reader's output.
        sigma (numpy.ndarray): the cross section at the lower edges of the bins.
        path (str or Path): destination file.
        title, xlabel, ylabel (str): labels of the figure and of the left axis.
        label (str): legend entry of the histogram.
        sigmaLabel (str): legend entry of the cross section.
        sigmaAxisLabel (str): label of the right axis.
        legendLoc (str or None): matplotlib legend location, None for the default.
        capsize (float): size of the error bar caps.

    Returns:
        str: the path written, as a string.
    """
    figure, ax1 = plt.subplots(1)
    errorbarHistogram(histogram, label=label, axes=ax1, capsize=capsize)
    ax1.plot([], [], label=sigmaLabel, color="C1")
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.grid(which="major")
    ax1.set_title(title)
    ax1.set_xlabel(xlabel)
    ax1.legend(loc=legendLoc)
    ax1.set_ylabel(ylabel, color="C0")
    ax1.tick_params(axis="y", colors="C0")

    ax2 = ax1.twinx()
    ax2.plot(histogram["lower"], sigma, color="C1")
    ax2.set_ylabel(sigmaAxisLabel, color="C1")
    ax2.set_yscale("log")
    ax2.tick_params(axis="y", colors="C1")

    written = savePdf(path, figure)
    plt.close(figure)
    return written
