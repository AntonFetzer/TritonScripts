"""Overlay GRAS histograms from several runs, for any histogram kind.

This replaces six scripts that were the same script six times:
CompareDoseHistograms, CompareLETHistograms, CompareFluenceHistograms,
CompareSourceHistograms, StackedDoseHistograms and StackedLETHistograms. Each
read a list of result directories with the reader matching its histogram kind,
looped over them drawing the same bar-plus-step or bar-plus-errorbar idiom, and
saved a log-log PDF. They differed only in the reader, the paths and the axis
labels, because all four readers return the same histogram dictionary.

The drawing idioms now live in Dependencies.HistogramPlots; this script selects
a reader, assembles the list of runs, and labels the axes.

    # explicit runs, the old "Compare" mode
    python3 Plotting/CompareHistograms.py --kind let --component let \
        --paths /path/A/Res /path/B/Res --labels "16 mm Al" "8 mm Al" \
        --output /path/LET-comparison.pdf

    # every subfolder of a parent, the old "Stacked" mode
    python3 Plotting/CompareHistograms.py --kind let --component let \
        --scan /path/Protons16mmAl-200micronSi --suffix MeV \
        --style errorbar --output /path/LET-Values.pdf

Paths are arguments rather than constants. The six scripts this replaces
hardcoded lists under ``/l/triton_work`` and ``/l/TritonPlots``, which are not
present on Triton any more: ``/l`` holds only ``lost+found``. Those lists are
recorded in the scripts' git history; none of them is reproduced here, because
a hardcoded path that cannot resolve is worse than no default.

Run from the repository root so Dependencies.* resolves.
"""

import argparse
import os
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

import matplotlib.pyplot as plt  # noqa: E402

from Dependencies.HistogramPlots import (  # noqa: E402
    colour, errorbarHistogram, logLogGrid, savePdf, stepHistogram)

# Readers are imported inside their component function, not at module level.
# They do not all have the same third-party dependencies -- the source reader
# needs pandas, which is absent from the Triton login-node interpreter -- and
# an eager import would make every kind fail because of one kind's dependency.


def doseComponents(path):
    """Dose histograms: by deposited dose, and by primary kinetic energy."""
    from Dependencies.TotalDoseHistograms import totalDoseHistograms
    dose, primary = totalDoseHistograms(path)
    return {"dose": dose, "primary": primary}


def letComponents(path):
    """LET histograms: the LET spectrum and the effective-LET spectrum."""
    from Dependencies.TotalLETHistos import totalLETHistos
    let, effective = totalLETHistos(path)
    return {"let": let, "eff": effective}


def fluenceComponents(path):
    """Fluence histograms, one per species."""
    from Dependencies.TotalFluenceHistograms import totalFluenceHistos
    electron, proton = totalFluenceHistos(path)
    return {"electron": electron, "proton": proton}


def sourceComponents(path):
    """Source histograms; the reader already returns a dict keyed by quantity.

    Needs pandas, which the Triton login-node interpreter does not have.
    """
    from Dependencies.TotalSourceHistograms import totalSourceHistos
    histograms = totalSourceHistos(path)
    return {name.lower(): histogram for name, histogram in histograms.items()}


# Per-kind reader and default axis labels. Labels are defaults only; --xlabel
# and --ylabel override them, because the right label depends on how the run
# was normalised and this script cannot know that.
KINDS = {
    "dose": {
        "reader": doseComponents,
        "default": "primary",
        "xlabel": "Kinetic energy [MeV]",
        "ylabel": "Relative ionising dose [a.u.]",
    },
    "let": {
        "reader": letComponents,
        "default": "let",
        "xlabel": "LET [MeV cm2 mg-1]",
        "ylabel": "Rate per LET bin [s-1]",
    },
    "fluence": {
        "reader": fluenceComponents,
        "default": "electron",
        "xlabel": "Kinetic energy [MeV]",
        "ylabel": "Fluence per bin [cm-2]",
    },
    "source": {
        "reader": sourceComponents,
        "default": "energy",
        "xlabel": "Energy [MeV]",
        "ylabel": "Counts per bin",
    },
}

# 'entries' is a raw count, so a label about rates or dose would be wrong.
ENTRIES_YLABEL = "Number of entries per bin"


def collectRuns(arguments):
    """Build the ordered list of (resultPath, label) pairs to overlay.

    Two ways in, matching the two families of script this replaces: an explicit
    list of result directories, or a parent directory whose subfolders each
    hold a Res/.
    """
    if arguments.paths:
        paths = [Path(path) for path in arguments.paths]
    else:
        parent = Path(arguments.scan)
        if not parent.is_dir():
            raise SystemExit(f"--scan directory not found: {parent}")
        folders = sorted(
            (child for child in parent.iterdir()
             if child.is_dir()
             and (arguments.suffix is None or child.name.endswith(arguments.suffix))
             and (child / "Res").is_dir()),
            key=lambda child: naturalKey(child.name),
            reverse=arguments.reverse,
        )
        if not folders:
            raise SystemExit(
                f"no subfolder of {parent} has a Res/ directory"
                + (f" and ends in {arguments.suffix!r}" if arguments.suffix else ""))
        paths = [folder / "Res" for folder in folders]

    if arguments.labels:
        if len(arguments.labels) != len(paths):
            raise SystemExit(
                f"--labels has {len(arguments.labels)} entries for {len(paths)} runs")
        labels = list(arguments.labels)
    else:
        # The run directory is the parent of Res/, which is what names the run.
        labels = [os.path.basename(os.path.dirname(str(path).rstrip("/")))
                  for path in paths]
    return list(zip(paths, labels))


def naturalKey(name):
    """Sort key that orders '8MeV' before '16MeV' rather than after it."""
    digits = "".join(character for character in name if character.isdigit())
    return (int(digits) if digits else 0, name)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--kind", required=True, choices=sorted(KINDS))
    parser.add_argument("--component",
                        help="which histogram of that kind; default depends on --kind")
    parser.add_argument("--field", default="value", choices=["value", "entries"],
                        help="per-bin quantity to draw (default: value)")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--paths", nargs="+", help="result directories to overlay")
    source.add_argument("--scan", help="parent directory whose subfolders hold Res/")
    parser.add_argument("--suffix", help="with --scan, only subfolders ending in this")
    parser.add_argument("--reverse", action="store_true",
                        help="with --scan, order subfolders descending")
    parser.add_argument("--labels", nargs="+", help="legend entries, one per run")
    parser.add_argument("--style", default="step", choices=["step", "errorbar"])
    parser.add_argument("--title")
    parser.add_argument("--xlabel")
    parser.add_argument("--ylabel")
    parser.add_argument("--xscale", default="log", choices=["log", "linear"])
    parser.add_argument("--yscale", default="log", choices=["log", "linear"])
    parser.add_argument("--output", type=Path, help="write a PDF here")
    parser.add_argument("--show", action="store_true", help="open the figure")
    arguments = parser.parse_args()

    if arguments.suffix and not arguments.scan:
        parser.error("--suffix only applies with --scan")
    if arguments.reverse and not arguments.scan:
        parser.error("--reverse only applies with --scan")
    if not arguments.output and not arguments.show:
        parser.error("nothing to do: pass --output, --show, or both")

    kind = KINDS[arguments.kind]
    component = arguments.component or kind["default"]
    runs = collectRuns(arguments)

    plt.figure()
    drawn = 0
    for index, (path, label) in enumerate(runs):
        components = kind["reader"](str(path))
        if component not in components:
            raise SystemExit(
                f"--component {component!r} is not one of "
                f"{sorted(components)} for --kind {arguments.kind}")
        histogram = components[component]
        if histogram is None:
            # The readers return None for a directory they could not use.
            print(f"WARNING: {path} returned no {component!r} histogram; skipped")
            continue
        draw = stepHistogram if arguments.style == "step" else errorbarHistogram
        draw(histogram, field=arguments.field, label=label, color=colour(index))
        drawn += 1

    if not drawn:
        raise SystemExit("no histogram could be drawn from any of the given runs")

    ylabel = arguments.ylabel
    if ylabel is None:
        ylabel = ENTRIES_YLABEL if arguments.field == "entries" else kind["ylabel"]
    logLogGrid(
        title=arguments.title or f"{arguments.kind} / {component}",
        xlabel=arguments.xlabel or kind["xlabel"],
        ylabel=ylabel,
        xscale=arguments.xscale,
        yscale=arguments.yscale,
    )

    if arguments.output:
        print(f"Wrote {savePdf(arguments.output)}")
    if arguments.show:
        plt.show()


if __name__ == "__main__":
    main()
