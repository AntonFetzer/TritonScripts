# Agent instructions for TritonScripts

## Repository purpose

This repository contains Python tools for analysing particle-radiation simulations and measurements, including GRAS/Geant4 dose, LET, fluence, shielding, and single-event-effect results. Treat numerical conventions, units, calibration factors, geometry mappings, and uncertainty propagation as scientifically significant. Do not change them silently.

## Project layout

- `Dependencies/`: aggregation and processing utilities shared by analysis scripts.
- `Read/`: parsers for GRAS and external radiation-data formats.
- `Carrington/`: Carrington-event spectra, shielding, dose, and SEE studies.
- `SEE/`: single-event-effect cross sections, fits, and rate calculations.
- `Plotting/`: plotting and comparison scripts.
- `Calibration/`, `SolarCells/`, and other study folders: domain-specific analyses.
- `Obsolete/`: retired implementations kept only for historical reference. Do not import from this directory in active code.

## Current data interfaces

Use dictionary-based data structures. Do not reintroduce positional column IDs or legacy dose arrays.

```python
from Dependencies.TotalDose import totalDose

result = totalDose(results_path)
dose = result["dose"]
error = result["error"]
entries = result["entries"]
non_zeros = result["non-zeros"]
```

`totalDose` accepts one GRAS result dataset directory. For a legacy directory containing multiple datasets, select one explicitly with `filename_contains`.

Histogram dictionaries use these fields unless the reader documents additional fields:

```python
required = ("lower", "upper", "mean", "value", "error", "entries")
```

Important return signatures:

```python
from Dependencies.TotalLETHistos import totalLETHistos
let_histogram, effective_histogram = totalLETHistos(path)

from Dependencies.TotalFluenceHistos import totalFluenceHistos
electron_histogram, proton_histogram = totalFluenceHistos(path)
```

Always unpack these tuples explicitly. LET readers already return the repository's documented converted units; do not apply a second density conversion in callers.

## Implementation conventions

- Run scripts from the repository root when they use absolute `Dependencies.*` or `Read.*` imports.
- Prefer `pathlib.Path` or `os.path` over constructing new paths with string concatenation.
- Derive array lengths, exponents, channel counts, and bin counts from input data when possible.
- Validate required dictionary keys and compatible shapes before combining datasets.
- Combine independent absolute uncertainties in quadrature unless the analysis documents a different statistical model.
- Preserve established plot appearance unless a task explicitly requests a style change.
- Keep generated data and plots outside this source repository unless they are intentional fixtures or documentation assets.
- Do not modify unrelated files in a dirty working tree.

## Scientific checks

- State or document assumptions that affect fluence, dose, detector position, particle energy, shielding, calibration, or units.
- Check for NaN and infinite values after divisions, fits, interpolations, and normalisation.
- Confirm histogram binning compatibility before summing histograms.
- Preserve raw measurements; write derived results to separate outputs.
- Do not infer missing facility or calibration parameters when they materially affect a reported result. Surface the uncertainty instead.

## Validation

After Python edits, run at minimum:

```bash
python3 -m compileall -q .
git diff --check
```

Also run the narrowest relevant analysis or test when its required input data is available. Report when end-to-end validation cannot be performed because external GRAS or measurement data is unavailable.

On Aalto Triton, keep login-node tasks below 120 CPU seconds and strictly single-threaded.
Use Slurm for longer or parallel computations.
Minimum job duration is 120s.
Maximum number of CPUs to be used concurrently 100.
Check how many CPUs are already in use before starting new jobs.
