# Uppsala proton-log analysis


Campaign inputs, caches, logs, and results remain under /scratch/work/fetzera1/Radiation Testing/2025-06 Proton/Uppsala Accelerator Logs/.
The version-controlled package contains a streaming parser and Slurm workflow for the Skandion
Kliniken IBA PBS data-recorder logs from 27–28 June 2025.

The input logs are never modified. Parsing first creates charge-domain products,
so calibration and plotting can be repeated without rescanning the 14.4 GiB of
CSV input.

## Campaign layout

- `raw/`: immutable recorder CSVs.
- `metadata/`: experiment log, conversion table, and file manifest.
- `cache/parsed/`: regenerable per-record charge products.
- `results/final/`: calibrated production results.
- `results/exploratory/`: cadence and time-spatial studies.
- `logs/slurm/`: retained scheduler logs.

## Workflow

1. Build `metadata/manifest.json` with `python3 -m UppsalaAcceleratorLogs.uppsala_analysis inventory` from the Python repository root.
2. Optionally run a small parse pilot with the `--session-filter` and `--limit` options.
3. Inspect charge conservation and delivered-versus-planned calibration checks.
4. From the campaign root, submit `/scratch/work/fetzera1/Python/UppsalaAcceleratorLogs/slurm/production_array.sbatch` (four byte-balanced, multi-minute tasks).
5. Submit `/scratch/work/fetzera1/Python/UppsalaAcceleratorLogs/slurm/production_reduce.sbatch` with an `afterok` dependency.

The facility file's extra values (`7.14E12` and `9.40E12`) equal the Experiment 1
MU totals multiplied by the 65 MeV and 85 MeV `Protons per MU` rows. Accordingly,
the campaign's nominal 64 MeV setting uses the facility's 65 MeV value
(`5.71E7 protons/MU`) without interpolation; 85 MeV uses `6.97E7 protons/MU`.

Delivered MU is reconstructed independently for every `part` record from the
logged incremental primary-monitor charge:

`MU = sum(DOSE_PRIM(C)) / (3.0E-9 C/MU * K_FACTOR)`.

The corresponding facility protons/MU value converts each record to protons
before time-series and spatial aggregation. `TOTAL_CHARGE` from `map_specif`
files is retained only as a QC field. It must not be summed as prescribed charge:
canceled, partial, and restarted attempts each repeat the full nominal
`TOTAL_CHARGE`, which would create a false apparent underdelivery.

All normally completed irradiations agree with their prescribed MU within a
fraction of a percent. `exp_09` is an intentionally early-canceled run; its
10,000 MU entry was only an estimate, while the primary monitor recorded
10,057.49 MU-equivalent.

One `exp_09` map-record part contains two complete numerical stream sections
covering the same time interval, separated by cancellation metadata. These are
two serializations of the same layer, not two delivered layers: the second is
slightly more complete (74,868 rows) and has the same measured charge and scan
trajectory. The parser therefore retains only the last numerical section of a
map-record. This removes the former artificial factor-of-two rate at the end of
`RadEx_Total` and reduces its delivered total from about 1.6679E13 to
1.6556E13 protons. It does not change the preceding irradiations.

## Principal outputs

- fluence_vs_time_1s.csv: one-second charge and delivered-particle series.
- fluence_vs_time_1s_RadEx_area.csv and fluence_vs_time_RadEx_area.png:
  mean incident fluence across the specified RadEx surface.
- calibration_qc.csv: per-irradiation calibration and consistency metrics.
- reduction_report.json: assumptions, completeness, and warnings.
- dose_prim_vs_sec.png and dose_prim_vs_sec_1s.csv: independent primary and
  secondary dose-monitor comparison.
- aggregates/RadEx_64MeV, aggregates/Radex_85MeV, and aggregates/RadEx_Total:
  summaries and time series for the RadEx irradiations.

Recorder-plane heatmaps are retained as diagnostics. They are not DUT-plane
fluence maps and are not used for the RadEx surface-fluence estimate.

## RadEx surface-fluence estimate

Only the RadEx irradiations are considered. The facility was explicitly asked
to irradiate the complete 154.0 mm by 66.8 mm instrument surface and assured
fluence uniformity within 5% across that surface. Its area is 102.872 cm2.

For every irradiation and time bin, mean incident surface fluence is calculated
as monitor-derived delivered protons divided by 102.872 cm2. The exp_12/exp_13
field-size fit and all inferred X/Y magnification factors have been removed.
Logged pencil-beam widths are not used in this estimate.

This is a mean incident fluence over the RadEx surface. The facility assurance
implies local values should be within approximately +/-5% of that mean, subject
to the precise facility definition of uniformity. Transport, scattering, and
energy deposition between the monitor reference and the dosimeters are handled
separately by the Geant4 model.

The diagnostic spatial heatmaps retain a configurable Gaussian spot-width
interpretation, but they do not establish DUT-plane coordinates or fluence.

## Calibration assumptions and QC

- `DOSE_PRIM(C)` is an incremental charge sample and is summed directly; the
  250 microsecond acquisition period is not applied again.
- `DOSE_SEC(C)` is retained as an independent monitor cross-check. Across all
  13 documented irradiations, its integrated charge is 1.6181% higher than
  `DOSE_PRIM`. For active one-second bins the median difference is +1.6060%,
  with a 5th--95th percentile range of +1.4308% to +1.7966%. The stable,
  positive offset is consistent with a relative monitor calibration/gain
  difference; fluence calibration remains based on the facility-specified
  primary monitor.
- `K_FACTOR` is applied per record, so changes such as the exp_13 setting are
  preserved.
- `part` records contribute to delivered MU and fluence. `tuning` records are
  reported separately and are negligible for these documented irradiations.
- All 367 discovered `map_record` files were parsed (29,241,769 retained acquisition
  rows, zero numerical parse errors). The single orphan specification is an
  undocumented tuning file from 27 June and has no associated record.
