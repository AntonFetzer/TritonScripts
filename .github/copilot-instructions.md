# AI Agent Instructions for TritonScripts Repository

This codebase contains Python scripts for analyzing particle radiation effects, shielding calculations, and radiation dose estimation based on simulation outputs from GRAS (Geant4 Radiation Analysis for Space).

## Project Structure

- **Carrington/**: Analysis of Carrington event radiation effects
  - Electron shielding curves, spectra analysis, and time dependence studies
  - SEE (Single Event Effects) rate calculations specific to the Carrington project
  
- **Dependencies/**: Core utilities for data processing
  - `TotalDose.py`, `MergeHistograms.py` - Essential histogram processing
  - Used across the codebase for radiation data analysis
  
- **SEE/**: Single Event Effects analysis
  - Cross section calculations and rate estimation
  - Weibull function fitting for SEE modeling

- **Read/**: Data input modules
  - Parsers for various file formats (Dose, LET, SPENVIS, etc.)
  - Example: `ReadDose.py` for reading dose measurement files

- **Plotting/**: Visualization modules
  - Supports various plot types: shielding curves, dose histograms, etc.
  - Use matplotlib defaults where possible

## Key Design Patterns

1. **Data Processing Pipeline**:
   ```python
   from Dependencies.TotalDose import totalDose

   # Point each call at a directory containing one GRAS result dataset.
   electron = totalDose(electron_results_path)
   proton = totalDose(proton_results_path)

   electron_dose = electron["dose"]
   electron_error = electron["error"]
   ```

   `totalDose` returns a dictionary with `dose`, `error`, `entries`, and
   `non-zeros` arrays. The retired function is not available. For legacy
   folders containing several datasets, use the optional `filename_contains`
   argument explicitly, for example
   `totalDose(path, filename_contains="Elec")`.

2. **Error Handling**:
   - Always check histogram consistency:
   ```python
   required_keys = ['lower', 'upper', 'mean', 'value', 'error', 'entries']
   lengths = [len(histogram[k]) for k in required_keys]
   if not all(l == lengths[0] for l in lengths):
       print("Inconsistent histogram lengths")
       continue
   ```

3. **File Path Conventions**:
   - Results typically stored in "Res/" subdirectories
   - Use absolute paths from `/l/triton_work/`
   - Common path structure: `{study_type}/{configuration}/{particle_type}/Res/`

## Common Operations

1. **Reading Data**:
   ```python
   from Dependencies.TotalDose import totalDose
   from Dependencies.TotalDoseHistos import totalGRASHistos
   from Dependencies.TotalLETHistos import totalLETHistos
   ```

2. **Error Propagation**:
   ```python
   # Standard error combination for total dose
   total_error = np.sqrt(electrons["error"]**2 + protons["error"]**2)
   ```

3. **Unit Conversions**:
   - Time: Usually in months (30.44 days) or hours
   - Dose: Usually in kRad/month
   - Remember to scale between flux and fluence appropriately

## Integration Points

1. **GRAS Dependencies**:
   - Repository package: `Dependencies`
   - Key modules: TotalDose, MergeHistograms, TotalLETHistos
   - Run from the repository root so absolute `Dependencies.*` imports resolve

2. **Data Formats**:
   - CSV files with specific column structures
   - Histogram data with required keys: lower, upper, mean, value, error, entries
   - SPENVIS compatibility for space radiation environment data

## Debugging Tips

1. **Common Issues**:
   - Check histogram consistency before merging
   - Verify unit conversions and scaling factors
   - Ensure proper error propagation in calculations

2. **Data Validation**:
   - Print shape of data arrays to verify dimensions
   - Check for NaN or infinity values in calculations
   - Verify histogram binning matches expectations
