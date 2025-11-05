# AI Agent Instructions for TritonScripts Repository

This codebase contains Python scripts for analyzing particle radiation effects, shielding calculations, and radiation dose estimation using GRAS (Geant4 Radiation Analysis for Space).

## Project Structure

- **Carrington/**: Analysis of Carrington event radiation effects
  - Electron shielding curves, spectra analysis, and time dependence studies
  - SEE (Single Event Effects) rate calculations
  
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
  - Uses matplotlib with project-specific formatting

## Key Design Patterns

1. **Data Processing Pipeline**:
   ```python
   from Dependencies.TotalDose import totalkRadGras
   # Read and process radiation data
   data = totalkRadGras(path, "Elec")  # Common pattern for electron data
   data = totalkRadGras(path, "Prot")  # Common pattern for proton data
   ```

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
   from Dependencies.TotalDose import totalkRadGras
   from Dependencies.TotalDoseHistos import totalGRASHistos
   from Dependencies.TotalLETHistos import totalLETHistos
   ```

2. **Error Propagation**:
   ```python
   # Standard error combination for total dose
   total[1] = np.sqrt(electrons[1]**2 + protons[1]**2)
   ```

3. **Unit Conversions**:
   - Time: Usually in months (30.44 days) or hours
   - Dose: Usually in kRad/month
   - Remember to scale between flux and fluence appropriately

## Integration Points

1. **GRAS Dependencies**:
   - Core package: `GRAS.Dependencies`
   - Key modules: TotalDose, MergeHistograms, TotalLETHistos
   - Import hierarchy matters - use absolute imports

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