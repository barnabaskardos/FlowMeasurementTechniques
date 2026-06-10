# Flow Measurement Techniques Lab Processing

This branch contains a reproducible processing pipeline for the AE4-180 Flow
Measurement Techniques laboratory exercise.

The workflow is intentionally split by responsibility:

- MATLAB performs calibration, HWA statistics, spectra, PIV data extraction,
  self-made PIV cross-correlation, and DaVis comparisons.
- Python reads the MATLAB-generated CSV tables and creates report-ready plots.

## Quick Start

From the repository root:

```powershell
matlab -batch "run(fullfile('analysis','process_data.m'))"
python analysis\plot_results.py
```

The top-level MATLAB files are also kept as convenient entry points:

```powershell
matlab -batch "PIV"
matlab -batch "HW_Post"
```

Both regenerate the same MATLAB processing tables.

## Dependencies

- MATLAB R2024a or compatible
- MATLAB Signal Processing Toolbox (`xcorr`, `xcorr2`, `pwelch`, `hann`)
- Python 3.10+
- Python packages in `requirements.txt`

Install Python packages with:

```powershell
python -m pip install -r requirements.txt
```

## Outputs

MATLAB writes CSV tables to:

```text
outputs/tables/
```

Python writes figures to:

```text
outputs/figures/
```

Key outputs include:

- `hwa_calibration_curve.*`
- `hwa_velocity_profiles.*`
- `hwa_autocorrelation.*`
- `hwa_energy_spectra.*`
- `piv_mean_fields.*`
- `piv_instantaneous_fields.*`
- `piv_processing_summary.*`
- `piv_davis_window_size_screenshots.*`
- `piv_self_window_fields.*`
- `piv_self_window_sensitivity.*`
- `piv_overlap_multipass_comparison.*`
- `piv_ensemble_size_comparison.*`
- `piv_delta_t_comparison.*`
- `piv_frame_brightness.*`
- `piv_self_vs_davis.*`
- `technique_profile_comparison.*`
- `frequency_resolution.*`

## Requirement Coverage

The pipeline supports the report sections by generating:

- HWA calibration points, PCHIP calibration curve, and 4th order polynomial fit.
- HWA mean and RMS wake profiles for AoA 0, 5, and 15 degrees.
- Autocorrelation, integral time scale, and sampling-time estimate from
  `correlationtest.txt`.
- HWA Welch spectra and dominant frequencies for the requested trailing-edge
  and shear-layer cases.
- Mean and instantaneous PIV fields from the DaVis exports.
- PIV processing sensitivity summaries for overlap, single-pass/multipass, and
  ensemble-size variants available in the processed data folders.
- Side-by-side parameter-study plots showing the effect of interrogation window
  size, overlap, pass strategy, ensemble size, sample spacing, and pulse
  separation on the resulting PIV field.
- Frame-brightness plots for the exported PIV image pair, including frame
  intensity distributions and a lower-minus-upper frame difference map.
- Self-made MATLAB PIV fields for 16, 32, and 64 px interrogation windows.
- Self-made MATLAB PIV versus DaVis comparison for the AoA 15 large-delta-t
  image pair.
- HWA/PIV profile comparison at the configured `x/c = 1.2` station.
- PIV and HWA resolvable frequency information.

## Assumptions To State In The Report

- DaVis `Vx` is negative for the streamwise flow direction in these exports, so
  streamwise velocity is reported as `u = -Vx`.
- The PIV profile comparison uses `x = 20 mm` in the DaVis coordinate system as
  the `x/c = 1.2` wake station. Update `cfg.pivProfileX_mm` in
  `analysis/process_data.m` if your lab notes define a different origin.
- PIV frequency resolution assumes an image-pair sampling rate of `8.3 Hz`,
  based on the camera limit stated in the manual. Update `cfg.pivSampleRate_Hz`
  if the actual acquisition rate was recorded differently.
- The self-made PIV calibration uses the 80 mm ruler segment marked in the
  calibration image. The endpoints are documented in `calibratePivImage`.

## Missing Input

The manual asks for comparison against pressure-probe mean velocity profiles.
No pressure-probe data were found in this repository. The pipeline records this
in:

```text
outputs/tables/missing_required_inputs.csv
```

If pressure-probe data become available, add:

```text
data/pressure_probe_profile.csv
```

with columns:

```text
aoa_deg,y_relative_mm,mean_velocity_ms
```

Then rerun the MATLAB and Python commands above.
