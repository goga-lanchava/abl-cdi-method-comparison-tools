# Changelog

## ABL-CDI Analyzer (renamed from ABG Analyzer Enhanced)
- Renamed application/class from `ABG_Analyzer_Enhanced` to
  `ABL_CDI_Analyzer`; window title simplified to "ABL-CDI Analyzer"
- LOO-CV model selection: eliminated a data-leakage issue for both the
  Deming error-variance ratio (λ) and the Hybrid method's Auto-Tune
  hyperparameters (τ, smoothing window, λ). These are now re-optimized by
  grid search strictly on each fold's N-1 training partition, rather than
  reusing the value optimized on the full dataset. This can shift the exact
  LOO-CV RMSE values, the winning candidate, and/or the tie-breaker outcome
  relative to earlier results — re-run `WALKTHROUGH.md` against this version
  before citing specific numbers.
- Corrected-series output: all six correction methods now consistently set
  MAD-filtered outlier points to `NaN` in the exported/plotted corrected
  series (previously only Passing–Bablok did this); before/after statistics
  computed from the corrected series may shift slightly as a result
- Before/after "Correction quality" verdict relabelled from clinical-sounding
  terms (GOOD/MODERATE/LIMITED) to neutral, software-defined descriptive
  labels (e.g. "Bias + SD improved"), with an explicit note that these are
  not clinical acceptance criteria
- Fixed a LaTeX export bug where the τ (tau) symbol in the Hybrid model
  formula used a literal Unicode character instead of the `\tau` command

## ABG Analyzer Enhanced v9.1 (State-Locked)
- Hybrid correction method: Time-Series + Deming with AutoTune grid search,
  MAD-based outlier filtering, and 95th-percentile derivative cap
- Auto (Leave-One-Out Cross-Validation) model selection across 6 candidate
  correction methods
- TimeShift / AutoShift for instrument clock offset alignment
- CorrView / BAView toggles; 4-panel Before/After comparison popup
- Pin correction; LaTeX export of the fitted correction formula
- Small-N warning label; continuous CDI line plot

## PatLogGUI v1.05
- Automatic delimiter/decimal-notation detection for CSV/Excel imports
- Two-stage cleaning (essential-column selection, patient-ID row validation)
- Three-dimensional filtering (patient, catheter location, time point)
- Trend Analysis, Multi-Parameter View, Location Comparison visualizations
- Statistical Analysis module with descriptive statistics and boxplots
- Cleaned/patient-filtered data export to Excel


