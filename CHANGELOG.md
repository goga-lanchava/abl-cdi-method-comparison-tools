# Changelog

## Unreleased

Fixes from an independent reproduction check of `WALKTHROUGH.md` against the
shipped code and example data.

### Documentation
- `WALKTHROUGH.md` §2: added the missing **Fit Window → Auto** step. Every
  statistic quoted in §2 depends on it, and following the previous steps
  literally produced none of them (68 pairs instead of 67, Bias 0.049,
  SD 0.096, and a different fitted model). Steps renumbered accordingly.
- `WALKTHROUGH.md` §3: the "several correction methods within 1% LOO-CV
  RMSE of one another" observation applies to the pH dataset in §2 (top two
  candidates 0.37% apart, tie-breaker engages), not to the pO2 dataset in
  §3 (34% apart). Moved to §2 and replaced with the actual §3 ranking.
  Also corrected a cross-reference ("repeat step 2" → step 3) and noted
  that §3 needs no fitting window.
- `WALKTHROUGH.md` §1 and PatLogGUI README: cleaning selects columns only.
  The previously documented removal of calibration/QC/metadata rows and the
  "row validation (patient-ID presence, plausible temperature range)" stage
  were never implemented; the description now matches the behaviour.
- Toolbox requirements corrected: `PatLogGUI` requires the Statistics and
  Machine Learning Toolbox (`boxplot`, in the Statistical Analysis module).
  `ABL_CDI_Analyzer` is core MATLAB only. Both were previously documented
  as core-MATLAB-only in the root README, `CODE_METADATA.md` and the
  per-tool READMEs. Added a note that on releases older than R2025b the
  Analyzer may also need it, because `prctile` moved into core only
  recently.
- `PatLogGUI` is a programmatic `uifigure` GUI, not an App Designer app;
  the root README and `CODE_METADATA.md` C5 said otherwise.

### Fixed
- **Export Results** in ABL-CDI Analyzer never produced a file, in any
  format: `struct2table(app.Stats)` always threw, because `Stats` mixes
  scalar statistics with the N-element `xABL`/`yCDI` vectors, and the error
  was swallowed into an "Export failed" alert. The summary table now holds
  the scalar fields only, with the paired vectors written alongside. The
  `.csv` option also passed `'Sheet'` to `writetable`, which CSV does not
  support; CSV exports are now written as `<name>_statistics.csv` and
  `<name>_paireddata.csv`.
- `PatLogGUI` forced UTF-8 when reading imports, but ABL Flex 800 exports
  are typically Latin-1, so the temperature column displayed as `T (?C)`.
  The encoding is now detected from the file's bytes.
- `.gitignore` matched `*.log`, silently ignoring CDI monitor recordings —
  one of the project's two primary input formats. Example logs are now
  exempt via `!examples/*.log`.

### Added
- `tests/WalkthroughTest.m`: regression tests asserting the exact figures
  quoted in `WALKTHROUGH.md` §2 and §3, the export fix, the import
  encoding, and the declared toolbox dependencies. Run with
  `runtests('tests')`.
- `ABL_CDI_Analyzer.runWorkflow(...)`: scripted equivalent of the
  walkthrough workflow, returning the statistics exactly as the panel shows
  them, for batch runs and testing.
- `ABL_CDI_Analyzer.exportResultsTo(path)`: the export writer, split out of
  the button callback so it can be called without a file dialog.
- `PatLogGUI(filePath)`: opens with that export already imported, skipping
  the file dialog. `PatLogGUI()` behaves as before.

## ABL-CDI Analyzer v1 (renamed from ABG Analyzer Enhanced)
- Added Weighted Deming (Linnet-style, error-variance ratio λ tuned by grid
  search) as a 7th correction method and Auto/LOO-CV candidate
- Hybrid correction method: Time-Series + Deming with AutoTune grid search,
  MAD-based outlier filtering, and 95th-percentile derivative cap; the
  smoothing time constant was split into independently Auto-Tuned rise/fall
  parameters (τ_rise, τ_fall) rather than a single symmetric τ
- Auto (Leave-One-Out Cross-Validation) model selection across 7 candidate
  correction methods
- TimeShift / AutoShift for instrument clock offset alignment
- CorrView / BAView toggles; 4-panel Before/After comparison popup
- Pin correction; LaTeX export of the fitted correction formula
- Small-N warning label; continuous CDI line plot
- Renamed application/class from `ABG_Analyzer_Enhanced` to
  `ABL_CDI_Analyzer`; window title simplified to "ABL-CDI Analyzer"
- LOO-CV model selection: eliminated a data-leakage issue for both the
  Deming error-variance ratio (λ) and the Hybrid method's Auto-Tune
  hyperparameters (τ_rise, τ_fall, smoothing window, λ). These are now
  re-optimized by grid search strictly on each fold's N-1 training
  partition, rather than reusing the value optimized on the full dataset.
  This can shift the exact LOO-CV RMSE values, the winning candidate,
  and/or the tie-breaker outcome relative to earlier results — re-run
  `WALKTHROUGH.md` against this version before citing specific numbers.
- Corrected-series output: all seven correction methods now consistently
  set MAD-filtered outlier points to `NaN` in the exported/plotted
  corrected series (previously only Passing–Bablok did this); before/after
  statistics computed from the corrected series may shift slightly as a
  result
- Before/after "Correction quality" verdict relabelled from clinical-sounding
  terms (GOOD/MODERATE/LIMITED) to neutral, software-defined descriptive
  labels (e.g. "Bias + SD improved"), with an explicit note that these are
  not clinical acceptance criteria
- Fixed a LaTeX export bug where the τ (tau) symbol in the Hybrid model
  formula used a literal Unicode character instead of the `\tau` command
- `WALKTHROUGH.md` §2 updated to distinguish the all-pairs pre-correction
  SD (0.095, 67/68 pairs) from the MAD-retained-pairs-only SD comparison
  (0.047 → 0.043, 61/67 pairs) the article now quotes, including where in
  the software's output each figure is found

## PatLogGUI v1
- Automatic delimiter/decimal-notation detection for CSV/Excel imports
- Two-stage cleaning (essential-column selection, patient-ID row validation)
- Three-dimensional filtering (patient, catheter location, time point)
- Trend Analysis, Multi-Parameter View, Location Comparison visualizations
- Statistical Analysis module with descriptive statistics and boxplots
- Cleaned/patient-filtered data export to Excel
