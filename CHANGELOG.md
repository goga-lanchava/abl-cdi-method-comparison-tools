# Changelog

## Unreleased

Changes on `main` since v1.1.0. None of them alters any statistic quoted in
`WALKTHROUGH.md`.

### Fixed

- The plotted and exported corrected series for an Auto-selected Hybrid model
  is now computed with the same derivative-clip threshold as the deployed
  coefficients. The fit used the threshold of the fitting-window span, while
  the full corrected series shown in the time-series panel, the composite
  figure and **Export Corrected** was recomputed with the whole-recording
  threshold. The deployed model now records the mask it was fitted with and the
  display path reuses it, so the series drawn is exactly the one the reported
  bias and SD describe. Manual Hybrid fits keep the whole-recording threshold
  they are fitted with, and a manual fit no longer inherits a mask left over
  from a previous Auto run.

  No reported statistic changes: the two thresholds produce identical values
  at every paired sample, so bias, SD, limits of agreement, correlation and
  LOO-CV RMSE are all unchanged. Only the drawn corrected trace of the pO2
  example shifts slightly between the paired samples (419 of 5,814 samples,
  by at most 8.4 mmHg and 0.45 mmHg on average).

### Documentation

- Descriptions corrected to match the software: τ_rise and τ_fall share one
  tuned value in both automatic paths; AutoShift estimates a temporal offset
  from the data rather than correcting a known clock offset; Auto model
  selection is single-level LOO-CV with in-fold tuning, not nested
  cross-validation; PatLogGUI detects the CSV delimiter itself; PatLogGUI
  cleaning selects columns and removes no rows.
- `examples/README.md`: the reproduction steps now include the required
  **Fit Window → Auto** step for Dataset 2026027, and the description of the
  example export states what its fields actually contain.
- `WALKTHROUGH.md`: added the limits of agreement on the retained pairs,
  exact values where approximate ones were given, and an explicit tolerance;
  corrected cross-references and the description of the exported
  Bland–Altman panels.
- `TESTING.md`: a step-by-step guide to running the test suite, including
  what the output looks like, how to read a failure, and a table of what each
  of the 23 tests checks.
- `WALKTHROUGH.md` now states the expected correlation values before and
  after correction, which the figures show but nothing previously let a reader
  check. They are asserted in `WalkthroughTest`.
- README: a fuller Tests section pointing to `TESTING.md`.

### Tests

- `WalkthroughTest` now asserts every value quoted in `WALKTHROUGH.md`,
  including the retained-pair statistics and limits of agreement, the LOO-CV
  RMSE values and LoA spans behind the tie-breaker, the values obtained
  without the fitting window, and PatLogGUI's column and row counts. The
  number of tests is unchanged (23).

## v1.1.0 — 2026-09-19

### Changed - leave-one-out cross-validation is now fully fold-independent

Hybrid's hyperparameters (W1, τ, λ) were already tuned inside each training
fold, and the MAD mask was already recomputed from the training pairs only.
Three further quantities were not:

- the 95th-percentile derivative clip was computed once over the **whole** CDI
  recording and baked into the filtered traces every fold was scored on;
- the raw roughness reference and the per-candidate roughness used in the
  Hybrid tuning penalty were computed once over the **whole** fitting window;
- the window those roughness statistics were taken over was bounded by the
  minimum and maximum of **all** paired times, including the held-out pair.

All three are now derived from the training pairs of the fold that uses them.
`computeAsymmetricFastCDI` takes an optional mask restricting which samples
define the clip threshold, and the Auto path precomputes one set of traces and
roughness statistics per distinct training span. Dropping an interior pair does
not move the span, so leave-one-out needs only three spans (interior,
earliest-dropped, latest-dropped) and the cost stays close to the previous
single precompute.

The held-out observation no longer influences parameter selection in any way.
Note that the clip threshold remains a within-span statistic rather than a
running causal one: it never uses the held-out pair, but it is not computed
strictly from samples preceding each time point.

Effect on the documented examples:

- **pH / 2026027** - selected model, deployed coefficients and all reported
  before/after statistics are unchanged. Hybrid's LOO-CV RMSE rose slightly
  (0.095654 to 0.095699), so the gap to the winner widened from 0.37% to
  0.42%; it stays inside the 1% band, the limits-of-agreement tie-breaker
  still engages, and Weighted Deming still wins.
- **pO2 / 2026007** - Hybrid still wins but retunes: RMSE 218.1773 to
  218.3203, τ_rise = τ_fall 7.0 to 8.0 min, λ 0.25 to 0.10, deployed fit
  `(CDI_fast - 2191.1501) / -7.1327`. After-correction bias 0.6979 to 0.1180,
  SD 161.9964 to 161.6657, SD reduction 41.7% to 41.8%, 95% LoA now
  [-316.7, 317.0]. W1 = 4 and 10/10 MAD-retained pairs are unchanged.

RMSE rises slightly in both cases, as expected: the previous values were
computed with tuning quantities that included the held-out observation, which
made them slightly lower.

### Changed - Deming λ is documented consistently as σ²(CDI)/σ²(ABL)

`fitWeightedDeming` implements the standard Deming form with x = ABL and
y = CDI, in which λ is the ratio of the y-error to the x-error variance - that
is, σ²(CDI)/σ²(ABL). The GUI label read
"Variance Ratio ABL/CDI", the inverse. The label, the function comment and the
Correction Report wording now all state λ = σ²(CDI)/σ²(ABL). λ = 1 remains
orthogonal regression. No numerical behaviour changed.

### Changed - correction quality verdict is purely descriptive

The verdict applied a fixed +/-5% band to the SD reduction and reported
"IMPROVED" / "WORSENED" / "NO IMPROVEMENT" with pass-fail marks. It now
compares the before and after values directly on the same MAD-retained pairs
and reports one of `BIAS + SD REDUCED`, `BIAS REDUCED`, `SD REDUCED` or
`NO REDUCTION`, in a neutral colour. No threshold is applied, since the
parameters differ in unit and scale, and the label makes no claim about
clinical acceptability. The SD and LoA percentages are still shown alongside.

### Added

- `tests/ComponentTest.m`: ten deterministic component tests covering ABL and
  CDI parsing of the known example files, rejection of malformed timestamps,
  tolerance of a missing patient-ID column, skipping of malformed CDI lines,
  temporal pairing at two tolerances, MAD filtering on artificial data
  (including the zero-MAD degenerate case), and a Bias Correction whose offset
  and corrected series are checked against an independently computed value.
- `WalkthroughTest/looCvTuningIsFoldIndependent`, which also fails if any
  candidate scores NaN.
- Public `readABL`, `readCDI` and `madRetainMask` wrappers so the parsers and
  the robust filter can be driven headlessly.
- The deployed Hybrid smoothing window (`w1`) and the candidate
  limits-of-agreement spans (`autoLoASpan`) are now recorded on the
  correction model, so the tie-breaker can be inspected and tested.

---

Documentation, fixes and tests, checked against the shipped code and example
data.

### Documentation
- `WALKTHROUGH.md` §2 includes the **Fit Window → Auto** step, on which
  every statistic quoted in §2 depends; without it the analysis uses 68
  pairs instead of 67 and gives Bias 0.049, SD 0.096 and a different fitted
  model. Steps renumbered accordingly.
- `WALKTHROUGH.md`: the observation that the top two correction methods lie
  within 1% LOO-CV RMSE of one another is stated for the pH dataset in §2
  (0.42% apart; the tie-breaker engages). §3 gives the pO2 ranking (34%
  apart), refers to step 3, and notes that no fitting window is needed.
- `WALKTHROUGH.md` §1 and PatLogGUI README: documentation no longer describes
  a row-removal or row-validation stage; cleaning selects columns only.
- Toolbox requirements are stated per tool in the root README,
  `CODE_METADATA.md` and the per-tool READMEs: `PatLogGUI` requires the
  Statistics and Machine Learning Toolbox (`boxplot`, in the Statistical
  Analysis module); `ABL_CDI_Analyzer` requires core MATLAB only. A note
  explains that on releases older than R2025b the Analyzer may also need
  the toolbox, because `prctile` moved into core MATLAB only recently.
- The root README and `CODE_METADATA.md` C5 describe `PatLogGUI` as a
  programmatic `uifigure` GUI rather than an App Designer app.

### Fixed
- **Export Results** in ABL-CDI Analyzer now writes the `.xlsx` and `.csv`
  files correctly (previously failed on mixed scalar/vector statistics). The
  summary table holds the scalar statistics, with the paired `xABL`/`yCDI`
  values written alongside; CSV exports are written as
  `<name>_statistics.csv` and `<name>_paireddata.csv`.
- `PatLogGUI` detects the file encoding from the file's bytes, so Latin-1
  ABL Flex 800 exports display the temperature column as `T (°C)`.
- Example CDI logs are no longer excluded by `.gitignore`
  (`!examples/*.log`).

### Added
- `tests/WalkthroughTest.m`: regression tests asserting the exact figures
  quoted in `WALKTHROUGH.md` §2 and §3, Export Results, the import
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
  smoothing time constant was split into separate rise/fall parameters
  (τ_rise, τ_fall), which can be set independently in manual mode; AutoTune
  and Auto (LOO-CV) tune a single shared value (τ_rise = τ_fall)
- Auto (Leave-One-Out Cross-Validation) model selection across 7 candidate
  correction methods
- TimeShift / AutoShift: manual or correlation-estimated temporal offset
  between the ABL and CDI streams
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
  and/or the tie-breaker outcome relative to earlier versions.
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
  (0.047 → 0.043, 61/67 pairs), including where in the software's output
  each figure is found

## PatLogGUI v1
- Automatic delimiter (semicolon/comma) and decimal-comma handling for CSV
  imports; Excel files read directly
- Cleaning by essential-column selection
- Three-dimensional filtering (patient, catheter location, time point)
- Trend Analysis, Multi-Parameter View, Location Comparison visualizations
- Statistical Analysis module with descriptive statistics and boxplots
- Cleaned/patient-filtered data export to Excel
