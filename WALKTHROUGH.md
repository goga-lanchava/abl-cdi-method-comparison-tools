# Reproducing the article's results

This walkthrough reproduces the representative analyses and figures shown
in the accompanying SoftwareX article, using the data in [`examples/`](examples/).

## 1. PatLogGUI — import and cleaning (article Fig. 3)

1. Launch `PatLogGUI` (source: `PatLogGUI/src/PatLogGUI.m`, or the compiled
   executable if you're using a release build).
2. **Import Patient Data File** → select `examples/PatLog_export.csv`.
3. **Process Clinical Data**. The status message should report the number
   of essential columns kept (18) and rows retained after removing
   calibration/QC/metadata rows.
4. Select **Patient ID** `2026027`, **Blood Gas Parameter** `sO2`, then
   click **Trend Analysis**. This reproduces the cleaned-data view shown in
   the article's Fig. 3.

## 2. ABG Analyzer Enhanced — method comparison (article Figs. 4–5)

This reproduces the pH result for Dataset 2026027 reported in the article
(uncorrected bias +0.047, SD 0.095; Hybrid-corrected bias 0.000, SD 0.042,
a 56% SD reduction).

1. Launch `ABG_Analyzer_Enhanced` (source or compiled executable).
2. **Load ABL Data** → `examples/PatLog_export.csv`.
3. **Load CDI Data** → `examples/2026027_cdi.log`.
4. Set **Patient ID** to `2026027`, **Select Parameter** to `pH`, **Time
   Tolerance** to `5` minutes.
5. Click **ANALYZE**. The Statistics panel should show approximately:
   `N Pairs: 67`, `Bias: 0.047`, `SD: 0.095`, `95% LoA: [-0.139, 0.236]` —
   this is the three-panel visualization in article Fig. 4.
6. Under **CDI Correction**, select **Auto (Best Model)** and click **Apply
   Correction**. The LOO-CV evaluation should select **Hybrid
   (Time-Series + Deming)** as the winner (τ ≈ 2.0 min, β0 ≈ −5.59,
   β1 ≈ 1.78), with the corrected statistics showing bias ≈ 0.000 and SD ≈
   0.042 (a ~56% SD reduction) — matching the numbers reported in the
   article's Illustrative examples section.
7. Click **Export Figures (SVG)** to reproduce the composite figure shown
   in article Fig. 5.

## 3. Dataset 2026007 (the small/unstable recording)

Repeat step 2 with `examples/2026007_cdi.log` and Patient ID `2026007` to
reproduce the N=10 result discussed in the article's Impact section (where
several correction methods were within 1% LOO-CV RMSE of one another).

## Notes

- Exact statistics may vary by a small amount (last-digit rounding)
  depending on MATLAB version and floating-point behavior, but should
  match the values above to within normal numerical tolerance.
- See [`ABG_Analyzer_Enhanced/docs/DATA_FORMATS.md`](ABG_Analyzer_Enhanced/docs/DATA_FORMATS.md)
  and [`PatLogGUI/docs/DATA_FORMAT.md`](PatLogGUI/docs/DATA_FORMAT.md) for
  the exact input file format specifications.
