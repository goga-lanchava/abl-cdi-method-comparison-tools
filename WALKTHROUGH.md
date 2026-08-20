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

## 2. ABL-CDI Analyzer — method comparison (article Fig. 4)

This reproduces the pH result for Dataset 2026027 reported in the article
(uncorrected bias +0.047, SD 0.095; after Auto correction, bias ≈ 0.000,
SD ≈ 0.043 — roughly a 55% SD reduction, narrowing the 95% limits of
agreement to approximately [-0.084, 0.083]).

> **Note:** the current source includes a fix to the LOO-CV procedure
> (removing a hyperparameter data-leakage issue for the Deming λ and Hybrid
> Auto-Tune parameters), adds Weighted Deming (Linnet) as a 7th Auto/LOO-CV
> candidate, and changes how outlier points are masked in the exported
> corrected series (see CHANGELOG.md). These can shift the exact figures
> below, and even which candidate wins, by a small amount; re-verify against
> your build rather than treating these as exact expected output.

1. Launch `ABL_CDI_Analyzer` (source or compiled executable).
2. **Load ABL Data** → `examples/PatLog_export.csv`.
3. **Load CDI Data** → `examples/2026027_cdi.log`.
4. Set **Patient ID** to `2026027`, **Select Parameter** to `pH`, **Time
   Tolerance** to `5` minutes.
5. Click **ANALYZE**. The Statistics panel should show approximately:
   `N Pairs: 67`, `Bias: 0.047`, `SD: 0.095`, `95% LoA: [-0.139, 0.236]` —
   this is the "before correction" panel shown in article Fig. 4 (top).
6. Under **CDI Correction**, select **Auto (Best Model)** and click **Apply
   Correction**. The LOO-CV evaluation (now run across all 7 candidates)
   should select **Weighted Deming (Linnet, tuned λ)** as the winner
   (λ ≈ 0.10; Raw_CDI ≈ 1.7395·ABL − 5.3238), with the corrected statistics
   showing Bias ≈ -0.0002 and SD ≈ 0.0426 (61 kept pairs after MAD
   filtering) — matching the "after correction" panel shown in article
   Fig. 4 (bottom).

## 3. Dataset 2026007 (article Fig. 5; the small/unstable recording)

Repeat step 2 with `examples/2026007_cdi.log` and Patient ID `2026007`,
selecting **pO2** as the parameter, to reproduce the N=10 result discussed
in the article's Impact section (where several correction methods were
within 1% LOO-CV RMSE of one another). Here, Auto (Best Model) selects
**Hybrid (Time-Series + Deming)**, reducing bias from ≈48.5 to ≈0.7 and
narrowing the 95% limits of agreement from roughly [-496, 593] to
[-317, 318] (~42% SD reduction). Click **Export Figures (SVG)** to
reproduce the composite figure shown in article Fig. 5.

## Notes

- Exact statistics may vary by a small amount (last-digit rounding)
  depending on MATLAB version and floating-point behavior, but should
  match the values above to within normal numerical tolerance.
- See [`ABL_CDI_Analyzer/docs/DATA_FORMATS.md`](ABL_CDI_Analyzer/docs/DATA_FORMATS.md)
  and [`PatLogGUI/docs/DATA_FORMAT.md`](PatLogGUI/docs/DATA_FORMAT.md) for
  the exact input file format specifications.
