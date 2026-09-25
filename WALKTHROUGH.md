# Reproducing the article's results

This walkthrough reproduces the representative analyses and figures shown
in the accompanying SoftwareX article, using the data in [`examples/`](examples/).

## 1. PatLogGUI — import and cleaning (article Fig. 3)

1. Launch `PatLogGUI` (source: `PatLogGUI/src/PatLogGUI.m`, or the compiled
   executable if you're using a release build).
2. **Import Patient Data File** → select `examples/PatLog_export.csv`.
3. **Process Clinical Data**. The panel should report `Kept 18 essential
   columns` and the number of non-essential columns removed (218 for this
   file). Cleaning selects columns only — no rows are removed at this
   stage. The status line below then reports the number of rows on display
   (1907 with **All Patients** selected).
4. Select **Patient ID** `2026027`, **Blood Gas Parameter** `sO2`, then
   click **Trend Analysis**. This reproduces the cleaned-data view shown in
   the article's Fig. 3.

## 2. ABL-CDI Analyzer — method comparison (article Fig. 4)

This reproduces the pH result for Dataset 2026027 reported in the article.
The article quotes two different "before" baselines, and this walkthrough
tells you where each one comes from so both are independently checkable:

- The **headline, all-pairs** pre-correction agreement (67 window pairs,
  68 total): bias +0.047, SD 0.095, 95% LoA [-0.140, 0.234]. This is what
  the live Statistics panel shows immediately after **ANALYZE**, before
  any correction is applied.
- The article's **"-10% SD reduction"** compares like-for-like: the SD on
  only the 61 of 67 pairs that survive MAD-based filtering, before (SD
  0.0475) and after (SD 0.0426) correction. That is a 10.3% reduction on
  unrounded values, which the software displays as `SD ▼10.3%` and the
  article quotes as -10% (recomputing from the rounded 3-decimal figures
  0.047 → 0.043 instead gives ≈8.5%). On the same 61 pairs the 95% limits
  of agreement narrow from [-0.020, 0.166] to [-0.084, 0.083]. These
  retained-pairs "before" figures are *not* shown on the live Statistics
  panel (which always reports the all-pairs figures); they appear in the
  "Bland-Altman Before" panel of the exported composite report, because
  `ExportFigureButtonPushed` recomputes the before-stats on the same pairs
  used for the after-stats (see step 8).


1. Launch `ABL_CDI_Analyzer` (source or compiled executable).
2. **Load ABL Data** → `examples/PatLog_export.csv`.
3. **Load CDI Data** → `examples/2026027_cdi.log`.
4. Set **Patient ID** to `2026027`, **Select Parameter** to `pH`, **Time
   Tolerance** to `5` minutes.
5. Tick the **Fit Window** checkbox, then click its **Auto** button. This
   selects the stable fitting window (`14.04.2026 16:19` to
   `14.04.2026 21:03`) and restricts the analysis to the 67 window pairs
   that every figure in this section refers to.

   > **This step is required.** Without it the analysis runs on all 68
   > pairs and none of the values below reproduce: you get Bias 0.049,
   > SD 0.096, 95% LoA [-0.139, 0.236], and a different fitted model
   > (`Raw_CDI = 1.8666*ABL -6.2502`, 68/68 pairs, 62 robust, SD down
   > 12.8%). The `(window)` wording in the N Pairs label below also only
   > appears once the Fit Window checkbox is ticked.
6. Click **ANALYZE**. The Statistics panel should show approximately:
   `N Pairs: 67 (window) / 68 total`, `Bias: 0.047`, `SD: 0.095`,
   `95% LoA: [-0.140, 0.234]`, `r = 0.1700` — this is the "before
   correction" panel
   shown in article Fig. 4 (top), and the all-pairs figures quoted in
   Section 3 of the article.
7. Under **CDI Correction**, select **Auto (Best Model)** and click **Apply
   Correction**. The LOO-CV evaluation (across all 7 candidates) should
   select **Weighted Deming (Linnet, tuned λ)** as the winner. This is the
   dataset where the top two candidates fall within the 1% RMSE band that
   triggers the limits-of-agreement tie-breaker: Weighted Deming scores
   RMSE 0.0953 against the Hybrid method's 0.0957 — a gap of 0.42% — so
   the narrower LoA span decides (0.3628 vs 0.3646), confirming Weighted
   Deming as the winner. The Correction Report's formula text should
   read:
   ```
   Model: Raw_CDI = 1.7395*ABL -5.3238  [Linnet Weighted Deming λ=0.10]
   Corrected = (Raw_CDI +5.3238) / 1.7395
   (Fitted on 67/68 window pairs, 61 robust [|diff - med| <= 4.5*MAD])
   ```
   confirming λ = 0.10, coefficients 1.7395 / −5.3238, and the 61-of-67
   retained-pair count quoted in the article. The corrected statistics
   should show `After Correction: Bias=-0.0002`, `SD=0.0426 (on 61 kept
   pairs)` and the descriptive label `BIAS + SD REDUCED | SD ▼10.3%
   LoA ▼10.3%` — matching the "after correction" panel shown in article
   Fig. 4 (bottom). With **Show corrected scatter** ticked, the correlation
   panel title becomes `After Correction (r=0.471)`, against `r=0.170`
   before correction.
8. Click **Export Figures (SVG)** and open the resulting `*_report.svg`
   (the composite shown as article Fig. 5 for the pO2 dataset — see §3 —
   but generated here for pH too). Its **"Bland-Altman Before"** panel
   reports the pre-correction statistics recomputed on the same 61 retained
   pairs used for the after-stats: title `Bias=0.0729  SD=0.0475`, legend
   `LoA=[-0.020,0.166]`. The **"Bland-Altman After"** panel's legend shows
   `Bias=-0.0002` and `LoA=[-0.084,0.083]`. Together with the
   after-correction SD of 0.0426 from step 7, this is the SD 0.047 → 0.043
   (10.3% reduction on unrounded values) and the narrowed limits of
   agreement quoted in the article's Fig. 4 discussion.

## 3. Dataset 2026007 (article Fig. 5; small recording)

Repeat step 3 with `examples/2026007_cdi.log` and Patient ID `2026007`,
selecting **pO2** as the parameter, to reproduce the N=10 result discussed
in Section 3 and the Conclusions of the article. Leave **Fit Window** unticked for this
dataset: all 10 pairs are used, and no window selection is needed. Unlike
the pH dataset in §2, the LOO-CV candidates here are widely separated —
the Hybrid method wins outright with RMSE 218.3203 against 292.9552 for
the runner-up (Bias Correction), a gap of 34.2% — so the 1% tie-breaker
never engages. Before correction: Bias=48.460,
SD=277.922, 95% LoA=[-496.267, 593.187], r=-0.1088. Auto (Best Model)
selects **Hybrid (Time-Series + Deming)** (RMSE=218.3203), with winning
parameters W1=4, τ_rise=τ_fall=8.0 min, Linnet λ=0.10, deployed as
`CDI_corrected = (CDI_fast - 2191.1501) / -7.1327`. Corrected statistics:
Bias=0.1180, SD=161.6657, and a correlation against the ABL draws that
turns from r=-0.109 before correction to r=0.201 after it (the
`Correlation Before` and `Correlation After` panels of the composite
export), giving the software's own displayed "SD ▼41.8%,
LoA ▼41.8%" reduction (reported as ~42% in the article, on unrounded
values) and narrowing the 95% limits of agreement to [-316.747, 316.983]
(legend of the exported "Bland-Altman After" panel). The Correction Report's formula text should confirm all
10 of 10 pairs are retained after MAD-based filtering (no exclusions) —
so, unlike the pH example in §2, the before/after SD comparison for pO2
uses the same N=10 pairs throughout without needing the retained-pairs
distinction. Click **Export Figures (SVG)** to reproduce the composite
figure shown in article Fig. 5.

> **Note:** the Patient ID dropdown may show as "All Patients" rather than
> `2026007` after loading `2026007_cdi.log` — because this CDI recording
> only overlaps in time with ABL reference draws from patient 2026007, the
> resulting pairing and statistics are identical either way. If you see a
> different N or different statistics with "All Patients" selected, ABL
> draws from another patient are being pulled into the window; explicitly
> select `2026007` in that case.

## Notes

- **Tolerance.** All values above are given to the decimals the software
  displays, and reproduce exactly on MATLAB R2025b. A difference of one
  unit in the last displayed digit can arise from floating-point
  differences between machines and is not a discrepancy; anything larger
  is. The automated tests ([TESTING.md](TESTING.md)) compare displayed
  values to the decimals shown here, and unrounded model coefficients and
  correlations to within an absolute tolerance of 5×10⁻⁴.
- Percentage reductions quoted in the article are computed from unrounded
  intermediate values, not from the 3-4 decimal figures printed in the UI
  or reported here; recomputing from the rounded figures alone can give a
  visibly different percentage (e.g. ≈8.5% vs. the displayed 10.3% for
  the pH example in §2) without indicating an error.
- See [`ABL_CDI_Analyzer/docs/DATA_FORMATS.md`](ABL_CDI_Analyzer/docs/DATA_FORMATS.md)
  and [`PatLogGUI/docs/DATA_FORMAT.md`](PatLogGUI/docs/DATA_FORMAT.md) for
  the exact input file format specifications.
