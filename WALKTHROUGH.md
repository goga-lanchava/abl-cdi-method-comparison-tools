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

This reproduces the pH result for Dataset 2026027 reported in the article.
The article quotes two different "before" baselines, and this walkthrough
tells you where each one comes from so both are independently checkable:

- The **headline, all-pairs** pre-correction agreement (67 window pairs,
  68 total): bias +0.047, SD 0.095, 95% LoA [-0.140, 0.234]. This is what
  the live Statistics panel shows immediately after **ANALYZE**, before
  any correction is applied.
- The article's **"-10% SD reduction"** claim compares like-for-like: the
  SD on only the 61 of 67 pairs that survive MAD-based filtering, both
  before (SD ≈ 0.047) and after (SD ≈ 0.043) correction — a ~10% reduction
  on unrounded values (≈8.5% if you recompute from the rounded 3-decimal
  figures; both are reported in the article for transparency). This
  retained-pairs-only "before" SD is *not* shown on the live Statistics
  panel (which always reports the all-pairs figure) — it appears in the
  "Bland-Altman Before" panel title of the exported composite report SVG,
  because `ExportFigureButtonPushed` recomputes the before-stats restricted
  to the same pairs used for the after-stats (see step 7).

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
   `N Pairs: 67 (window) / 68 total`, `Bias: 0.047`, `SD: 0.095`,
   `95% LoA: [-0.140, 0.234]` — this is the "before correction" panel
   shown in article Fig. 4 (top), and the all-pairs figures quoted in the
   article's abstract and Fig. 4 text.
6. Under **CDI Correction**, select **Auto (Best Model)** and click **Apply
   Correction**. The LOO-CV evaluation (across all 7 candidates) should
   select **Weighted Deming (Linnet, tuned λ)** as the winner. The
   Correction Report's formula text should read:
   ```
   Model: Raw_CDI = 1.7395*ABL -5.3238  [Linnet Weighted Deming λ=0.10]
   Corrected = (Raw_CDI +5.3238) / 1.7395
   (Fitted on 67/68 window pairs, 61 robust [|diff - med| <= 4.5*MAD])
   ```
   confirming λ = 0.10, coefficients 1.7395 / −5.3238, and the 61-of-67
   retained-pair count quoted in the article. The corrected statistics
   should show Bias ≈ -0.0002 and SD ≈ 0.0426 (61 kept pairs after MAD
   filtering) — matching the "after correction" panel shown in article
   Fig. 4 (bottom).
7. Click **Export Figures (SVG)** and open the resulting `*_report.svg`
   (the composite shown as article Fig. 5 for the pO2 dataset — see §3 —
   but generated here for pH too). Its **"Bland-Altman Before"** panel
   title reports the pre-correction bias/SD recomputed on the same 61
   retained pairs used for the after-stats: SD ≈ 0.047. Compared against
   the **"Bland-Altman After"** panel's SD ≈ 0.043, this is the SD 0.047 →
   0.043 (≈10% reduction on unrounded values) comparison quoted in the
   article's Fig. 4 discussion.

## 3. Dataset 2026007 (article Fig. 5; the small/unstable recording)

Repeat step 2 with `examples/2026007_cdi.log` and Patient ID `2026007`,
selecting **pO2** as the parameter, to reproduce the N=10 result discussed
in the article's Impact section (where several correction methods were
within 1% LOO-CV RMSE of one another). Before correction: Bias=48.460,
SD=277.922, 95% LoA=[-496.267, 593.187], r=-0.1088. Auto (Best Model)
selects **Hybrid (Time-Series + Deming)** (RMSE=218.1773), with winning
parameters W1=4, τ_rise=τ_fall=7.0 min, Linnet λ=0.25, deployed as
`CDI_corrected = (CDI_fast - 2184.4811) / -7.0966`. Corrected statistics:
Bias=0.6979, SD=161.9964, giving the software's own displayed "SD ▼41.7%,
LoA ▼41.7%" reduction (reported as ~42% in the article, on unrounded
values) and narrowing the 95% limits of agreement to approximately
[-316.8, 318.2]. The Correction Report's formula text should confirm all
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

- Exact statistics may vary by a small amount (last-digit rounding)
  depending on MATLAB version and floating-point behavior, but should
  match the values above to within normal numerical tolerance.
- Percentage reductions quoted in the article are computed from unrounded
  intermediate values, not from the 3-4 decimal figures printed in the UI
  or reported here; recomputing from the rounded figures alone can give a
  visibly different percentage (e.g. ≈8.5% vs. the reported ≈10% for the
  pH example in §2) without indicating an error.
- See [`ABL_CDI_Analyzer/docs/DATA_FORMATS.md`](ABL_CDI_Analyzer/docs/DATA_FORMATS.md)
  and [`PatLogGUI/docs/DATA_FORMAT.md`](PatLogGUI/docs/DATA_FORMAT.md) for
  the exact input file format specifications.
