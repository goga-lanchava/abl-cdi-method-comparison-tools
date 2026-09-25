# Example data

Real device exports from anesthetized porcine experiments, used as the
illustrative examples in the accompanying thesis and SoftwareX article.
These are animal experimental data (not human patient data).

**Contents and identifiers:** `PatLog_export.csv` is a complete ABL Flex 800
export of 1,907 measurements from many experiments, not only the two
recordings used below. In this laboratory the analyzer's `Last Name` and
`First Name` fields are used to record the protocol time point (e.g. `BL1`,
`MAP1 60`) and the sampling site (e.g. `CVC`, `ECMO *v(`), and PatLogGUI
reads them as such (see
[`../PatLogGUI/docs/DATA_FORMAT.md`](../PatLogGUI/docs/DATA_FORMAT.md)); they
contain no personal names. The `Patient Id` field holds the experiment
number, in some rows with a free-text sample label. The birthdate,
patient-note and physician fields are empty, and the operator field holds
the device default `Anonymous`.

See [`../WALKTHROUGH.md`](../WALKTHROUGH.md) for step-by-step instructions
that reproduce the article's representative figures and reported statistics
using this data.

| File | Source device | Use |
|---|---|---|
| `PatLog_export.csv` | ABL Flex 800 | Load into **PatLogGUI** directly. Also the ABL-side input for **ABL-CDI Analyzer** (contains multiple patients, including 2026007 and 2026027 below). |
| `2026007_cdi.log` | Terumo CDI 500 | CDI-side input for **ABL-CDI Analyzer**, pairs with patient **2026007** in `PatLog_export.csv`. Corresponds to "Dataset 2026007" in the article (small, unstable recording, N=10 paired measurements). |
| `2026027_cdi.log` | Terumo CDI 500 | CDI-side input for **ABL-CDI Analyzer**, pairs with patient **2026027** in `PatLog_export.csv`. Corresponds to "Dataset 2026027" in the article (larger, stable recording, N=67 paired measurements). |

## Reproducing the PatLogGUI example

1. Launch `PatLogGUI`, click **Import Patient Data File**, select
   `PatLog_export.csv`.
2. Click **Process Clinical Data** to clean it.
3. Select a **Patient ID** (e.g. `2026027`) and a **Blood Gas Parameter**,
   then use **Trend Analysis**, **Multi-Parameter View**, or **Location
   Comparison**.

## Reproducing the ABL-CDI Analyzer example

1. Launch `ABL_CDI_Analyzer`.
2. **Load ABL Data** → `PatLog_export.csv`.
3. **Load CDI Data** → `2026027_cdi.log` (or `2026007_cdi.log`).
4. Select **Patient ID** `2026027` (or `2026007`), a **Parameter** (`pH` for
   2026027, `pO2` for 2026007), and a **Time Tolerance** of 5 minutes.
5. **For 2026027 only:** tick the **Fit Window** checkbox and click its
   **Auto** button. This step is required: without it the analysis uses all
   68 pairs instead of the 67 window pairs, and the article's values do not
   reproduce. Leave Fit Window unticked for 2026007.
6. Click **Analyze**.
7. Select **Auto (Best Model)** under Correction Method and click **Apply
   Correction** to reproduce the auto-selected correction results reported
   in the article (e.g. for Dataset 2026027, pH: Linnet weighted-Deming
   correction, λ=0.10, bias reduced to -0.0002, SD 0.0475 → 0.0426
   (▼10.3%) on the 61/67 MAD-retained pairs — see `../WALKTHROUGH.md` §2
   and §3 for the exact figures and where each one appears in the
   software's output).

See `../ABL_CDI_Analyzer/docs/DATA_FORMATS.md` and
`../PatLogGUI/docs/DATA_FORMAT.md` for the exact file-format specifications
these examples follow.
