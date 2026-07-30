# Example data

Real device exports from anesthetized porcine experiments, used for
validation in the accompanying thesis and SoftwareX article. These are
animal experimental data (not human patient data).

| File | Source device | Use |
|---|---|---|
| `PatLog_export.csv` | ABL Flex 800 | Load into **PatLogGUI** directly. Also the ABL-side input for **ABG Analyzer Enhanced** (contains multiple patients, including 2026007 and 2026027 below). |
| `2026007_cdi.log` | Terumo CDI 500 | CDI-side input for **ABG Analyzer Enhanced**, pairs with patient **2026007** in `PatLog_export.csv`. Corresponds to "Dataset 2026007" in the article (small, unstable recording, N=10 paired measurements). |
| `2026027_cdi.log` | Terumo CDI 500 | CDI-side input for **ABG Analyzer Enhanced**, pairs with patient **2026027** in `PatLog_export.csv`. Corresponds to "Dataset 2026027" in the article (larger, stable recording, N=67 paired measurements). |

## Reproducing the PatLogGUI example

1. Launch `PatLogGUI`, click **Import Patient Data File**, select
   `PatLog_export.csv`.
2. Click **Process Clinical Data** to clean it.
3. Select a **Patient ID** (e.g. `2026027`) and a **Blood Gas Parameter**,
   then use **Trend Analysis**, **Multi-Parameter View**, or **Location
   Comparison**.

## Reproducing the ABG Analyzer Enhanced example

1. Launch `ABG_Analyzer_Enhanced`.
2. **Load ABL Data** → `PatLog_export.csv`.
3. **Load CDI Data** → `2026027_cdi.log` (or `2026007_cdi.log`).
4. Select **Patient ID** `2026027` (or `2026007`), a **Parameter** (e.g.
   `pH`), and click **Analyze**.
5. Select **Auto (Best Model)** under Correction Method and click **Apply
   Correction** to reproduce the auto-selected correction results reported
   in the article (e.g. for Dataset 2026027, pH: Hybrid correction, bias
   reduced to ~0.000, SD reduced by ~56%).

See `../ABG_Analyzer_Enhanced/docs/DATA_FORMATS.md` and
`../PatLogGUI/docs/DATA_FORMAT.md` for the exact file-format specifications
these examples follow.
