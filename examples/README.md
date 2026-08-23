# Example data

Real device exports from anesthetized porcine experiments, used for
validation in the accompanying thesis and SoftwareX article. These are
animal experimental data (not human patient data).

**De-identification:** the ABL Flex 800 export's name/birthdate/note
fields are blank in the original device output for these recordings — no
personal or institutional identifiers were present to remove. Patient
identity is limited to the device-assigned numeric Patient ID.

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
4. Select **Patient ID** `2026027` (or `2026007`), a **Parameter** (e.g.
   `pH`), and click **Analyze**.
5. Select **Auto (Best Model)** under Correction Method and click **Apply
   Correction** to reproduce the auto-selected correction results reported
   in the article (e.g. for Dataset 2026027, pH: Linnet weighted-Deming
   correction, λ=0.10, bias reduced to ~0.000, SD reduced by ~10% on the
   61/67 MAD-retained pairs — see `../WALKTHROUGH.md` §2 for the exact
   figures and where each one appears in the software's output).

See `../ABL_CDI_Analyzer/docs/DATA_FORMATS.md` and
`../PatLogGUI/docs/DATA_FORMAT.md` for the exact file-format specifications
these examples follow.
