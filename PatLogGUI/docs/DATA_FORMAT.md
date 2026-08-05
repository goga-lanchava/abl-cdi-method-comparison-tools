# Input data format

PatLogGUI ingests ABL Flex 800 PatLog exports (`.csv` or `.xlsx`). This spec
is derived directly from the parsing logic in `src/PatLogGUI.m`
(`cleanDataFunc`, `enableAnalysis`).

## Column selection

From the imported file's header row, PatLogGUI keeps only these columns
(case-insensitive substring match, except where noted):

`Time`, `Sample #`, `Patient Id`, `Last Name`, `First Name`, `Sample type`,
`Patient sex`, `pH` (exact match on `"pH"`/`pH`), `pO2` (excludes any header
containing `pO2(v)`), `pCO2`, `T` (read from a fixed column position —
column 26 — to accommodate non-standard header formatting), `tHb`, `sO2`
(excludes any header containing `sO2(v)`), `Lac`, `K+`, `Na+`, `Cl-`,
`FIO2`.

All other columns (calibration/QC metadata, venous-side duplicates, etc.)
are dropped. The app reports how many essential columns were kept and how
many were removed.

## Row / field semantics

- **Patient ID** — the column matching `Patient Id`/`PatientId`/`Patient
  ID`/`PatientID` (case-insensitive). Values are coerced to numeric; rows
  with a non-numeric or zero Patient ID do not populate the patient
  dropdown.
- **Catheter placement location** — read from the `First Name` column
  field.
- **Clinical protocol time point** — read from the `Last Name` column field
  (e.g. `MAP 40 NOR ANO`, `BL 1`).
- **Parameter columns** — matched by case-insensitive substring against a
  fixed list of eleven blood gas parameter names; Temperature is matched by
  a stricter pattern (`T (` … `C)`) after stripping surrounding quotes.

## Notes

- Delimiter/decimal-notation handling is performed by MATLAB's import
  functions when reading the CSV/Excel file; no manual configuration is
  required by the user.
