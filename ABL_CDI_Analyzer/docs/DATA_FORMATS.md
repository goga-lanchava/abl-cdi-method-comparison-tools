# Input data formats

ABL-CDI Analyzer ingests two file types. This spec is derived directly
from the parsing logic in `src/ABL_CDI_Analyzer.m` (`parseABL`,
`parseCDI`).

## 1. ABL export (reference / gold-standard measurements)

- **Format:** delimited text (`.csv`), delimiter auto-detected as comma or
  semicolon (whichever occurs more often in the header row)
- **Header row required.** Recognized/expected columns:
  - A timestamp column (first column), format `d.M.yyyy HH:mm`
  - A patient identifier column — any of `Patient Id`, `PatientId`,
    `Patient_Id`, or any header containing both "patient" and "id"
  - Any number of numeric parameter columns (e.g. pH, pCO2, pO2, etc.)
- **Automatically excluded** non-measurement columns: `Error code`, patient
  name/status/sex/birthdate/notes/department/accession/site/draw
  time/physician/operator/approval/report layout/measuring mode/errors
  detected fields, and any column with no numeric values
- Column names are cleaned (parenthetical units stripped, quotes removed);
  duplicate column names are de-duplicated by keeping the version with more
  non-zero numeric values
- Decimal commas (`,`) are converted to decimal points before parsing

## 2. CDI monitor log (continuous signal)

- **Format:** plain text log, UTF-8 (BOM tolerated), lines beginning with a
  bracketed date, e.g. `[yyyy-MM-dd HH:mm:ss.SSS] HH:mm:ss,val1,val2,...`
- Delimiter auto-detected between comma and tab based on the first bracketed
  line
- Fixed column set (in order after the device time field):
  `pH, pCO2, pO2, TEMP, HCO3, BE, sO2, K+, VO2, Q, BSA, pH_v, pCO2_v, pO2_v,
  TEMP_v, SO2_v, HCT, tHb`
- Placeholder values (`%`, `---`, `--`, `-.-`) are treated as missing
- Rows without a valid device time (`HH:mm:ss`) or without at least one
  valid numeric field are discarded
- A parameter column is only kept in the parsed table if at least one row
  has a non-missing value for it

## Notes for the SoftwareX illustrative example

For the manuscript's illustrative example, use a small, de-identified or
synthetic sample of each file type (a handful of patients/timepoints is
sufficient) — see `examples/`.
