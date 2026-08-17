# ABL-CDI Analyzer

A MATLAB App Designer tool for aligning, comparing, and correcting continuous
blood-gas monitor readings (CDI) against intermittent arterial blood gas
(ABG) reference measurements (ABL), with multiple statistical
correction/calibration methods and automated model selection.

## Motivation

Continuous in-line blood gas monitors (e.g. CDI systems) drift relative to
intermittent, gold-standard ABG samples (e.g. ABL systems). This tool
time-aligns the two data streams, quantifies the agreement/bias between them
(Bland–Altman, correlation, regression), and applies a correction model so
the continuous signal tracks the reference measurements more closely.

## Features

- ABL export CSV parser and CDI monitor log parser
- Nearest-neighbour time alignment with MAD-based outlier filtering and
  optional auto time-shift detection for clock offset correction
- Six correction methods (Bias, OLS, Proportional, Deming (error-variance
  ratio λ, default λ=1, user-configurable), a simplified Passing–Bablok fit
  (median pairwise-slope regression; does not include the confidence-interval
  or linearity-test procedures of the full standard method), and a Hybrid
  Time-Series+Deming method with Auto-Tune) plus Auto (nested Leave-One-Out
  Cross-Validation) model selection, with the narrower Bland–Altman
  limits-of-agreement span used as a tie-breaker when the top two
  candidates' RMSE differ by less than 1%
- Fitting-window controls with a stability score
- Time trend, correlation, Bland–Altman, and before/after comparison
  visualizations
- Export of corrected data, figures (SVG), and LaTeX-formatted correction
  formulas
- Before/after summary uses a neutral, software-defined descriptive label
  (e.g. "Bias + SD improved") rather than a clinical-sounding verdict

## Requirements

- MATLAB R2025b (developed and tested in this version; compatibility with
  other MATLAB releases has not been verified)
- No additional toolboxes identified as required by static inspection;
  please confirm against your MATLAB installation before first use.

## Usage

1. **Load ABL** and **Load CDI** data files.
2. Select a **Patient ID** and **Parameter**.
3. Adjust **time tolerance**/**time shift** (or use **Auto Shift**) to align
   the two streams, then click **Analyze**.
4. Choose a **Correction Method** (or Auto for LOO-CV selection), optionally
   restrict to a **fitting window**, and click **Apply Correction**.
5. **Export** the corrected trend, figures, or formula.

See [docs/DATA_FORMATS.md](docs/DATA_FORMATS.md) for input file format
details and [../CODE_METADATA.md](../CODE_METADATA.md) for the SoftwareX
submission metadata table.
