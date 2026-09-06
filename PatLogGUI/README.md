# PatLogGUI

A MATLAB GUI tool for importing, cleaning, exploring, and exporting
ABL Flex 800 blood gas data exports.

## Motivation

A single ABL Flex 800 export can contain several hundred columns spanning
patient measurements, calibration cycles, quality-control runs, and
instrument metadata, with multiple patients, catheter sites, and time points
interleaved. PatLogGUI automates the cleaning and provides systematic,
filterable exploration of the resulting dataset.

## Features

- Automatic delimiter (comma/semicolon) and decimal-notation (EU/US)
  detection for CSV and Excel imports
- Essential-column selection (with special handling for pH, temperature,
  and arterial/venous pO2 and sO2); cleaning drops columns only, and no
  rows are removed at this stage
- Three-dimensional filtering: patient identity (single or multi-patient),
  catheter placement location, and clinical protocol time point
- Single- and multi-parameter temporal visualization, with single- and
  multi-patient overlay modes (simple or aligned-timeline)
- Catheter-location comparison plots
- Descriptive statistics (mean, SD, median, min, max) with boxplots for any
  subset of the eleven detected blood gas parameters
- Export of cleaned or patient-filtered datasets to Excel

## Requirements

- MATLAB R2025b (developed and tested in this version; compatibility with
  other MATLAB releases has not been verified)
- **Statistics and Machine Learning Toolbox** — required by the Statistical
  Analysis module, which draws box plots with `boxplot`. All other features
  run on core MATLAB.

## Usage

1. **Import Patient Data File** (CSV or Excel export from the ABL Flex 800).
2. **Process Clinical Data** to clean the import (drops non-essential
   columns; all data rows are kept).
3. Select a **Patient ID** (or enter multiple, comma-separated), **Sample
   Location**, and **Blood Gas Parameter**.
4. Use **Trend Analysis**, **Multi-Parameter View**, or **Location
   Comparison** to visualize; use **Statistical Analysis** for descriptive
   statistics and boxplots.
5. **Export Cleaned Data** or **Export Patient Data** to Excel.

See [docs/DATA_FORMAT.md](docs/DATA_FORMAT.md) for the expected input file
format.
