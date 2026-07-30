# ABL–CDI Method Comparison Tools

Two MATLAB App Designer applications for working with blood gas data from the
ABL Flex 800 (intermittent electrochemical analyzer) and the Terumo CDI 500
(continuous optical monitor):

| Tool | Purpose |
|---|---|
| [`PatLogGUI`](PatLogGUI/) | Import, clean, filter, visualize, and export ABL Flex 800 exports |
| [`ABG_Analyzer_Enhanced`](ABG_Analyzer_Enhanced/) | Temporally align ABL and CDI data, compute Bland–Altman method-comparison statistics, and fit/select correction models |

See each tool's own README for details, requirements, and usage instructions.

## Requirements

- MATLAB R2025b or later (App Designer)
- No additional toolboxes required beyond core MATLAB, based on static review
  of the source — please confirm in a clean MATLAB install before relying on
  this for your own environment.

Both tools are also distributed as standalone Windows executables (see
[Releases](../../releases)) that bundle the MATLAB Runtime, so they can be
run without a MATLAB license.

## Repository structure

```
.
├── PatLogGUI/
│   ├── src/PatLogGUI.m
│   └── docs/
├── ABG_Analyzer_Enhanced/
│   ├── src/ABG_Analyzer_Enhanced.m
│   └── docs/DATA_FORMATS.md
├── examples/            # sample/illustrative ABL + CDI data
├── CHANGELOG.md
├── CODE_METADATA.md     # SoftwareX submission metadata table (both tools)
├── LICENSE
└── README.md
```

## Citation

If you use this software in your research, please cite the accompanying
SoftwareX article (citation details to be added upon publication).

## License

See [LICENSE](LICENSE). *(Currently set to MIT as a placeholder — confirm
the intended license before publishing.)*

## Contact

George (Goga) Lanchava — issues and questions via this repository's issue
tracker.
