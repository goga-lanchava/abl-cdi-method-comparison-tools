# ABG analyzes and ABL–CDI measurement methods comparison and correction Tools

Two MATLAB App Designer applications for working with blood gas data from the
ABL Flex 800 (intermittent electrochemical analyzer) and the Terumo CDI 500
(continuous optical monitor):

| Tool | Purpose |
|---|---|
| [`PatLogGUI`](PatLogGUI/) | Import, clean, filter, visualize, and export ABL Flex 800 exports |
| [`ABL_CDI_Analyzer`](ABL_CDI_Analyzer/) | Temporally align ABL and CDI data, compute Bland–Altman method-comparison statistics, and fit/select correction models |

See each tool's own README for details, requirements, and usage instructions.

## Requirements

- MATLAB R2025b (developed and tested in this version; compatibility with
  other MATLAB releases has not been verified)
- No additional toolboxes required beyond core MATLAB

Compiled standalone Windows executables (no
MATLAB license is needed to run them) are published under the
[Executables release](../../releases/tag/v1).

## Reproducing the article's results

See [WALKTHROUGH.md](WALKTHROUGH.md) for step-by-step instructions that use
the data in `examples/` to reproduce the representative figures and
statistics reported in the accompanying SoftwareX article.

## Repository structure

```
.
├── PatLogGUI/
│   ├── src/PatLogGUI.m
│   └── docs/
├── ABL_CDI_Analyzer/
│   ├── src/ABL_CDI_Analyzer.m
│   └── docs/DATA_FORMATS.md
├── examples/            # sample/illustrative ABL + CDI data
├── WALKTHROUGH.md       # reproduces the article's reported figures/statistics
├── CHANGELOG.md
├── CODE_METADATA.md     # SoftwareX submission metadata table (both tools)
├── LICENSE
└── README.md
```

## Citation

If you use this software in your research, please cite the accompanying
SoftwareX article (citation details to be added upon publication).

## License

See [LICENSE](LICENSE). 
