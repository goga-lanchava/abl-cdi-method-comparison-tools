# ABG analyzes and ABL–CDI measurement methods comparison and correction tools

Two MATLAB GUI applications for working with blood gas data from the
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
- `ABL_CDI_Analyzer` — core MATLAB only
- `PatLogGUI` — additionally requires the **Statistics and Machine Learning
  Toolbox**, used by the Statistical Analysis module's box plots
  (`boxplot`). Every other feature runs on core MATLAB.

On MATLAB releases older than R2025b, `ABL_CDI_Analyzer` may also require
the Statistics and Machine Learning Toolbox: it calls `prctile`, which
moved into core MATLAB only in recent releases.

Standalone Windows installers are published under the
[v1.1.0 release](../../releases/tag/v1.1.0). No MATLAB license is needed to
run the applications. These are MATLAB Compiler *web installers*: each is a
few megabytes and downloads the matching MATLAB Runtime (~1.2 GB, shared
between the two tools) during a one-time installation, so internet access is
required the first time you install.

## Reproducing the article's results

See [WALKTHROUGH.md](WALKTHROUGH.md) for step-by-step instructions that use
the data in `examples/` to reproduce the representative figures and
statistics reported in the accompanying SoftwareX article.

## Tests

Open this folder in MATLAB and run:

```matlab
runtests('tests')
```

That is the whole setup — the tests put the source folders on the path
themselves. Expect 23 passing tests in roughly three to four minutes.

Two suites:

- **`tests/ComponentTest.m`** — ten deterministic component tests: parsing of
  the known ABL and CDI example files, rejection of malformed timestamps,
  tolerance of a missing patient-ID column, skipping of malformed CDI log
  lines, temporal pairing at two tolerances, MAD filtering on hand-checkable
  artificial data (including the degenerate zero-MAD case), and a Bias
  Correction whose fitted offset and corrected series are checked against a
  value the test recomputes independently.
- **`tests/WalkthroughTest.m`** — asserts the exact statistics quoted in
  [WALKTHROUGH.md](WALKTHROUGH.md) against the data in `examples/`, so the
  software and the article cannot drift apart unnoticed. It also covers the
  export fix, the import encoding, and the declared toolbox dependencies.

To run one suite or one test:

```matlab
runtests('tests/ComponentTest.m')
runtests('tests/WalkthroughTest.m', 'ProcedureName', 'section2_autoCorrection')
```

Notes for anyone reproducing the results:

- The expected values are those of MATLAB R2025b. Other releases have not
  been verified, and small numerical differences would show up here first.
- Application windows open and close during the run. That is normal, and no
  interaction is needed — the tests drive the real callbacks directly.
- `toolboxDependenciesMatchDocumentation` asserts that `PatLogGUI` requires
  the Statistics and Machine Learning Toolbox and that `ABL_CDI_Analyzer`
  does not. Without that toolbox installed this test is expected to fail,
  which is itself the correct signal about the environment.

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
├── tests/               # regression tests for the WALKTHROUGH.md figures
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
