# Running the tests

The repository ships 23 automated tests. They check that the software still
produces every value [WALKTHROUGH.md](WALKTHROUGH.md) says it reports, which
includes every result the accompanying article quotes for the two example
recordings, and that the file parsers, temporal pairing, MAD filtering and
correction arithmetic behave as documented.

Running them takes about four minutes and needs no setup.

## What you need

- **MATLAB R2025b.** The expected values are those of this release; other
  versions have not been verified.
- Nothing else. The tests add the source folders to the MATLAB path
  themselves, and use only the example data already in `examples/`.

---

## Step 1 — Get the code

In a terminal:

```
git clone https://github.com/goga-lanchava/abl-cdi-method-comparison-tools
```

Or, without git: open the
[v1.1.0 release page](../../releases/tag/v1.1.0), click **Source code (zip)**
at the bottom, and unzip it.

Either way you end up with a folder containing `ABL_CDI_Analyzer`,
`PatLogGUI`, `examples` and `tests`.

## Step 2 — Open that folder in MATLAB

Start MATLAB and navigate the **Current Folder** panel into the folder from
Step 1 — the one that contains `tests`. You should see `README.md`,
`WALKTHROUGH.md`, `examples` and `tests` listed.

From the Command Window instead, if you prefer:

```matlab
cd 'C:\path\to\abl-cdi-method-comparison-tools'
```

## Step 3 — Run them

In the Command Window:

```matlab
runtests('tests')
```

That is the whole command. There is nothing to configure first.

## Step 4 — While it runs

A dot appears for each test that finishes:

```
Running ComponentTest
..........
Done ComponentTest
```

**Application windows will open and close on their own during the run.**
This is expected: the tests click the real buttons in both programs rather
than calling internal functions behind the scenes. Do not close them or
click anything — they close themselves.

## Step 5 — Read the result

The run ends with:

```
Totals:
   23 Passed, 0 Failed, 0 Incomplete.
   206.03 seconds testing time.
```

`23 Passed, 0 Failed` means every figure quoted in `WALKTHROUGH.md` was
reproduced, along with all the parsing, pairing, filtering and correction
checks.

---

## If a test fails

The output names the exact value that disagrees, for example:

```
Verification failed in WalkthroughTest/section2_autoCorrection.
    Test Diagnostic:
    SD of paired differences no longer matches WALKTHROUGH.md
    Framework Diagnostic:
    verifyEqual failed.
    --> The error was not within absolute tolerance.
    --> Failure table:
            Actual    Expected    Error     AbsoluteTolerance
            0.0431     0.0426     0.0005          5e-05
    Stack Information:
    In ...\tests\WalkthroughTest.m (WalkthroughTest.section2_autoCorrection) at 118
```

so you can see which figure moved and by how much. A failure does not stop
the remaining tests — you always get the complete picture from one run.

### One failure that is expected on some machines

`toolboxDependenciesMatchDocumentation` asserts that `PatLogGUI` requires the
Statistics and Machine Learning Toolbox (it uses `boxplot`) and that
`ABL_CDI_Analyzer` does not. **Without that toolbox installed this single
test is expected to fail** — that is the test correctly reporting your
environment, not a fault in the code. The other 22 should still pass.

---

## Running only part of the suite

```matlab
runtests('tests/ComponentTest.m')                        % ~1 minute
runtests('tests/WalkthroughTest.m')                      % ~3 minutes
runtests('tests/WalkthroughTest.m', 'ProcedureName', 'section2_autoCorrection')
```

---

## What the tests check

### `tests/ComponentTest.m` — 10 deterministic component tests

| Test | Checks |
|---|---|
| `ablParserReadsKnownFile` | `PatLog_export.csv` yields 1907 rows, keeps `Time`/`PatientID`/`pH`/`pO2`/`pCO2`, lists both example patients, and parses every timestamp |
| `cdiParserReadsKnownFile` | `2026027_cdi.log` yields 3235 rows, the exact 9-column set, and the exact first and last timestamps |
| `cdiParserDropsEntirelyMissingColumns` | of the 18 fixed CDI fields only the 8 carrying data are kept; `VO2` is absent |
| `cdiParserSkipsMalformedLines` | lines with no bracketed date, an invalid device time, or only placeholders (`%`, `---`, `--`, `-.-`) are discarded; valid lines survive |
| `ablParserRejectsMalformedTimestamps` | a file whose first column is not a timestamp raises an error |
| `ablParserToleratesMissingPatientColumn` | a file with no patient-ID column still parses, returning an empty patient list |
| `temporalPairingRespectsTolerance` | 68 pairs at a 5-minute tolerance, 67 at 1 minute |
| `madFilterExcludesGrossOutlierOnly` | on differences `1…10` plus `1000` (median 5.5, MAD 2.5, band ±11.25) only the outlier is dropped |
| `madFilterKeepsEverythingWhenMadIsZero` | identical differences give MAD = 0, and no pair is rejected |
| `biasCorrectionMatchesHandComputedOffset` | the fitted offset, the corrected series, the `NaN`-ing of excluded pairs and the zero residual bias all match a value the test recomputes independently |

### `tests/WalkthroughTest.m` — 13 reproduction and regression tests

| Test | Checks |
|---|---|
| `section2_beforeCorrection` | §2 pre-correction panel: 67/68 pairs, bias 0.047, SD 0.095, LoA, r; the fitting window Auto selects |
| `section2_fitWindowIsRequired` | without the fitting window the §2 figures do **not** appear, and the values the walkthrough warns you would get instead (bias 0.049, the 1.8666/−6.2502 model, 62 robust pairs, SD ▼12.8%) do |
| `section2_autoCorrection` | Weighted Deming wins; coefficients, λ, 61/67 retained pairs, after-correction bias, SD and r; on the 61 retained pairs, SD 0.0475 → 0.0426 (▼10.3%) and LoA [−0.020, 0.166] → [−0.084, 0.083] |
| `section2_tieBreakerEngages` | the top two candidates fall inside the 1% RMSE band: RMSE 0.0953 vs 0.0957 (0.42%), LoA spans 0.3628 vs 0.3646 |
| `section3_beforeCorrection` | §3 pre-correction panel for pO2 |
| `section3_autoCorrection` | Hybrid wins; RMSE, W1, τ, λ, deployed fit, after-correction bias, SD, r, SD ▼41.8% and LoA [−316.747, 316.983] |
| `section3_noTieBreaker` | the candidates are well separated (runner-up RMSE 292.9552, 34.2% behind), so the tie-breaker does not engage |
| `section3_allPatientsGivesSameResult` | the walkthrough's note that "All Patients" gives identical statistics |
| `looCvTuningIsFoldIndependent` | the LOO-CV result is reproducible and no candidate scores `NaN` |
| `exportResultsWritesFiles` | Export Results writes both the `.xlsx` and the two `.csv` files |
| `patLogKeepsEighteenEssentialColumns` | cleaning keeps exactly 18 essential columns, removes 218, and leaves all 1907 rows on display |
| `patLogReadsLatin1Encoding` | the Latin-1 example import renders `T (°C)` rather than a replacement character |
| `toolboxDependenciesMatchDocumentation` | the declared per-tool toolbox requirements are accurate |

---

## Why the expected values are fixed numbers

`WALKTHROUGH.md` commits the software to exact statistics for the two example
recordings, and the article quotes them. That makes the walkthrough a
specification, so the tests assert those values literally rather than
checking that the code merely runs. If a change shifts a result — as the
leave-one-out fix in v1.1.0 did for the pO2 example — the suite fails and
names the figure that moved, which is the intended behaviour.
