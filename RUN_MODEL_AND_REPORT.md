# Run the weekly model and build the report

The entry point is `run_weekly_model_report.m`. It reads the current configured CSV, rebuilds the weekly panel and empirical-Bayes priors, runs one joint 14-variable correlation model, generates the Bayesian correlation paths, and compiles `output/pdf/weekly_correlation_report.pdf` with those paths included.

## Requirements

- MATLAB R2024b with the Statistics and Machine Learning and Econometrics toolboxes.
- Python 3 with NumPy, pandas, Matplotlib, and Pillow.
- A TeX installation containing `latexmk` and `pdflatex`.
- Poppler for the report checks (`pdfinfo`, `pdftotext`, and `pdftoppm`).

Build the accelerated correlation kernel once for the installed MATLAB version and platform:

```matlab
build_dsc_mex()
```

The default `CorrelationBackend="auto"` uses that MEX file when it is available and otherwise uses the slower MATLAB implementation.

## Basic run

Start MATLAB in the repository root and run:

```matlab
result = run_weekly_model_report();
```

The defaults request 10 warm-up iterations followed by 100 retained draws, keep every draw, allow two hours for the sampler, use four native correlation threads, run the MATLAB tests, generate all correlation-path figures, compile the report, and render-check the report.

The same run from a terminal is:

```sh
/Applications/MATLAB_R2024b.app/bin/matlab -batch "result=run_weekly_model_report();"
```

## Choosing warm-up, draws, and thinning

Warm-up iterations are discarded. `RetainedDraws` is the number of saved post-warm-up draws. `Thin` saves one draw after every specified number of completed post-warm-up iterations. The requested number of sampler sweeps is

```text
WarmupIterations + RetainedDraws * Thin
```

For example, 100 warm-up iterations and 1,000 retained draws with no thinning require 1,100 completed sweeps:

```matlab
result = run_weekly_model_report( ...
    WarmupIterations=100, ...
    RetainedDraws=1000, ...
    Thin=1, ...
    MaxHours=14);
```

Based on the measured 100-iteration run on this machine, one sweep took about 41 seconds. Therefore, 1,100 sweeps would take roughly 12.5 hours before figure and report generation. Runtime varies with the machine, backend, thread count, and proposal behavior; set `MaxHours` above the expected sampler time.

Useful configurations are:

| Purpose | Warm-up | Retained draws | Thin | Suggested time cap |
|---|---:|---:|---:|---:|
| Pipeline smoke test | 2 | 3 | 1 | 0.25 hours |
| Preliminary path figure | 10 | 100 | 1 | 2 hours |
| Initial longer chain | 100 | 1,000 | 1 | 14 hours |

The smoke test verifies that the complete pipeline runs, but its posterior summaries are not substantively useful. The 10/100 configuration is suitable for checking the appearance and behavior of the paths. A fixed draw count cannot guarantee convergence. For substantive posterior claims, run independently initialized chains and assess rank-normalized split R-hat, effective sample size, Monte Carlo error, and trace plots. Increase warm-up or retained draws when those diagnostics require it.

Thinning reduces stored output but does not reduce sampler work. Keep `Thin=1` unless storage is a demonstrated problem. Each retained draw stores the 91 distinct entries of the weekly correlation matrices, plus the mean and log-volatility states.

## Parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `WarmupIterations` | `10` | Completed initial sweeps that are discarded. |
| `RetainedDraws` | `100` | Saved draws after warm-up; must be at least 2 to create bands. |
| `Thin` | `1` | Save every nth completed post-warm-up sweep. |
| `MaxHours` | `2` | Sampler wall-time limit. Figure and report time is additional. |
| `ChunkSize` | `10` | Retained draws per MAT-v7.3 posterior chunk. |
| `Seed` | `20260914` | Base random-number seed. |
| `ChainID` | `1` | Chain identifier; the effective seed is `Seed + ChainID - 1`. |
| `CorrelationBackend` | `"auto"` | `"auto"`, `"mex"`, or `"matlab"`. |
| `CorrelationThreads` | `4` | Native date-block threads used by the MEX correlation calculation. |
| `SourceFile` | configured CSV | Optional path to another source CSV with the expected layout. |
| `OutputRoot` | `outputs/weekly_research` | Numerical data, run directories, tests, and latest-run pointer. |
| `RunTests` | `true` | Run the MATLAB acceptance suite before estimation. |
| `BuildReport` | `true` | Build the PDF report after generating the paths. |
| `VerifyReport` | `true` | Inspect report structure and render every PDF page. |
| `PythonExecutable` | `"python3"` | Python executable used by the publication scripts. |

Example with explicit compute settings:

```matlab
result = run_weekly_model_report( ...
    WarmupIterations=10, ...
    RetainedDraws=100, ...
    MaxHours=2, ...
    CorrelationBackend="mex", ...
    CorrelationThreads=4, ...
    ChunkSize=10, ...
    Seed=20260914, ...
    ChainID=1);
```

Use `CorrelationThreads=4` on the machine used for the current benchmark. Raising it is not automatically faster; the best value depends on physical CPU cores and memory bandwidth. When running several chains simultaneously, reduce the threads per chain so the total number of native threads does not substantially exceed the available cores.

## Outputs

Each run is preserved under

```text
outputs/weekly_research/runs/TIMESTAMP-chainN-wW-rR-thinK/
```

The run directory contains:

- `run_manifest.json`: data hash, configuration, requested draw counts, MATLAB version, and code provenance.
- `checkpoint.mat`: complete sampler state and random-number state.
- `posterior_chunk_*.mat`: retained draws.
- `pilot_summary.json`: actual iterations, draws, timing, stop reason, and implementation details.
- `posterior_correlation_bands.mat`: pointwise 16th, 50th, and 84th percentiles of the actual correlations.
- `figures/`: one selected-pair PDF and 11 pages covering all 91 pairs.

The report is written to `output/pdf/weekly_correlation_report.pdf`. It includes the selected Bayesian paths and all 91 Bayesian paths from the same run. The vertical background colors reproduce Figure 1's continuous scale: blue for negative correlation, white near zero, and red for positive correlation, with intensity determined by the posterior median.

`outputs/weekly_research/latest_model_report.json` records the latest run directory, summary, figure directory, and report path.

If `MaxHours` stops the sampler after at least two retained draws, the pipeline still creates the figures and report and records the shortfall. If fewer than two draws are retained, it preserves the checkpoint and stops with an error because percentile bands cannot be computed.

## Rebuild the report without rerunning the model

To rebuild from a saved run, use the run's own summary and path figures:

```sh
python3 scripts/build_publication.py \
  --results outputs/weekly_research \
  --pilot outputs/weekly_research/runs/RUN_NAME/pilot_summary.json \
  --posterior-run outputs/weekly_research/runs/RUN_NAME \
  --report-only

python3 scripts/check_publication.py --report-only
```

This publication command reads saved results only. It does not start MATLAB or rerun MCMC.
