# DSC-AC Weekly Model Guide

This file is the practical map of the repository: what each file group does,
how to run the model, how to run the new four-variable Model 1, and which
options control the workflow.

## Main Entry Points

Run commands from the repository root in MATLAB, or use MATLAB's `-batch`
option from a terminal.

```matlab
result = run_weekly_model_report();
```

This is the default end-to-end workflow. It prepares the weekly panel, runs
acceptance checks, estimates or loads parameters, samples the latent paths,
builds Bayesian correlation paths, and optionally builds/checks the report.

```matlab
result = run_weekly_model1_report();
```

This runs the same workflow on Model 1 only:

```matlab
["NKYTR Index", "SPXT Index", "SX5T Index", "TUKXG Index"]
```

These are variables 10, 11, 12, and 13 in the canonical 14-variable data set.
The Model 1 wrapper writes to `outputs/weekly_research_model1` by default and
sets `BuildReport=false` and `VerifyReport=false`, because the current
publication report is written for the 14-variable baseline. It also defaults to
`ParameterSmoothing="draws"` so the reported paths use the joint
parameter/state draws rather than a second fixed-parameter smoothing pass.

```matlab
main
```

This is a lightweight analysis entry point. It does not run the full configurable
report pipeline.

```matlab
results = run_weekly_research('analysis');
results = run_weekly_research('pilot');
```

`analysis` prepares the panel, calibrates priors, and writes descriptive output.
`pilot` runs the bounded low-iteration sampler configured in `weekly_config.m`.

## Basic Runs

Default full 14-variable run:

```matlab
result = run_weekly_model_report();
```

Four-variable Model 1 run:

```matlab
result = run_weekly_model1_report( ...
    WarmupIterations=10, ...
    RetainedDraws=100, ...
    MaxHours=2);
```

Run the same four variables without the wrapper:

```matlab
result = run_weekly_model_report( ...
    VariableNames=["NKYTR Index","SPXT Index","SX5T Index","TUKXG Index"], ...
    OutputRoot="outputs/weekly_research_model1", ...
    BuildReport=false, ...
    VerifyReport=false);
```

Run from the terminal:

```sh
/Applications/MATLAB_R2024b.app/bin/matlab -batch "result=run_weekly_model_report();"
```

## Variables

The default is all supported variables in canonical order:

```matlab
VariableNames="all"
```

The canonical full set is:

```text
1   AUDUSD Index
2   BCOMCOT Index
3   BCOMGCTR Index
4   EURUSD Index
5   GBPUSD Index
6   I02981JP Index
7   LEATTREU Index
8   LSG1TRGU Index
9   LUATTRUU Index
10  NKYTR Index
11  SPXT Index
12  SX5T Index
13  TUKXG Index
14  USDJPY Index
```

You can pass any unique supported subset with at least two variables:

```matlab
result = run_weekly_model_report( ...
    VariableNames=["SPXT Index","SX5T Index"]);
```

The requested order is preserved. Saved parameter files are checked against the
variable order, so changing the variable set usually requires a fresh estimation
run rather than loading an old 14-variable parameter file.

## Public Options

These options belong to `run_weekly_model_report`. The Model 1 wrapper forwards
the same options except that it fixes `VariableNames` to the four Model 1
tickers and defaults `ParameterSmoothing` to `"draws"`.

| Option | Default | Meaning |
|---|---:|---|
| `WarmupIterations` | `10` | Initial sweeps discarded in each chain. |
| `RetainedDraws` | `100` | Saved post-warm-up draws per chain. |
| `Thin` | `1` | Save every nth post-warm-up sweep. |
| `MaxHours` | `12` | Wall-time cap per sampling stage. |
| `ChunkSize` | `10` | Retained draws per posterior chunk file. |
| `Seed` | `20260914` | Base random seed. |
| `ChainID` | `1` | First chain identifier. |
| `NumChains` | `1` | Number of independent chains. Use at least 4 for convergence assessment. |
| `ParallelChains` | `false` | Run chains simultaneously using Parallel Computing Toolbox. |
| `ConvergenceMode` | `"report"` | `"report"` computes diagnostics and continues; `"off"` skips them. |
| `AutoExtend` | `false` | Continue adding retained draws until diagnostics pass or a cap/time limit stops the run. |
| `CheckEvery` | `250` | Retained draws per chain between diagnostic checks. |
| `MaxRetainedDraws` | `2000` | Per-chain cap when `AutoExtend=true`. |
| `MinDiagnosticDraws` | `100` | Minimum retained draws per chain for diagnostics. |
| `RhatThreshold` | `1.01` | Maximum allowed R-hat. |
| `MinESS` | `400` | Minimum bulk and tail effective sample size. |
| `MaxMCSERatio` | `0.05` | Maximum mean MCSE divided by posterior standard deviation. |
| `CorrelationBackend` | `"mex"` | `"mex"`, `"matlab"`, or `"auto"`. |
| `CorrelationThreads` | `4` | Native date-block threads per chain for the MEX backend. |
| `SourceFile` | `""` | Optional alternative CSV path. Empty uses `weekly_config.m`. |
| `OutputRoot` | `""` | Optional alternative output directory. Empty uses `weekly_config.m`. |
| `VariableNames` | `"all"` | `"all"` or a string array of supported ticker names. |
| `EstimateParameters` | `true` | Estimate parameters before smoothing. If `false`, load saved parameters. |
| `EstimationEndDate` | `""` | Optional estimation cutoff date. Empty uses the latest available return. |
| `ParameterDirectory` | `""` | Directory containing saved parameter files. Empty uses `OutputRoot/parameters`. |
| `ParameterFile` | `""` | Exact saved parameter file to load. Empty selects the latest. |
| `KeepNewestParameterFiles` | `3` | Keep only this many newest saved parameter files at the end of the run. Set `0` to disable cleanup. |
| `KeepIntermediateArtifacts` | `false` | Delete regenerable per-chain summaries, runtime JSON handoffs, cached correlation bands, checkpoints, and Finder metadata after the run. Run provenance and convergence JSON are retained. Set `true` only for debugging. |
| `ParameterSmoothing` | `"mean"` | `"mean"` uses posterior mean parameters; `"draws"` propagates saved parameter draws. Model 1 defaults this option to `"draws"`. |
| `RunTests` | `true` | Run acceptance checks before the model. |
| `BuildReport` | `true` | Build the publication report after sampling. |
| `VerifyReport` | `true` | Verify the generated PDF report. Requires `BuildReport=true`. |
| `PythonExecutable` | `"python3"` | Python executable used by publication scripts. |

## Common Run Recipes

Smoke test:

```matlab
result = run_weekly_model_report( ...
    WarmupIterations=2, ...
    RetainedDraws=3, ...
    MaxHours=0.25, ...
    RunTests=false, ...
    BuildReport=false, ...
    VerifyReport=false);
```

Model 1 smoke test:

```matlab
result = run_weekly_model1_report( ...
    WarmupIterations=2, ...
    RetainedDraws=3, ...
    MaxHours=0.25, ...
    RunTests=false);
```

Estimate through an earlier cutoff, then smooth through the latest data:

```matlab
result = run_weekly_model_report( ...
    EstimateParameters=true, ...
    EstimationEndDate="2020-12-25", ...
    ParameterSmoothing="mean", ...
    WarmupIterations=10, ...
    RetainedDraws=100, ...
    MaxHours=2);
```

Load the latest saved parameters and smooth only:

```matlab
result = run_weekly_model_report( ...
    EstimateParameters=false, ...
    ParameterSmoothing="mean", ...
    WarmupIterations=10, ...
    RetainedDraws=100, ...
    MaxHours=2);
```

Longer four-chain diagnostic run:

```matlab
result = run_weekly_model_report( ...
    NumChains=4, ...
    ParallelChains=false, ...
    WarmupIterations=100, ...
    RetainedDraws=1000, ...
    ConvergenceMode="report", ...
    MaxHours=60);
```

Use the compiled MEX correlation and historical-mean smoothing kernels:

```matlab
build_dsc_mex();
result = run_weekly_model_report(CorrelationBackend="mex");
```

Both kernels are already built on this device. Rebuild them with
`build_dsc_mex` after cloning on another device. The Model 1 benchmark
reduced 100 retained draws from 252.8 to 33.1 seconds; see
[the performance report](research/sampler_performance_100_draws.md).

Convergence reporting checks the derived correlation paths `P_pairs` and
log-volatility paths `h` at every modeled date. The latent correlation
coordinates `r`, transition paths `B`, and estimated hyperparameters `V`,
`sig2h`, and `sig2r` remain available where saved but are excluded from the
convergence pass/fail decision.

## File Map

### Root MATLAB Files

| File | Purpose |
|---|---|
| `weekly_config.m` | Central defaults: source file, output folder, sampler controls, parameter mode, backend settings, and default variable set. |
| `run_weekly_model_report.m` | Main public workflow for the model and report pipeline. |
| `run_weekly_model1_report.m` | Convenience wrapper for variables 10-13: NKYTR, SPXT, SX5T, TUKXG. |
| `run_weekly_research.m` | Lower-level analysis/pilot workflow used by the main runner. |
| `run_weekly_pilot.m` | Bounded sampler pilot helper. |
| `run_weekly_acceptance.m` | Acceptance checks for the MATLAB workflow. |
| `build/build_dsc_mex.m` | Builds the native correlation and historical-mean smoothing MEX kernels. |

### Core Model Functions

| File | Purpose |
|---|---|

### Toolbox

| File | Purpose |
|---|---|
| `toolbox/prepare_weekly_panel.m` | Reads levels, creates weekly Friday as-of returns, applies variable selection, and builds pair indices. |
| `toolbox/analyze_weekly_panel.m` | Descriptive analysis, stationarity checks, serial-dependence checks, rolling correlations, and panel summaries. |
| `toolbox/dsc_prepare_inference.m` | Splits calibration/estimation/smoothing periods and prepares priors/inference metadata. |
| `toolbox/dsc_run_inference.m` | Coordinates estimation and smoothing stages. |
| `toolbox/dsc_run_chains.m` | Runs one or more chains, sequentially or in parallel. |
| `toolbox/dsc_sample.m` | Main sampler loop. |
| `toolbox/dsc_mean_mex.cpp` | Native Kalman filter and backward mean-path sampler; MATLAB supplies random noise. |
| `toolbox/dsc_volatility_likelihood.m` | Exact cached Gaussian likelihood change for a volatility coordinate. |
| `toolbox/dsc_save_parameters.m` | Saves and loads timestamped parameter files with identity checks. |
| `toolbox/dsc_load_parameters.m` | Loads saved parameter artifacts. |
| `toolbox/dsc_calibrate_priors.m` | Empirical-Bayes prior calibration helper. |
| `toolbox/dsc_mcmc_diagnostics.m`, `toolbox/dsc_diagnose_chains.m`, `toolbox/dsc_validate_convergence.m` | Convergence and diagnostic calculations. |
| `toolbox/dsc_credible_bands.m` | Computes pointwise 16/50/84 percent posterior bands from retained draws. |
| `toolbox/dsc_corr_matrix.m`, `toolbox/dsc_corr_batch_mex.cpp` | Correlation matrix construction and optional native acceleration. |
| `toolbox/dsc_random_walk_draw.m` | Random-walk proposal helper. |
| `toolbox/dsc_panel_slice.m` | Slices panel structs over date ranges. |
| `toolbox/dsc_parallel_available.m` | Checks whether process-based parallel chains can run. |
| `toolbox/dsc_sampler_signature.m` | Hash/signature metadata for sampler provenance. |

### Scripts

| File | Purpose |
|---|---|
| `scripts/plot_bayes_correlation_paths.m` | Builds Bayesian correlation path figures from posterior chunks. |
| `scripts/build_publication.py` | Builds publication tables, figures, and report assets from saved evidence. |
| `scripts/check_publication.py` | Checks generated publication PDF structure/rendering. |

### Tests

| File | Purpose |
|---|---|
| `tests/test_weekly_panel.m` | Tests weekly panel construction, calendar handling, and variable subset selection. |
| `tests/test_sampler.m` | Tests sampler behavior on synthetic inputs. |
| `tests/test_chain_runner.m` | Tests chain execution and checkpoints. |
| `tests/test_convergence_diagnostics.m`, `tests/test_convergence_metadata.m` | Tests diagnostic calculations and metadata. |
| `tests/test_correlation_kernel.m` | Tests correlation backend behavior. |
| `tests/test_credible_bands.m` | Tests posterior band calculations. |
| `tests/test_inference_workflow.m` | Tests estimation/smoothing workflow logic. |
| `tests/test_parameter_files.m` | Tests saved parameter schemas and identity checks. |
| `tests/test_bayes_paths.m` | Tests Bayesian path plotting behavior. |
| `tests/test_publication_inference.py` | Tests publication inference metadata handling. |
| `tests/verify_exports.py` | Independent export reconciliation for saved outputs. |
| `tests/test_model_report.m` | Public entry-point integration test on synthetic data. |
| `tests/test_parallel_chains.m` | Parallel-chain integration test. |

### Data, Research, and Outputs

| Path | Purpose |
|---|---|
| `data/yad_tickers_no_usdcnh_from_20030804.csv` | Main 14-variable source CSV used by default. |
| `data/yad_tickers.csv` | Source CSV including USDCNH. |
| `research/` | LaTeX report/slides sources, bibliography, derivations, and performance notes. |
| `outputs/weekly_research/` | Default numerical output directory. |
| `outputs/weekly_research_model1/` | Default numerical output directory for `run_weekly_model1_report`. |
| `output/pdf/` | Publication PDF destination used by the report builder. |

## Outputs Created by a Run

A normal run creates a timestamped directory such as:

```text
outputs/weekly_research/runs/YYYYMMDD-HHMMSS-SSS-chain1-w10-r100-thin1/
```

Important files inside or near that run are:

| File | Meaning |
|---|---|
| `run_manifest.json` inside each run directory | Configuration, data hash, requested draws, MATLAB version, Git/source provenance for that run. |
| `pilot_summary.json` | Actual iterations, retained draws, stop reason, timings, backend details. |
| `checkpoint.mat` | Saved sampler state and RNG state. |
| `posterior_chunk_*.mat` | Retained posterior draws. |
| `posterior_correlation_bands.mat` | 16th/50th/84th percentile correlation bands. |
| `figures/` | Pairwise Bayesian correlation path PDFs. |
| `convergence_diagnostics.*` | Diagnostic JSON/CSV/MAT files when convergence reporting is enabled. |
| `inference_metadata.json` | Estimation/smoothing mode and parameter-file metadata. |
| `parameters/dsc_parameters_*.mat` | Saved parameter files, when parameters are estimated. |
| `latest_model_report.json` | Pointer to the latest run, including its manifest, under an output root. |

## Requirements

- MATLAB R2024b.
- MATLAB Statistics and Machine Learning Toolbox.
- MATLAB Econometrics Toolbox.
- Parallel Computing Toolbox only if `ParallelChains=true`.
- Python 3 with NumPy, pandas, Matplotlib, and Pillow for publication scripts.
- TeX with `latexmk` and `pdflatex` for report building.
- Poppler for PDF checks.

## Development Checks

Run MATLAB tests:

```matlab
checks = runtests('tests');
assertSuccess(checks);
```

Run the public integration test:

```matlab
checks = runtests('tests/test_model_report.m');
assertSuccess(checks);
```

Run only the variable-subset tests:

```matlab
checks = runtests('tests/test_weekly_panel.m', ...
    'ProcedureName','testSelectedVariablesUseCanonicalSubset');
assertSuccess(checks);
```

The current local data hash may differ from an older expected hash in
`tests/test_weekly_panel.m`. That is separate from the new variable-selection
logic.

## Practical Notes

- `RetainedDraws` is per chain. Four chains with 1,000 retained draws produce
  up to 4,000 retained draws before any time cap stops the run.
- `MaxHours` applies per sampling stage. A separate estimation stage and
  smoothing stage can each use that budget.
- `ParameterSmoothing="mean"` conditions smoothing on posterior mean
  parameters and is the default for the general runner.
- `ParameterSmoothing="draws"` propagates saved parameter draw variation.
- `run_weekly_model1_report` defaults to `ParameterSmoothing="draws"`.
- Full-report publication currently assumes the 14-variable baseline.
- Model 1 is ready for model output and figures, but its report text should be
  generalized before using `BuildReport=true`.
