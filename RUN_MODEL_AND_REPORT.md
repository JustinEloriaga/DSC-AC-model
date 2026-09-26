# Run the weekly model and build the report

The entry point is `run_weekly_model_report.m`. It reads the current configured CSV, builds the weekly panel, estimates or loads the model parameters, generates the Bayesian correlation paths, and compiles `output/pdf/weekly_correlation_report.pdf` with those paths included.

## Estimate parameters or load a saved model

`EstimateParameters` is the switch. Its default is `true`.

Estimate using observations through a chosen date, save the parameter estimates, and smooth correlations through the latest available week:

```matlab
result = run_weekly_model_report( ...
    EstimateParameters=true, ...
    EstimationEndDate="2020-12-25", ...
    WarmupIterations=10, RetainedDraws=100, MaxHours=2);
```

The first 104 weekly returns are reserved for initial-prior calibration and **excluded from both parameter estimation and smoothing**. Estimation begins at week 105 and ends at the selected cutoff. For the current data, the calibration period is 15 August 2003–5 August 2005, and the remaining **1,100 weekly returns, from 12 August 2005 to 4 September 2026**, form the model sample when no earlier cutoff is selected. Descriptive tables still summarize all 1,204 returns.

A date between Fridays selects the last available Friday on or before it; cutoffs within the reserved calibration period or outside the available return period are rejected. Leave `EstimationEndDate=""` to estimate through the latest available week. The 104-week rolling calculation for the evolution priors still needs at least 106 total weekly returns and sufficient complete observations: the first 104 are excluded from the likelihood, leaving at least two modeled weeks. This is a minimum data requirement, not the actual sample size.

Initial-state moments use the reserved 104 weeks. The existing evolution-prior calibration continues to use rolling windows through the estimation cutoff, including some observations that also enter estimation. Thus the initial 104-week period is held out of the likelihood, while evolution priors retain their empirical-Bayes calibration. No observations after the cutoff enter either prior calibration or parameter estimation.

The saved parameter file contains the retained posterior draws of `V` (the covariance of mean-state innovations), `sig2h` (log-variance innovation variances), and `sig2r` (correlation-coordinate innovation variances). It also stores their posterior means as a convenience summary, plus the calibrated priors, the exact training dates and returns, configuration, and run provenance. It is written atomically with a UTC date/time and unique suffix, for example:

```text
outputs/weekly_research/parameters/dsc_parameters_20260926T143012123Z_6d4a92e1-8b73-4f26-a015-9c72e680b3fd.mat
```

After adding new observations to the data in `data/yad_tickers_no_usdcnh_from_20030804.csv`, smooth with the latest saved parameters:

```matlab
result = run_weekly_model_report( ...
    EstimateParameters=false, ...
    ParameterSmoothing="mean", ...
    WarmupIterations=10, RetainedDraws=100, MaxHours=2);
```

The program prints the loaded file, the UTC date/time when its parameters were estimated, and the last weekly return used in that estimation. "Latest" is determined by the saved estimation timestamp, not the file's modification time. Smoothing runs do not create new parameter files. An empty parameter directory produces an error explaining that an estimation run is needed first. Legacy checkpoint or prior MAT files are not parameter files for this interface.

To use a particular saved model or share parameters between different output directories:

```matlab
result = run_weekly_model_report( ...
    EstimateParameters=false, ...
    ParameterDirectory="/absolute/path/to/parameters", ...
    ParameterFile="/absolute/path/to/parameters/dsc_parameters_TIMESTAMP_UNIQUE.mat");
```

`ParameterFile` is optional. The loader checks the variable order, pair order, parameter dimensions, and the dates, returns, and observation masks in both the reserved calibration block and the subsequent estimation block. Additional weeks are allowed. Changed or missing historical observations require a fresh estimation run. The saved calibration length determines where smoothing starts, even if the current configuration differs. The original saved priors are reused without calibrating them on later observations. A corrupt latest artifact is reported visibly instead of silently selecting an older model. Parameter files from the earlier version that included calibration weeks in the likelihood are rejected with an instruction to re-estimate; their original files are preserved.

`ParameterSmoothing` controls how saved parameters enter smoothing. The default `ParameterSmoothing="mean"` uses the posterior means of `V`, `sig2h`, and `sig2r` as fixed plug-in values, so the plotted bands condition on those means and exclude parameter uncertainty. `ParameterSmoothing="draws"` explicitly reuses the retained parameter draws, so the reported correlation bands include variation across those saved parameter draws. Parameter files retain their full draws in either mode. In both cases the latent mean, volatility, and correlation paths are sampled jointly from the first week after the reserved calibration period through the latest week. This includes the requested period from the estimation cutoff through the latest week, and earlier modeled estimates can also change. The calibration weeks never re-enter the likelihood. This is historical smoothing, not a real-time filter. Smoothing still needs warm-up and retained draws, and most correlation calculations remain, so it need not be substantially faster.

With an earlier estimation cutoff, the workflow runs **two stages**: parameter estimation through the cutoff, followed by conditional smoothing through the latest week. Each stage uses the specified warm-up, retained-draw target, and `MaxHours` budget; the combined sampling budget can therefore be twice `MaxHours`. With no earlier cutoff and the default `ParameterSmoothing="mean"`, `EstimateParameters=true` estimates parameters first, then runs a second smoothing stage using their posterior means. Explicit `ParameterSmoothing="draws"` uses a single joint estimation and its joint posterior draws for the figures when no earlier cutoff is set. If estimation stops at the time cap with at least two retained draws, the partial parameter draw set can be saved and its actual draw count is recorded; convergence is not established by saving a file. Numerical failures and warm-up-only runs cannot publish parameter files.

The lower-level configuration equivalents are `cfg.estimate_parameters`, `cfg.estimation_end_date`, `cfg.parameter_dir`, `cfg.parameter_file`, and `cfg.parameter_smoothing` in `weekly_config.m`. `cfg.prior_weeks` controls the number of initial calendar weeks reserved and excluded (104 by default). Use the configurable report entry point above for estimation and saved-parameter smoothing; the bounded `run_weekly_research('pilot',cfg)` also excludes those initial weeks.

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

Thinning reduces stored output but does not reduce sampler work. Keep `Thin=1` unless storage is a demonstrated problem. Each retained draw stores the 91 distinct entries of the weekly correlation matrices, the mean and log-volatility states, and `V`, `sig2h`, and `sig2r`.

## Parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `EstimateParameters` | `true` | Estimate and save parameters; `false` loads saved parameters for conditional smoothing. |
| `EstimationEndDate` | `""` | Estimation cutoff (`yyyy-MM-dd`); empty uses the latest week. Leave empty when loading. |
| `ParameterDirectory` | `OutputRoot/parameters` | Directory for timestamped parameter files. |
| `ParameterFile` | `""` | Optional exact saved file to load; empty selects the latest. Only valid with estimation off. |
| `ParameterSmoothing` | `"mean"` | `"mean"` uses posterior-mean parameters; `"draws"` explicitly uses saved parameter draws in smoothing. |
| `WarmupIterations` | `10` | Completed initial sweeps that are discarded. |
| `RetainedDraws` | `100` | Saved draws after warm-up; must be at least 2 to create bands. |
| `Thin` | `1` | Save every nth completed post-warm-up sweep. |
| `MaxHours` | `2` | Wall-time limit per sampler stage. Figure and report time is additional. |
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
- `inference_metadata.json`: estimation/smoothing mode, loaded or saved parameter file, estimation timestamp, training cutoff, and smoothing end date.
- `calibrated_priors.mat` and `prior_summary.json`: the selected run's own prior snapshot, used when rebuilding its report even after a different cutoff has been run.
- `estimation/`: separate parameter-estimation run when its cutoff precedes the latest week.
- `posterior_correlation_bands.mat`: pointwise 16th, 50th, and 84th percentiles of the actual correlations.
- `figures/`: one selected-pair PDF and 11 pages covering all 91 pairs.

The report is written to `output/pdf/weekly_correlation_report.pdf`. It includes the selected Bayesian paths and all 91 Bayesian paths from the same run. The vertical background colors reproduce Figure 1's continuous scale: blue for negative correlation, white near zero, and red for positive correlation, with intensity determined by the posterior median.

`outputs/weekly_research/latest_model_report.json` records the latest run directory, summary, figure directory, and report path.
It also records the inference metadata. In MATLAB, inspect `result.inference` and `result.paths.parameters` for the selected parameter file and estimation date.

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

## Development checks

`run_weekly_acceptance()` runs the regular regression tests. The more expensive synthetic CSV-to-figures check is separate so normal runs do not repeatedly export its 24 test PDFs:

```matlab
checks = runtests('tests/integration/test_model_report.m');
assertSuccess(checks);
```

It verifies both estimation and loading through the public entry point using temporary synthetic data. It does not estimate the actual market sample or overwrite the published report.
