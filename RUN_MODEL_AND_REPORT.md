# Run the weekly model and build the report

The entry point is `run_weekly_model_report.m`. It reads the current configured CSV, builds the weekly panel, estimates or loads the model parameters, generates the Bayesian correlation paths, and compiles the publication outputs. The main report is `output/pdf/weekly_correlation_report.pdf`; the complete correlation-path figure set is written separately to `output/pdf/weekly_correlation_paths.pdf`.

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

The saved parameter file contains the retained posterior draws of `V` (the covariance of mean-state innovations), `sig2h` (log-variance innovation variances), and `sig2r` (correlation-coordinate innovation variances). It also stores their posterior means as a convenience summary, plus the calibrated priors, the exact training dates and returns, configuration, and run provenance. New files use schema 4: each draw retains its original chain ID and iteration, alongside per-chain counts, seeds, and the estimation diagnostic report. Posterior means weight every retained draw equally, including when chain lengths differ. Schema-3 draw files remain loadable, but their convergence status is `not_checked`; older mean-only formats require re-estimation. It is written atomically with a UTC date/time and unique suffix, for example:

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

The current report build keeps the methodology report short: appendices and sample path pages are excluded from `weekly_correlation_report.pdf`. The full set of Bayesian correlation paths is written separately to `weekly_correlation_paths.pdf`.

Latest saved run used in the report: four smoothing chains were requested with 20 warm-up iterations and 500 retained draws per chain. The stage reached the 12-hour time cap after 908 completed iterations across chains and 828 pooled retained draws. The offline convergence diagnostic is reported by default in Section 5 and currently fails: retained counts differ by chain (`209`, `206`, `204`, `209`), so the aligned diagnostic uses 204 draws per chain and convergence is not established.

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

With an earlier estimation cutoff, the workflow runs **two stages**: parameter estimation through the cutoff, followed by conditional smoothing through the latest week. Each stage uses the specified warm-up, retained-draw target, and `MaxHours` budget; the combined sampling budget can therefore be twice `MaxHours`. With no earlier cutoff and the default `ParameterSmoothing="mean"`, `EstimateParameters=true` estimates parameters first, then runs a second smoothing stage using their posterior means. Explicit `ParameterSmoothing="draws"` uses a single joint estimation and its joint posterior draws for the figures when no earlier cutoff is set. In report/off mode, a time-capped partial parameter draw set can be saved if every requested chain has retained chunks and there are at least two retained draws overall; actual counts are recorded. Convergence is not established by saving a file. Numerical failures and warm-up-only runs cannot publish parameter files.

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

The defaults request **one sequential chain**, with 10 warm-up iterations followed by 100 retained draws, keep every draw, allow 12 hours per sampling stage, use four native correlation threads, run the MATLAB tests, generate all correlation-path figures, compile the report, and render-check the report. `ConvergenceMode="report"` records diagnostics without blocking the workflow, and `AutoExtend=false` keeps the requested draw target. This short single-chain default cannot receive a passing convergence status.

The same run from a terminal is:

```sh
/Applications/MATLAB_R2024b.app/bin/matlab -batch "result=run_weekly_model_report();"
```

## Run with or without parameter estimation

Estimate parameters and then smooth conditional on posterior-mean parameters:

```matlab
result = run_weekly_model_report( ...
    EstimateParameters=true, ...
    ParameterSmoothing="mean", ...
    NumChains=4, ParallelChains=true, ...
    WarmupIterations=20, RetainedDraws=500, ...
    MaxHours=12, ConvergenceMode="report");
```

Reuse the latest saved parameter file and run only the conditional smoothing stage:

```matlab
result = run_weekly_model_report( ...
    EstimateParameters=false, ...
    ParameterSmoothing="mean", ...
    NumChains=4, ParallelChains=true, ...
    WarmupIterations=20, RetainedDraws=500, ...
    MaxHours=12, ConvergenceMode="report");
```

`ConvergenceMode="report"` is the default and should be left on for report-producing runs. It records rank-normalized split/folded R-hat, bulk and tail ESS, and MCSE/SD diagnostics in the run folder, then prints the status in Section 5 of the report. Use `ConvergenceMode="off"` only for smoke tests or debugging runs where convergence reporting is intentionally not needed.

## Choosing warm-up, draws, and thinning

Warm-up iterations are discarded **separately in each chain**. `RetainedDraws` is the number of saved post-warm-up draws **per chain**. `Thin` saves one draw after every specified number of completed post-warm-up iterations. The requested number of sampler sweeps per chain is

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

## Independent chains and convergence checks

Use `NumChains=4` to compare four independently initialized chains. Each chain has its own seed, warm-up, checkpoint, and retained draws. Chains are kept separate for diagnostics and pooled only afterward for parameter means and correlation bands. The same options apply to parameter estimation and, by default, the subsequent posterior-mean-parameter smoothing stage.

For example, this requests four chains **one after another**, without Parallel Computing Toolbox:

```matlab
result = run_weekly_model_report( ...
    EstimateParameters=true, ParameterSmoothing="mean", ...
    NumChains=4, ParallelChains=false, ...
    WarmupIterations=100, RetainedDraws=1000, ...
    ConvergenceMode="report", CheckEvery=250, ...
    AutoExtend=false, MaxHours=60);
```

The warm-up and draw numbers are starting choices, not evidence of convergence. Four chains with 1,000 retained draws give 4,000 pooled draws per stage. A longer single chain can be inspected with `NumChains=1`, but it cannot assess agreement between independent chains and cannot pass this workflow's four-chain requirement. Split-chain halves do not count as additional independent chains.

To run those chains simultaneously, change `ParallelChains=true`. **This option requires MATLAB Parallel Computing Toolbox and an available license.** The program uses a process pool with at least `NumChains` workers, creating one if none exists, following [MathWorks' worker execution API](https://www.mathworks.com/help/parallel-computing/parallel.pool.parfeval.html). It rejects an existing thread pool or an undersized process pool; close or resize that pool before retrying. Without the toolbox, use `ParallelChains=false`; the statistical diagnostics are the same. There is no silent fallback from parallel to sequential execution.

Each process holds its own model state. `CorrelationThreads` is the number of native calculation threads **per chain**, not across all chains. For example, four chains with `CorrelationThreads=1` use four native correlation threads in total; with `CorrelationThreads=4`, they may use sixteen, in addition to MATLAB's own work. Parallel execution consumes more memory and is not guaranteed to be four times faster. The existing process pool is left available after the run.

### What is checked

The code records rank-normalized split/folded R-hat, bulk and tail effective sample sizes (ESS), and Monte Carlo standard error (MCSE) for the derived correlation paths `P_pairs` and log-volatility paths `h`. The latent correlation coordinates `r`, transition paths `B`, and hyperparameters `V`, `sig2h`, and `sig2r` are excluded from the convergence verdict. The default passing criteria are:

- At least four original chains and at least 100 retained draws in every chain.
- R-hat strictly below `1.01` for every checked quantity.
- Bulk and tail ESS each at least `400` across the chains.
- Mean MCSE divided by the sample standard deviation at most `0.05`.
- Equal retained counts across chains, with no unavailable or failing diagnostics.

During estimation and smoothing, checks cover `P(i,j,t)` and `h(j,t)` at **every modeled week**. The transition paths, hyperparameters, and latent `r` are not part of the convergence verdict. Diagnostics follow the [Stan posterior definitions](https://mc-stan.org/posterior/reference/diagnostics.html) and [rank-normalized diagnostic methodology](https://doi.org/10.1214/20-BA1221).

These are diagnostic criteria, not a mathematical guarantee that the posterior has been explored. Inspect chain behavior and the quantities relevant to your application as well. Numerical validity of the correlation matrices is a separate issue. With unequal lengths, diagnostics use the earliest common retained length for exploratory comparison, but the run cannot pass because pooled summaries also use the extra draws.

### Report diagnostics or switch them off

| `ConvergenceMode` | Effect |
|---|---|
| `"off"` | Skip diagnostic calculations; save status `not_checked`. |
| `"report"` (default) | Save diagnostic results and continue with available valid retained draws, even if criteria fail. Figures/report disclose the status. |

In report mode, failed or unavailable diagnostics do not block estimation, parameter loading, smoothing, or report generation. A loaded file's estimation status is retained in the output, and the new smoothing chains receive their own diagnostics. A saved file is selected by its estimation timestamp; specify `ParameterFile` explicitly to choose another.

`ParameterSmoothing="mean"` remains the default. When smoothing conditional on previously estimated `"draws"`, the sampler cycles parameter values between sweeps, so that stage is marked `not_checked`; the workflow still continues. This limitation does **not** apply to a single-stage joint estimation through the latest week using `EstimateParameters=true, ParameterSmoothing="draws"`.

### Automatic extension and time limits

`CheckEvery` specifies retained draws **per chain** between checkpoints for convergence assessment. An assessment is made after all chains reach a batch boundary and have at least `MinDiagnosticDraws`, and again at the end when feasible. Warm-up is not repeated between batches. The program always attempts the initial `RetainedDraws` target even if an early check passes.

To continue beyond that target when checks do not pass, add:

```matlab
% Include these options in a run_weekly_model_report(...) call:
% AutoExtend=true, MaxRetainedDraws=4000, CheckEvery=250
```

It stops when checks pass after the initial target, the per-chain `MaxRetainedDraws` cap is reached, or `MaxHours` expires. `MaxRetainedDraws` must be at least `RetainedDraws` when extending; otherwise it has no effect. `AutoExtend` requires checking to be enabled. With fewer than four chains, extension cannot produce a passing status and will reach a cap instead. It adds retained draws only; it does not increase warm-up automatically. If initial transient behavior remains, start a new run with more warm-up.

`MaxHours` is a **shared wall-time budget for an entire stage**, not a fresh budget for every chain or batch. Pool startup, sequential chains, and diagnostic calculations use that budget. Parameter estimation and conditional smoothing have separate stage budgets; with the default mean mode, total sampling can approach twice `MaxHours`. Time checks occur at safe boundaries; a calculation or file write already underway can finish after the limit. A diagnostic assessment interrupted by the time budget is marked `not_checked`, never passed. Figure/report generation is additional.

All state-path checks can take appreciable time and disk space. Diagnostics stage one posterior chunk at a time into temporary MAT arrays and process bounded blocks; temporary staging needs approximately another copy of the retained data and is cleaned up afterward. The current figure builder pools correlation draws in memory, so allow memory for `modeled weeks × 91 pairs × total retained draws × 8 bytes`, plus working arrays.

Automatic batches resume the complete saved state and RNG stream. For manual continuation, use the lower-level `dsc_run_chains(panel,priors,cfg,run_dir)` with the same modeled panel, priors, chain IDs, seed, warm-up, and thinning, `cfg.resume=true`, and a larger total `cfg.max_seconds` and/or draw target. `cfg.max_iterations` is the per-chain sweep target. The recorded stage time is not replenished on resume. The high-level `run_weekly_model_report` always creates a new timestamped run; invoking it again does not resume the last run.

## Parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `EstimateParameters` | `true` | Estimate and save parameters; `false` loads saved parameters for conditional smoothing. |
| `EstimationEndDate` | `""` | Estimation cutoff (`yyyy-MM-dd`); empty uses the latest week. Leave empty when loading. |
| `ParameterDirectory` | `OutputRoot/parameters` | Directory for timestamped parameter files. |
| `ParameterFile` | `""` | Optional exact saved file to load; empty selects the latest. Only valid with estimation off. |
| `ParameterSmoothing` | `"mean"` | `"mean"` uses posterior-mean parameters; `"draws"` explicitly uses saved parameter draws in smoothing. |
| `NumChains` | `1` | Independent original chains per stage; use at least 4 for a passing convergence assessment. |
| `ParallelChains` | `false` | Run chains simultaneously with Parallel Computing Toolbox; false runs sequentially. |
| `ConvergenceMode` | `"report"` | `"off"` skips diagnostic calculations; `"report"` records results and continues regardless of diagnostic status. |
| `AutoExtend` | `false` | Continue beyond the initial retained target until diagnostic success or a cap. |
| `CheckEvery` | `250` | Retained draws per chain between diagnostic batch boundaries. |
| `MaxRetainedDraws` | `2000` | Per-chain retained cap when `AutoExtend=true`. |
| `MinDiagnosticDraws` | `100` | Minimum retained draws per chain for a passing assessment; supported minimum is 6. |
| `RhatThreshold` | `1.01` | Every checked R-hat must be strictly below this value. |
| `MinESS` | `400` | Required bulk and tail ESS for each quantity, combining the chains. |
| `MaxMCSERatio` | `0.05` | Maximum mean MCSE divided by the sample standard deviation for each quantity. |
| `WarmupIterations` | `10` | Completed initial sweeps discarded in each chain. |
| `RetainedDraws` | `100` | Saved draws per chain after warm-up; must be at least 2 to create bands. |
| `Thin` | `1` | Save every nth completed post-warm-up sweep. |
| `MaxHours` | `12` | Shared wall-time limit per stage, including all chains and diagnostics. Figure/report time is additional. |
| `ChunkSize` | `10` | Retained draws per MAT-v7.3 posterior chunk. |
| `Seed` | `20260914` | Base random-number seed. |
| `ChainID` | `1` | First chain identifier; subsequent IDs increase by one. Effective seed is `Seed + ID - 1`. |
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
- `checkpoint.mat` and `posterior_chunk_*.mat`: complete sampler/RNG state and retained draws for a single-chain stage. For multiple chains these live under `chains/chain_001/`, `chains/chain_002/`, and so on, each with its own summary/result.
- `pilot_summary.json`: actual iterations, draws, timing, stop reason, and implementation details.
- `convergence_diagnostics.json`, `.csv`, and `.mat`: overall status/thresholds/provenance and per-quantity R-hat, ESS, MCSE, availability, and pass/fail details. Parameter estimation has its own copies under `estimation/` when a separate estimation stage is used.
- `inference_metadata.json`: estimation/smoothing mode, loaded or saved parameter file, estimation timestamp, training cutoff, and smoothing end date.
- `calibrated_priors.mat` and `prior_summary.json`: the selected run's own prior snapshot, used when rebuilding its report even after a different cutoff has been run.
- `estimation/`: separate parameter-estimation stage whenever conditional smoothing follows, including the default posterior-mean mode.
- `posterior_correlation_bands.mat`: pointwise 16th, 50th, and 84th percentiles of the actual correlations.
- `figures/`: one selected-pair PDF and 11 pages covering all 91 pairs.

The main report is written to `output/pdf/weekly_correlation_report.pdf`. The complete set of Bayesian correlation paths from the same run is written to `output/pdf/weekly_correlation_paths.pdf`. The path-file vertical background colors reproduce Figure 1's continuous scale: blue for negative correlation, white near zero, and red for positive correlation, with intensity determined by the posterior median.

`outputs/weekly_research/latest_model_report.json` records the latest run directory, summary, figure directory, and report path.
It also records the inference metadata. In MATLAB, inspect `result.inference` and `result.paths.parameters` for the selected parameter file and estimation date.

In `"report"`/`"off"` mode, a time-limited run can produce preliminary output if enough valid draws exist. Parameter saving additionally requires retained chunks from every requested chain; a chain that never started cannot be silently omitted. A numerical failure stops the workflow. Completed checkpoints are preserved in all these cases.

Inspect `result.convergence` for the final state-chain assessment, `result.inference.estimation_convergence` for the saved estimation assessment, and `result.saved_draws_per_chain` for actual state-chain counts. The report records estimation and state-smoothing statuses separately.

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

To check real parallel-worker execution separately (requires Parallel Computing Toolbox):

```matlab
checks = runtests('tests/integration/test_parallel_chains.m');
assertSuccess(checks);
```

This compares two synthetic parallel chains with their sequential equivalents, including retained draws and checkpoint RNG states. It closes only a pool created by the test. Ordinary acceptance tests do not start a process pool; the missing-toolbox branch is skipped on machines where the toolbox is installed.
