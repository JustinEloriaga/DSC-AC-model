# Weekly dynamic correlation research

Reproducible weekly analysis of 14 market series, excluding USDCNH. The first-stage package contains 52/104-week rolling correlations, data diagnostics, a bounded joint Bayesian sampler pilot, a LaTeX report and Beamer slides.

## Quick start

For the configurable end-to-end workflow that runs the model, creates the Bayesian correlation paths, and compiles the report, see [RUN_MODEL_AND_REPORT.md](RUN_MODEL_AND_REPORT.md). The entry point separates warm-up iterations from retained draws:

```matlab
result = run_weekly_model_report(WarmupIterations=10,RetainedDraws=100,MaxHours=2);
```

MATLAB R2024b with Statistics and Machine Learning and Econometrics toolboxes was verified locally. From the repository root:

```matlab
main                                      % data analysis and prior calibration only
results = run_weekly_research('analysis'); % same workflow, explicit entrypoint
checks = run_weekly_acceptance();          % tests plus saved JSON evidence
results = run_weekly_research('pilot');    % at most 20 sweeps or 30 minutes
pilot_checks = verify_weekly_pilot(results.pilot.run_dir); % no further MCMC
```

The pilot is explicitly **not a converged posterior estimate**. Its default 20 iterations are all warm-up. A full multi-chain production run is not part of the first-stage authorization.

Defaults are in [weekly_config.m](weekly_config.m). For a closure sensitivity analysis, use a separate output directory so the baseline remains intact:

```matlab
cfg = weekly_config();
cfg.closure_policy = 'mask';
cfg.output_root = fullfile(pwd,'outputs','weekly_research_masked');
sensitivity = run_weekly_research('analysis',cfg);
```

On this Mac, MATLAB is installed at `/Applications/MATLAB_R2024b.app/bin/matlab` but is not necessarily on PATH. Example:

```sh
/Applications/MATLAB_R2024b.app/bin/matlab -batch "checks=runtests('tests'); assertSuccess(checks);"
```

## Data and meaning

The original source and existing Excel workbook are preserved. The configured three-header CSV contains **levels**, not returns. Each series uses its own last published price on or before each completed Friday. Percent log returns are computed as `100 * diff(log(weekly_levels))`.

The supplied snapshot produces 1,205 Friday levels and **1,204 weekly returns from 2003-08-15 to 2026-09-04**, across 14 variables and 91 distinct pairs. The unfinished week ending 2026-09-11 is excluded. Native currency quotes remain unchanged.

Every endpoint retains its source date and quality flags. Routine carry age is at most four calendar days. NKYTR on 2019-05-03 is the single documented seven-day exception, during Japan's extended market closure. Baseline last-published marks produce zero/catch-up returns. Masked sensitivity excludes the two affected NKYTR returns without deleting calendar weeks. No future backfilling or interpolated prices are used.

Rolling correlations are trailing sample correlations. They begin only after the complete 52/104-observation window is available. The Bayesian target is a different estimand: historical, full-sample-smoothed **conditional innovation correlations** after time-varying means.

## Estimator and reproducibility

The new path uses the existing DSC-SV-AH model family with `p=0`: 14 time-varying means, 14 log-variance states, 91 unrestricted matrix-log correlation states, and a joint positive-definite correlation matrix. It does not estimate 91 separate bivariate models.

Initial empirical-Bayes moments use the first 104 returns with 5% diagonal covariance shrinkage. Evolution scales use rolling covariance estimates over the full historical panel, comparing 52/104/156 weeks and choosing 104. This data reuse is deliberate and documented. All 1,204 dates remain in estimation.

The new sampler uses linear-memory random-walk draws, observed-subspace likelihoods, cached correlation factorizations for volatility updates, explicit numerical failures and atomic completed-sweep checkpoints. It stores chain-specific post-burn-in packed correlations, not full matrix histories. Warm-up diagnostics are labeled as such.

Each analysis saves configuration, source SHA-256, MATLAB version, Git revision/dirty status and a hash of the MATLAB source tree. Each sampler run has its own directory, full state/RNG checkpoint and timing/proposal diagnostics. Historical runs are not deleted. To resume, load that run's saved panel/configuration/priors and call `dsc_sample` with `cfg.resume=true` and the same directory. Fixed settings must match. The pilot wrapper continues to enforce its cap.

Legacy `tvsvar_modified_msv2_gam2_gen` and original MARS support functions remain for historical reproducibility. They are not used by `main` or the weekly workflow. Do not run the old convergence scripts on new pilot files or interpret old MARS outputs as this dataset.

### Optional native correlation acceleration

The correlation stage supports a batched C++ MEX kernel using MATLAB's own LAPACK/BLAS. It evaluates the same matrix-log inverse and observed-subspace Gaussian likelihood, with the same convergence tolerance and positive-definiteness checks. Within the fixed-point solve, it computes only the diagonal of the matrix exponential; it builds the full correlation matrix once after convergence. Standardized residuals are reused across correlation proposals. The default four native threads evaluate separate contiguous blocks of weekly dates, with deterministic summation in date order; this does not run multiple chains. Set `cfg.correlation_threads=1` to use a single native thread.

```matlab
build_info = build_dsc_mex();              % compile locally for this MATLAB/platform
checks = run_weekly_acceptance();          % includes native/reference and resume checks
benchmark = benchmark_dsc_kernel();        % fixed saved state; no MCMC
cfg = weekly_config();
cfg.correlation_backend = 'mex';           % require the compiled implementation
results = run_weekly_research('pilot',cfg); % fresh run, same 20-sweep/30-minute cap
verify_weekly_pilot(results.pilot.run_dir);
verify_dsc_likelihood(results.pilot.run_dir); % independent final likelihood check
```

The default `auto` backend uses MEX when available and MATLAB otherwise; `matlab` explicitly selects the portable implementation. Native binaries are rebuilt locally and excluded from Git. Checkpoints record the resolved backend, MATLAB version, and a hash of sampler sources and the native binary. A changed implementation requires a fresh run; historical checkpoints are preserved. Native and MATLAB calculations agree within numerical tolerance, but floating-point differences can eventually change individual MCMC accept/reject decisions. Invalid correlation proposals are counted in run diagnostics.

The 13 September rerun completed all 20 warm-up sweeps in 10.83 minutes. Comparing the same first 13 sweeps, mean time fell from 132.58 to 30.69 seconds (**4.32x faster**), with matching proposal counts. Across all 20 new sweeps the mean was 32.46 seconds; peak MATLAB process RSS was 1.149 GiB. All 33 tests, 1,204 final-state correlation matrices and the independent final-likelihood check passed. See [the performance analysis](research/sampler_performance.md) and the run-local `performance_comparison.json` for measurements and projection limits. These warm-up draws do not provide posterior estimates or credible bands.

## Outputs and publication

Derived numerical artifacts live under `outputs/weekly_research/`:

- `data/weekly_panel.mat`, weekly returns/levels and per-cell quality CSVs.
- Descriptive, stationarity, serial-dependence and prior-calibration outputs.
- All 91 full-sample and 52/104-week rolling correlations.
- A run manifest and timestamped chain directories under `runs/`.

Editable LaTeX sources and bibliography are under `research/`. The publication builder reads saved outputs and never reruns estimation. It produces shared vector figures, the complete 91-pair report appendix, and a 12-15-slide Beamer presentation.

```sh
python3 scripts/build_publication.py --help
python3 scripts/build_publication.py --results outputs/weekly_research --pilot outputs/weekly_research/runs/20260913-184612-504-chain1/pilot_summary.json
python3 scripts/check_publication.py
python3 tests/verify_exports.py
```

Final PDFs are `output/pdf/weekly_correlation_report.pdf` and `output/pdf/weekly_correlation_slides.pdf`. Python needs numpy, pandas and matplotlib for publication, plus Pillow for the visual-review contact sheets. TeX needs latexmk, pdflatex and standard article/Beamer packages. The PDF checks use Poppler. MATLAB does all statistical estimation; Python renders saved evidence. `tests/verify_exports.py` independently reconciles every exported level, return and rolling correlation for this supplied snapshot.

The command above rebuilds the accelerated 13 September 2026 snapshot. Its run directory preserves the checkpoint, numerical validation, acceptance tests, kernel benchmark, final-likelihood check and performance comparison. The earlier `20260911-175022-407-chain1` pilot, which completed 13 iterations before its 30-minute guard, remains available for comparison. For a future run, use its own saved pilot-summary path and run-local evidence.

## Skills and agents

Repository-local skills:

- `dsc-data`: prepare/validate data and run descriptive analysis.
- `dsc-matlab`: calibrate, test and pilot the joint estimator.
- `dsc-publication`: rebuild research outputs from saved evidence.

Project agent roles under `.codex/agents/` divide data, MATLAB and publication ownership, with a separate read-only reviewer. Model settings are inherited. Skills and agent roles do not authorize additional computation or external actions.

## Before production

### Requested 68% posterior bands

Bayesian correlation plots should use the posterior median as their center line and the 16th/84th percentiles as pointwise, equal-tailed 68% credible bounds. Compute these across retained draws of **actual** `P(i,j,t)`, not across time, rolling-window estimates, matrix-log coordinates, or state shocks. They are neither simultaneous bands for the entire path nor highest-posterior-density intervals. They condition on the model and the empirical-Bayes prior calibration.

`toolbox/dsc_credible_bands.m` provides the calculation without starting MCMC. Supply one dates-by-pairs-by-draws array per chain, the corresponding original iteration IDs, and each chain's warm-up cutoff. Use matched dates/pair ordering across chains; use date/pair blocks for large outputs. Original iteration IDs prevent double burn-in removal when using already-trimmed sampler chunks. The result preserves retained IDs/counts by chain and is explicitly not publication-ready until convergence has been assessed separately. The current pilot has zero retained posterior draws, so no empirical posterior bands are available and the function rejects that input. Existing rolling-correlation PDFs are unchanged.

```matlab
addpath('toolbox');
% P_by_chain and iteration_ids_by_chain must come from saved posterior chunks.
bands = dsc_credible_bands(P_by_chain, iteration_ids_by_chain, burnin_by_chain);
% bands.lower, bands.median, bands.upper have dimensions dates x pairs.
```

The pilot measures feasibility only. Production needs separate approval and a measured compute budget, four chains, rank-normalized split R-hat below 1.01 and bulk/tail ESS at least 400 for reported quantities. Evaluate serial dependence, a parsimonious AR mean extension, holiday-mask sensitivity and prior sensitivity before treating posterior correlations as final.

## References

- Arias, Rubio-Ramirez and Shin (2023), *Macroeconomic forecasting and variable ordering in multivariate stochastic volatility models*, Journal of Econometrics 235(2), 1054-1086. DOI: [10.1016/j.jeconom.2022.04.013](https://doi.org/10.1016/j.jeconom.2022.04.013).
- Archakov and Hansen (2021), *A New Parametrization of Correlation Matrices*, Econometrica 89(4), 1699-1715. DOI: [10.3982/ECTA16910](https://doi.org/10.3982/ECTA16910).
- Vehtari et al. (2021), *Rank-Normalization, Folding, and Localization: An Improved R-hat for Assessing Convergence of MCMC*. [Paper](https://arxiv.org/abs/1903.08008).
- [JPX closure notice](https://www.jpx.co.jp/english/news/1030/20190115-01.html).
