# Weekly Dynamic Correlation Model

The current MATLAB model uses four weekly equity return series: **NKYTR Index,
SPXT Index, SX5T Index, and TUKXG Index**. Each series has one estimated constant
mean. Volatility and the joint correlation matrix remain stochastic.

The observation equation is

\[
y_t=\mu+D_tP_t^{1/2}z_t,\qquad z_t\overset{\mathrm{iid}}{\sim}N(0,I),
\qquad D_t=\operatorname{diag}(\exp(h_t/2)).
\]

The independent shocks are i.i.d.; returns have time-varying conditional
covariance. The random walks for log variances and matrix-log correlation
coordinates have estimated innovation **variances** `sig2h` (four coordinates)
and `sig2r` (six coordinates). The evolving states you call `delta` and `gamma`
are named `h` and `r` in the code; they are distinct from these variances.
The Archakov–Hansen inverse map produces a coherent positive-definite `P_t`.
The sampler estimates four entries of `mu`, without sampling a mean path
`B_t` or its innovation covariance `V`.

## Build

From the repository root in MATLAB:

```matlab
addpath("toolbox","build","scripts");
build_dsc_mex
assert(exist("dsc_corr_batch_mex_new","file")==3);
assert(exist("dsc_constant_mean_mex","file")==3);
```

On macOS the builder invokes Apple Command Line Tools directly and writes
artifacts into the existing `toolbox/` directory. If the tools are missing,
install them with `xcode-select --install`. On other platforms, configure a
C++ compiler with `mex -setup C++` first. Rebuild binaries on each computer.

The build also includes the classical `dsc_corr_batch_mex` for numerical
regressions and timing comparisons. The historical dynamic-mean FFBS kernel
is not required or built by the current workflow.

## Run

The default configuration selects the four series and constant means. A
bounded pilot performs at most 20 iterations or 30 minutes, with every
iteration treated as warm-up:

```matlab
cfg=weekly_config();
cfg.output_root=fullfile(pwd,"outputs");
pilot=run_weekly_research("pilot",cfg);
```

`run_weekly_research("analysis",cfg)` prepares descriptive results and priors
without starting MCMC. Both modes require existing output directories.

For an explicitly budgeted estimation run, the production entry point retains
its existing name, `run_weekly_model_report`:

```matlab
result=run_weekly_model_report( ...
    WarmupIterations=100,RetainedDraws=1000, ...
    NumChains=4,MaxHours=12);
```

This example requests a fixed draw target and an overall time cap. The normal
research pilot remains capped at 20 warm-up iterations or 30 minutes. The
production function estimates parameters and saves joint draws and parameter
summaries.
It uses Friday-as-of weekly log returns from
`data/yad_tickers_no_usdcnh_from_20030804.csv`. The first 104 returns calibrate
initial priors and are excluded from the likelihood. Evolution priors retain
the documented rolling empirical-Bayes calibration.

Useful options are:

| Option | Default | Meaning |
| --- | --- | --- |
| `WarmupIterations` | 10 | Discarded iterations per chain |
| `RetainedDraws` | 100 | Saved draws per chain after warm-up and thinning |
| `NumChains` | 1 | Independent chains |
| `MaxHours` | 12 | Overall run time cap |
| `Thin` | 1 | Save every selected post-warm-up iteration |
| `ChunkSize` | `RetainedDraws` | Draws per saved chunk; specify only to override the default |
| `Seed` | 20260914 | Base deterministic seed; chain identity is retained |
| `CorrelationBackend` | `"mex_new"` | Analytic inverse-map solver |
| `CorrelationThreads` | `"auto"` | Native thread count selected within the CPU allocation; integer 1–64 overrides |
| `MeanBackend` | `"auto"` | Native constant-mean Gaussian draw when available |
| `CreatePDF` | `true` | Generate a pooled mean correlation-path PDF with 16th/84th percentile 68% pointwise credible bands; with multiple chains also generate a per-chain mean overlay for visual convergence diagnosis |
| `OutputDirectory` | Existing `outputs/` | Destination for prefixed files |
| `EstimateParameters` | `true` | Estimate `mu`, `sig2h` and `sig2r` jointly with states |
| `BuildReport`, `VerifyReport` | `false` | Historical report options; unsupported for the current model |

`ParallelChains=true` requires an existing MATLAB process pool. Sequential
chains avoid creating MATLAB job directories. Retained targets are per chain;
checkpoints retain the original chain identifiers, seeds and RNG states.
Each invocation starts a fresh run and rejects an existing output prefix.
Checkpoints are written at initialization, after each full retained-chunk
boundary, and on every return after flushing any short buffer, including time
limits and caught failures. An unexpected MATLAB or process termination can
lose work since the last chunk boundary. Public resume remains unavailable.

`CorrelationThreads` is automatic unless explicitly set. Let `C` be the CPU
core count reported by MATLAB and `K` the number of simultaneous chains:
`K=NumChains` with `ParallelChains=true`, otherwise `K=1`. The per-chain upper
limit is `max(1,min(64,floor(C/K)))`. For example, an
8-core computer gives a limit of 8 for sequential chains or 2 per chain for
four parallel chains; automatic selection may choose fewer. With the default
`mex_new` backend, after data and priors are prepared, a brief deterministic
native benchmark compares candidate counts within this limit using a
representative observed-mask workload. It uses no RNG or MCMC draws. The
selected count is fixed across chains and batches and saved in the run
configuration. An explicit integer such as `CorrelationThreads=2` bypasses
timing. Other backends or runs without native context use the CPU allocation
limit directly. Remote pools retain a limit of one thread per chain because
local core count does not describe their CPU allocation; use an explicit
override there. If core detection is unavailable, the default limit is one.
The single-chain research pilot uses the sequential allocation. These threads
apply to the MEX correlation backends; MATLAB backends retain MATLAB's own
thread management.

With `EstimateParameters=false`, the loader selects a compatible constant-mean
parameter snapshot. `ParameterFile` can select one explicitly; historical
dynamic-mean snapshots are a different model. `EstimationEndDate` optionally
limits the training sample when estimating new parameters.

With no earlier cutoff, estimation returns the joint parameter/state draws
directly. Separate loaded-parameter or extended-history smoothing fixes
`sig2h` and `sig2r` at their saved posterior means and resamples `mu` and the
state paths. There is no public smoothing-method option.

## Outputs and checks

The workflow writes files inside the existing `outputs/` directory, using an
automatically generated `constant_mu_<timestamp>_` prefix. The timestamp uses
Europe/Madrid time in `yyyyMMddTHHmmssSSS` format; callers cannot override it.
After a successful production run, each type of information has one home:

| File | Contents |
| --- | --- |
| `<prefix>run_manifest.json` | Model, configuration, data source, priors and calibration, inference metadata, run outcome and timings |
| `<prefix>parameters.mat` | Posterior means and standard deviations for `mu`, `sig2h`, `sig2r`, and their prior/posterior comparison table |
| `<prefix>chain_001_posterior_chunk_000001.mat` | Complete retained joint draws: `mu`, `sig2h`, `sig2r`, the full `h` paths and packed correlation paths `P_pairs` |
| `<prefix>correlation_paths.pdf` | Generated by default; pooled mean correlations across all retained draws, weighted by actual retained counts, with 16th/84th percentile 68% pointwise credible bands |
| `<prefix>convergence_correlation_paths.pdf` | Generated by default only for multiple chains; separate mean paths by original chain with dynamic colors and a legend, for visual diagnosis |

Each chunk contains matching parameter and state draws, identified by chain
and iteration. `P_pairs` stores all distinct off-diagonal correlations at every
date and draw; the full symmetric `P` is recovered by reflecting those entries
and setting the diagonal to one. There is no need to store a second dense copy.
Warm-up is excluded from these retained draws, and thinning follows `Thin`.

Raw parameter draws live only in the chunks. The parameter report references
them and the manifest, so keep each run's files together when reusing saved
parameters. Loading reconstructs the parameter draws and checks the supplied
calibration/training sample against saved fingerprints. Source priors are read
from the manifest. A loaded-parameter smoothing run writes its own summary,
marked as conditional smoothing rather than a new parameter estimate.

The separate `prior_summary.json`, `calibrated_priors.mat` and
`inference_metadata.json` are consolidated into the manifest. There is no final
`result.mat`, `posterior_correlation_bands.mat`, per-chain summary or checkpoint
file. Intermediate files remain if the workflow fails before successful
finalization.

The production run generates the pooled correlation-path PDF by default. Its
mean uses all retained draws, weighted by each chain's actual retained count;
the band is the pointwise 16th and 84th percentile, not mean plus or minus a
standard deviation. For multiple chains, it also creates a separate
`convergence_correlation_paths.pdf` with one mean line per original chain and
a dynamic-color legend to compare chain means visually. Each PDF can contain
multiple pages when needed. A one-chain run
produces only the pooled PDF. Set `CreatePDF=false` to skip PDF generation.
The standalone `plot_bayes_correlation_paths(output_directory,prefix)`
function remains available for figures from saved chunks and the manifest; it
does not create a bands MAT file.

Use this call to skip PDF generation:

```matlab
run_weekly_model_report(CreatePDF=false)
```

New estimation runs save a prior-versus-posterior comparison table as
`parameters.prior_posterior_summary` in `<prefix>parameters.mat`. It has one row per
`mu` and `sig2h` series and per `sig2r` pair, with the parameter name, series
or pair label, prior distribution, prior mean, prior standard deviation,
posterior mean and posterior standard deviation. For example:

```matlab
saved=load("outputs/<prefix>parameters.mat","parameters");
disp(saved.parameters.prior_posterior_summary)
```

Posterior means and standard deviations use all retained parameter draws
pooled across the original chains, after warm-up and thinning. Standard
deviations use the sample convention (`N-1` denominator) and describe posterior
spread. The comparison retains actual draw counts in
`parameter_summary_metadata`.

Prior moments are analytic: `mu` has mean `mu0` and marginal standard
deviations `sqrt(diag(mu_cov0))`. For an inverse-gamma variance prior with
shape `a` and scale `b`, the mean is `b/(a-1)` and the standard deviation is
`b/((a-1)*sqrt(a-2))` when `a>2`. For `1<a<=2`, the standard deviation is
infinite and is saved as `Inf`. These are the actual empirical-Bayes priors
calibrated for the estimation sample, so they incorporate calibration data.
No prior simulation is required for this comparison.

The default correlation solver reuses its eigensystem and an analytic Newton
system. Proposal warm starts change only initialization: the inverse map
still uses the same tolerance, 200-iteration cap, finite-value checks and
positive-definiteness checks. Masked likelihoods use the observed subspace
without deleting dates.

Output retention protects only active runs and the explicitly selected stopped
run. Incomplete, failed and unclassified runs are eligible for cleanup. The
newest-three rule applies to complete prefixed run groups and removes all
artifacts belonging to older eligible prefixes together.

Historical 14-series fits, dynamic-mean code, report sources and prior timing
results remain reference material. They do not describe the current production
model.
