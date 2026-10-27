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

The default configuration selects the four series and constant means. Run
estimation through the production entry point, `run_weekly_model_report`:

```matlab
result=run_weekly_model_report( ...
    NumChains=4,PerDraw=10000,WarmupIterations=100,MaxHours=12);
```

This invocation requests a fixed draw target and an overall time cap. The
production function estimates parameters and saves joint draws and parameter
summaries.
It uses Friday-as-of weekly log returns from
`data/yad_tickers_no_usdcnh_from_20030804.csv`. The first 104 returns calibrate
initial `mu`, `h` and `r` levels and are excluded from the likelihood. These
weeks supply 103 adjacent-week pairs for AR(1) analysis and are distinct from
the posterior period. Evolution variances use coordinate-specific rolling
empirical-Bayes calibration through the estimation cutoff.

The public options are:

| Option | Default | Meaning |
| --- | --- | --- |
| `WarmupIterations` | 10 | Discarded iterations per chain |
| `NumChains` | 1 | Independent chains |
| `PerDraw` | 100 | Saved draws per chain after warm-up and thinning |
| `MaxHours` | 12 | Overall run time cap |
| `Thin` | 1 | Save every selected post-warm-up iteration |
| `ChunkSize` | `PerDraw` | Draws per saved chunk; specify only to override the default |
| `Seed` | 20260914 | Base deterministic seed; chain IDs are fixed at `1` through `NumChains` |
| `CorrelationBackend` | `"mex_new"` | Analytic inverse-map solver |
| `CorrelationThreads` | `"auto"` | Native thread count selected within the CPU allocation; integer 1–64 overrides |
| `MeanBackend` | `"auto"` | Use the compiled C++ (`"native"`) constant-mean Gaussian draw when available, otherwise the MATLAB implementation; this selects an implementation, not a different mean model |
| `CreatePDF` | `true` | Generate a pooled mean correlation-path PDF with 16th/84th percentile 68% pointwise credible bands; with multiple chains also generate a per-chain mean overlay for visual convergence diagnosis |
| `OutputDirectory` | Existing `outputs/` | Destination for prefixed files |
| `EstimateParameters` | `true` | Estimate `mu`, `sig2h` and `sig2r` jointly with states |
| `ParallelChains` | `false` | Run chains simultaneously using an existing process pool |
| `SourceFile` | Repository CSV | Market-level input file |
| `VariableNames` | Four tickers listed above | Series to include, in model order |

`ParallelChains=true` requires an existing MATLAB process pool. If the
Parallel Computing Toolbox is not installed or licensed, the run warns and
runs the same chains sequentially; draws are unchanged, only slower. The
manifest then records `ParallelChains=false`. Sequential
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
selected count is fixed across chains and batches and saved in the parameter
MAT report. An explicit integer such as `CorrelationThreads=2` bypasses
timing. Other backends or runs without native context use the CPU allocation
limit directly. Remote pools retain a limit of one thread per chain because
local core count does not describe their CPU allocation; use an explicit
override there. If core detection is unavailable, the default limit is one.
These threads apply to the MEX correlation backends; MATLAB backends retain
MATLAB's own thread management.

With `EstimateParameters=false`, the loader automatically selects the latest
compatible constant-mean parameter snapshot; historical dynamic-mean snapshots
are a different model. Fresh estimates use every available week after the
initial calibration period; there is no public estimation-cutoff option.

Estimation returns the joint parameter/state draws directly.
Separate loaded-parameter or extended-history smoothing fixes
`sig2h` and `sig2r` at their saved posterior means and resamples `mu` and the
state paths. There is no public smoothing-method option.

## Outputs and checks

The workflow writes files inside the existing `outputs/` directory, using an
automatically generated `constant_mu_<timestamp>_` prefix. The timestamp uses
Europe/Madrid time in `yyyyMMddTHHmmssSSS` format; callers cannot override it.
After a successful production run, each type of information has one home:

| File | Contents |
| --- | --- |
| `<prefix>run_manifest.json` | Public invocation options, model equations, three prior-setting notes, estimation/calibration dates, AR(1) ticker/p-value rows, and timings |
| `<prefix>parameters.mat` | Posterior means and standard deviations for `mu`, `sig2h`, `sig2r`, their prior/posterior comparison table, and the calibrated priors and internal context needed to reuse the run |
| `<prefix>chain_001_posterior_chunk_000001.mat` | Complete retained joint draws: `mu`, `sig2h`, `sig2r`, the full `h` paths and packed correlation paths `P_pairs` |
| `<prefix>correlation_paths.pdf` | Generated by default; pooled mean correlations across all retained draws, weighted by actual retained counts, with 16th/84th percentile 68% pointwise credible bands |
| `<prefix>convergence_correlation_paths.pdf` | Generated by default only for multiple chains; separate mean paths by original chain with dynamic colors and a legend, for visual diagnosis |

Each chunk contains matching parameter and state draws, identified by chain
and iteration. `P_pairs` stores all distinct off-diagonal correlations at every
date and draw; the full symmetric `P` is recovered by reflecting those entries
and setting the diagonal to one. There is no need to store a second dense copy.
Warm-up is excluded from these retained draws, and thinning follows `Thin`.

Raw parameter draws live only in the chunks. The parameter report references
them, so keep each run's MAT files together when reusing saved
parameters. Loading reconstructs the parameter draws and checks the supplied
calibration/training sample against saved fingerprints. Source priors are read
from `parameters.mat`. A loaded-parameter smoothing run writes its own summary,
marked as conditional smoothing rather than a new parameter estimate.

The separate `prior_summary.json`, `calibrated_priors.mat` and
`inference_metadata.json` are temporary staging files. Their required numerical
context is kept in `parameters.mat`, while the manifest contains only the
compact description above. There is no final
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
For multi-chain runs, the convergence PDF appends pages with prior and
posterior means and standard deviations, followed by weekly AR(1) two-sided
p-values for `H0: rho = 0`. Each AR(1) result row records its ticker and
p-value. The regression uses the initial 104
weekly returns reserved for prior calibration, yielding 103 consecutive-week
pairs when fully observed. This period is separate from the posterior
likelihood sample. Parameter summaries retain their estimation or
conditional-smoothing scope, shown with their dates. The ordinary pooled
correlation-path PDF remains unchanged. Each correlation panel has 0.1 padding
beyond the displayed extrema. The pooled limits include the mean and both 68%
band edges, while convergence limits include every original-chain mean.
The correlation color scale and colorbar remain fixed to `[-1, 1]`.
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

In the convergence PDF, posterior means and standard deviations include the
signed percent change from their corresponding prior moments in parentheses,
computed as `100*(posterior/prior - 1)` from unrounded stored values; undefined
changes are shown as `(N/A)`.

Prior moments are analytic: `mu` has mean `mu0` and marginal standard
deviations `sqrt(diag(mu_cov0))`. For an inverse-gamma variance prior with
shape `a` and scale `b`, the mean is `b/(a-1)` and the standard deviation is
`b/((a-1)*sqrt(a-2))` when `a>2`. For `1<a<=2`, the standard deviation is
infinite and is saved as `Inf`. These are the actual empirical-Bayes priors
calibrated for the estimation sample, so they incorporate calibration data.
Future estimation runs use independent coordinate-specific inverse-gamma
priors for the innovation variances. Calibration uses the existing 104-week
rolling covariance and correlation histories over the full sample through the
estimation cutoff, with the existing 5% covariance shrinkage. For each
log-variance coordinate, the raw scale `v_hj` is the sample variance of
`diff(log rolling variance_j)`; for each matrix-log correlation coordinate,
`v_rk` is the sample variance of `diff(matrix-log correlation coordinate_k)`.
Each coordinate is calibrated separately; there is no pooling or median
across series or pairs. The initial 104-week level calibration remains
unchanged.

The fixed default shape is `a=4`, and each coordinate's prior mean is
`2*max(v_coordinate,1e-8)`, with inverse-gamma scale
`b_coordinate=(a-1)*mean`. Thus the prior standard deviation is
`mean/sqrt(a-2)=mean/sqrt(2)`. Relative to shape `a=10` centered on that same
coordinate's raw `v`, the new mean is twice as large and its standard
deviation is four times as large. No pooling option is exposed. Each run has
four `h` scales and six `r` scales, with the latter in strict-lower-triangle
pair order. The parameter MAT report records the coordinate-specific prior vectors;
prior/posterior tables and PDF summaries use each coordinate's own moments.
Older saved scalar priors remain loadable using their original values
broadcast across coordinates.

Manifest prior notes record the actual innovation-mean multiplier (2), inverse-
gamma shape (4), and Normal `mu` standard-deviation multiplier (2). The `mu`
multiplier of 2 corresponds to a covariance factor of 4.
No prior simulation is required for this comparison.

The manifest's `config` section records the public PascalCase options accepted
by `run_weekly_model_report`, with `CorrelationThreads="auto"` retained as the
requested option. `ChainID` and `resume` are not public options: chain IDs are
fixed at `1` through `NumChains`, and every invocation starts a fresh run. The
manifest does not expose internal `weekly_config` fields or derived
values such as the sweep limit and time limit. Removed metadata fields are
excluded whenever a manifest is written, including intermediate snapshots.
The selected backends, resolved thread count, actual draw counts, outcome and
calibrated numerical priors remain in the internal MAT report. The mean model is fixed to
constant means; the weekly-return AR(1) results are a separate diagnostic.
The estimation sample end is recorded once as
`inference.parameter_estimation_end` rather than repeated in `config`.

The default correlation solver reuses its eigensystem and an analytic Newton
system. Correlation warm starts are always enabled on a fresh run and on
continued numerical updates. Each inverse-map solve uses the last solver
solution as its initial guess; this is separate from MCMC warm-up and does not
resume a stopped run. The inverse map still uses the same tolerance,
200-iteration cap, finite-value checks and positive-definiteness checks.
Masked likelihoods use the observed subspace without deleting dates.

The repository layout is compact: `scripts/` contains the basic launcher and
the current saved-run correlation-path PDF generator; `toolbox/` contains the
estimator, numerical kernels, reporting helpers, and standalone
`analyze_weekly_panel.m` data-quality helper. `build/` contains only
`build_dsc_mex.m`; `data/` contains model inputs; and `outputs/` stores flat,
timestamp-prefixed run groups. The obsolete `research/` tree and Python
publication workflow have been removed.

After successful final publication, including requested PDFs, retain exactly
the three newest completed run groups by generated timestamp, counting the
newly completed run among the three. Always protect active runs, stopped or
time-limited runs, and source parameter groups needed by protected runs. Do
not prune at startup or after failed publication. Inactive failed, incomplete,
and unclassified prefixed groups are eligible for cleanup; delete all prefixed
artifacts in an eligible group together, including MAT, JSON, PDF and temporary
files. Do not create backups, inventory sidecars, or folders. Internal
`<prefix>run_state.json` records process ownership while a run is active or
stopped; it is removed after successful completed publication. Public
manifests contain no `lifecycle` section. Protect a run with an active process
or an unverifiable owner; a stale inactive owner is not automatically
protected. Existing runs are registered before applying this policy.
Retention is fixed at three; there is
no public `KeepNewestParameterFiles` option. Routine cleanup must not create
pre-cleanup backup archives or inventory JSON sidecars; retain only the
documented run outputs.
