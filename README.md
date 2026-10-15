# Weekly Dynamic Correlation Model

This repository runs the weekly Bayesian dynamic-correlation model in MATLAB.
The current applied example is Model 1, using NKYTR Index, SPXT Index, SX5T
Index, and TUKXG Index.

## First-Time MEX Setup On This Mac

The default correlation backend is the compiled MEX backend. In MATLAB R2024b
on this Mac, run the following once from the repository root:

    addpath("toolbox","build");
    mex -setup C++
    build_dsc_mex

Choose the configured C++ compiler when MATLAB asks. If MATLAB cannot find a
compiler, install Apple's Command Line Tools in Terminal with
`xcode-select --install`, restart MATLAB, and run the commands again. The
builder detects whether MATLAB is running on Apple Silicon or Intel and creates
the matching files in `toolbox/`:

    dsc_corr_batch_mex.<mex extension>
    dsc_corr_batch_mex_new.<mex extension>
    dsc_mean_mex.<mex extension>

Verify the build with:

    assert(exist("dsc_corr_batch_mex","file") == 3);
    assert(exist("dsc_corr_batch_mex_new","file") == 3);
    assert(exist("dsc_mean_mex","file") == 3);

This compilation is computer-specific. If the repository is copied to another
computer, rebuild the MEX files there rather than copying the compiled files.

## Run Model

From the repository root in MATLAB:

    result = run_weekly_model_report( ...
        WarmupIterations=100, RetainedDraws=50000, ...
        NumChains=4, MaxHours=12, ...
        VariableNames=["NKYTR Index","SPXT Index","SX5T Index","TUKXG Index"], ...
        ParameterSmoothing="draws", CorrelationBackend="mex");

The model reads data/yad_tickers_no_usdcnh_from_20030804.csv, prepares weekly
Friday-as-of returns, estimates the model, saves posterior draws, checks
convergence using correlations P and log-volatility states h, and creates
figures.

## Options

    WarmupIterations       Initial iterations discarded per chain. DEFAULT: 10.
    RetainedDraws          Target saved draws per chain. DEFAULT: 100.
    NumChains              Number of independent chains. DEFAULT: 1.
    ChainID                First chain identifier. DEFAULT: 1.
    MaxHours               Maximum runtime per inference stage. DEFAULT: 12.
                            Prior mode shares one cap across all parameter sets.
    Thin                   Save every Thin-th post-warm-up iteration. DEFAULT: 1.
    ChunkSize              Draws per chunk and convergence-check interval. DEFAULT: 1000.
    Seed                   Random-number seed. DEFAULT: 20260914.
    ParallelChains         Run chains concurrently. DEFAULT: false.
    ConvergenceMode        "off" or "report". DEFAULT: "report".
    AutoExtend             Extend sampling until convergence or cap. DEFAULT: false.
    MaxRetainedDraws       Maximum draws per chain when extending. DEFAULT: 2000.
    RhatThreshold          Maximum accepted R-hat. DEFAULT: 1.01.
    MinESS                 Minimum bulk and tail ESS. DEFAULT: 400.
    MaxMCSERatio            Maximum MCSE/posterior-SD ratio. DEFAULT: 0.05.
    ParameterSource        "posterior" or "prior". DEFAULT: "posterior".
    EstimateParameters     Estimate posterior parameters or load saved ones.
                            Applies to ParameterSource="posterior". DEFAULT: true.
    SmoothStates           Run the smoothing/output stage. DEFAULT: true.
                            false estimates/saves parameters only; latent states
                            are still sampled inside parameter estimation.
    PriorParameterDraws    Independent prior parameter draws. DEFAULT: 10.
                            Applies to ParameterSource="prior" only.
    EstimationEndDate      Optional estimation cutoff date. DEFAULT: "".
    SourceFile             Input CSV. DEFAULT: configured Model 1 CSV.
    VariableNames          Tickers to estimate. DEFAULT: "all".
    ParameterSmoothing     "draws" or "mean". DEFAULT: "mean".
                            Applies to posterior parameter smoothing only.
    CorrelationBackend     "mex", "matlab", "mex_new", or "matlab_new".
                            DEFAULT: "mex". The *_new backends use the
                            analytic log-diagonal Newton solver while preserving
                            the existing model and coordinates. It reuses one
                            eigensystem per step and a small Cholesky solve.
    CorrelationThreads     Threads for correlation calculations. DEFAULT: "auto".
                            "auto" requires RunTests=true. If RunTests=false,
                            MATLAB writes a warning and enables RunTests.
                            A fixed integer from 1 to 64 may be supplied instead.
    RunTests               Run acceptance tests before estimation. DEFAULT: true.
    BuildReport            Build the publication report. MODEL 1 DEFAULT: false.
    VerifyReport           Verify the publication report. MODEL 1 DEFAULT: false.
    PythonExecutable       Python used for report building. DEFAULT: "python3".
    KeepNewestParameterFiles
                            Parameter snapshots to retain. DEFAULT: 3.

When `ParameterSource="posterior", EstimateParameters=false`, the model
automatically loads the newest compatible saved parameter file and performs
only smoothing of the latent states. See the three workflows below.
Short test:

    result = run_weekly_model_report( ...
        WarmupIterations=10, RetainedDraws=100, NumChains=1, ...
        MaxHours=1, BuildReport=false, VerifyReport=false, ...
        VariableNames=["NKYTR Index","SPXT Index","SX5T Index","TUKXG Index"]);

## Three inference workflows

The editable [workflow note](research/inference_workflows.tex) and its
[short PDF](research/inference_workflows.pdf) describe the algorithms and
the interpretation of their uncertainty bands. Examples below illustrate
the options; they are not saved estimation results.

### 1. Estimate parameters, then smooth states

    result = run_weekly_model_report( ...
        ParameterSource="posterior", EstimateParameters=true, ...
        SmoothStates=true, ParameterSmoothing="mean", ...
        WarmupIterations=100, RetainedDraws=100, NumChains=4, ...
        ParallelChains=true, MaxHours=12, CorrelationBackend="mex_new");

First the Metropolis-within-Gibbs sampler jointly updates the full latent paths `B`, `h`, `r`
and their innovation variances `V`, `sig2h`, `sig2r`, using returns through
`EstimationEndDate` (the last available week by default). It saves posterior
parameter draws. With the default `ParameterSmoothing="mean"`, a second
sampler freezes these parameters at their posterior means and smooths states
through the last available week. Its bands describe conditional state
uncertainty and exclude parameter uncertainty. Each stage has its own
`MaxHours` limit.

With `ParameterSmoothing="draws"` and identical estimation/smoothing dates,
the code uses the joint estimation output directly; its posterior meaning
requires adequate convergence. When loading saved draws or extending the
history, the existing draw-cycling smoother is a heuristic and its
convergence is not checked as a common-target MCMC chain.

To reuse posterior parameters later, set
`ParameterSource="posterior", EstimateParameters=false, SmoothStates=true`.
This loads the saved prior calibration and parameter snapshot. An appended
history can revise earlier smoothed states; parameters are not re-estimated.

### 2. Estimate and save parameters only

    result = run_weekly_model_report( ...
        ParameterSource="posterior", EstimateParameters=true, ...
        SmoothStates=false, WarmupIterations=100, RetainedDraws=100, ...
        NumChains=4, ParallelChains=true, MaxHours=12, ...
        CorrelationBackend="mex_new");

This performs the same joint parameter-estimation sampler and diagnostics,
and saves the parameter snapshot for later smoothing. It skips the second
conditional smoother and state figures. "Estimation only" still needs
latent-state draws internally to update the parameters.

### 3. Draw parameters from their priors, then smooth

    result = run_weekly_model_report( ...
        ParameterSource="prior", PriorParameterDraws=10, ...
        SmoothStates=true, WarmupIterations=100, RetainedDraws=100, ...
        NumChains=4, ParallelChains=true, MaxHours=12, ...
        CorrelationBackend="mex_new");

This mode skips posterior parameter estimation and does not load a posterior
parameter file. It draws `K=PriorParameterDraws` independent parameter sets
from the calibrated inverse-Wishart/inverse-gamma priors. For each fixed set,
it runs `C=NumChains` separate conditional state chains, each with its own
warm-up and `R=RetainedDraws` target. The example targets 10 parameter sets
and 4,000 conditional state draws (`10*4*100`), not 4,000 independent
parameter draws. Parameters stay fixed throughout each conditional chain;
they are not replaced after every MCMC sweep.

The target is the equally weighted prior mixture
`q(theta,states | y) = prior_EB(theta) * p(states | y,theta)`.
Observed returns inform the smoothed states but do not update the prior
weights of the parameters. This is a prior-sensitivity calculation, neither
the joint posterior nor prior-predictive simulation. The priors remain
empirical Bayes: the first 104 reserved weeks set initial moments, and
104-week rolling windows through the calibration cutoff set evolution
scales. Therefore it is not a strictly data-free analysis.

`MaxHours` is one shared budget for the entire prior-smoothing stage, across
all `K` parameter sets. A time stop with incomplete, unequal component quotas
does not produce an equally weighted pooled summary. Diagnostics compare
chains only within a fixed parameter set; they do not compare different
parameter sets or establish convergence of the outer prior integral.
`AutoExtend=true` and resume are unsupported in this initial prior mode.
Prior mode requires `SmoothStates=true`. `BuildReport=true` is currently
supported only for posterior mode with smoothing; the standalone note above
covers all three workflows. Run completion alone never establishes
convergence.

## New correlation solver and timing comparison

`mex_new` and `matlab_new` use the analytic exponential derivative described
in [Archakov and Hansen's 2026 reference implementation](https://github.com/reinhardhansen/GFT/blob/v1.2.0/src/GFT.jl).
For the small DSC matrices, they assemble the dense Newton system and solve it
by Cholesky rather than matrix-free conjugate gradients. An initial phase uses
the existing AH fixed-point update until the largest log-diagonal residual is
at most 0.25, then switches once to full Newton steps. Exact trace normalization
reuses the same eigensystem. There are no finite-difference Jacobians,
backtracking evaluations, step-size heuristics, or failed-step fallbacks.
Input validation, finite-value checks, the 200-iteration cap, and the final
positive-definiteness check follow the old backend's contracts. This lean
variant does not include the official solver's global-convergence safeguards.

To time the same 100 synthetic 14-variable correlation matrices with both MEX
backends, without starting estimation or creating output directories:

    addpath("scripts");
    timings = benchmark_correlation_backends;

The benchmark includes the masked likelihood and all five MEX outputs, tests
one and four threads, and reports the median of seven alternating timing
trials. These are times for 100 correlation transformations, not 100 complete
MCMC iterations. An optional output filename writes JSON into an existing
directory.

Measured on 2026-10-12 with MATLAB R2024b (Apple Silicon), including masked
likelihoods and all five outputs:

| Threads | Old MEX | Optimized `mex_new` | Old/new speed ratio |
| --- | ---: | ---: | ---: |
| 1 | 18.32 ms | 7.85 ms | 2.33x |
| 4 | 5.24 ms | 2.40 ms | 2.18x |

The maximum correlation difference was below `6e-9`; 15 kernel regression
tests and the MATLAB reference round-trip test passed. Full timing trials and
the prior finite-difference prototype comparison are saved in
`outputs/correlation_backend_benchmark.json`.

## Outputs

All model results go to the repository's outputs/ folder. This location is fixed;
there is no output-directory option. The layout is:

    outputs/
        parameters/               Saved parameter snapshots; newest three kept by default.
        runs/<run-id>/            Draws, provenance, convergence results, and figures.
        latest_model_report.json  Summary of the newest run.
        matlab_test_results.json  Results of the latest acceptance tests.

Parameter snapshots use names like
`parameters_20261006-213950-chain1-w1000-r10000-thin1.mat`.
The filename is the matching `runs/<run-id>/` folder name; each run manifest
also records the parameter-file path.

The data/ folder is created temporarily inside outputs/ during preparation
and removed after the run. In posterior mode, EstimateParameters=false loads
parameters from outputs/parameters/ automatically. Prior mode skips this load
and saves its own parameter draws and component metadata with its run.
After a successful run, only the three newest completed run folders are kept.

Permanent results include:

    chains/                                  Posterior draw chunks.
    run_manifest.json                        Data and code provenance.
    convergence/convergence_diagnostics.json Convergence results.
    convergence/chain_mean_paths.pdf         Chain means for P and h entries.
    figures/                                 Bayesian correlation-path PDFs.

Temporary panel-analysis files and regenerable run handoffs are always removed
after the run. There is no option to retain them.

BuildReport=true requires Python with NumPy, pandas, and Matplotlib, plus a
TeX installation with latexmk. VerifyReport=true also requires Poppler and
Pillow. Generated .tex files, copied charts, LaTeX build files, and verification
images are temporary and deleted automatically. The three newest reports and
three newest figures are kept directly in report/, with no subfolder:

    report/report_<run-id>.pdf
    report/figure_<run-id>.pdf

The run ID matches the corresponding folder under outputs/runs/.

`run_weekly_model_report` is the single entry point. Select the variables with
`VariableNames`; Model 1 is the four-ticker example above.
