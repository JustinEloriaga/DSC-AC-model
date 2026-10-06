# Weekly Dynamic Correlation Model

This repository runs the weekly Bayesian dynamic-correlation model in MATLAB.
The current applied example is Model 1, using NKYTR Index, SPXT Index, SX5T
Index, and TUKXG Index.

## Run Model 1

From the repository root in MATLAB:

    result = run_weekly_model1_report( ...
        WarmupIterations=100, RetainedDraws=50000, ...
        NumChains=4, MaxHours=12, ...
        ParameterSmoothing="draws", CorrelationBackend="mex");

The model reads data/yad_tickers_no_usdcnh_from_20030804.csv, prepares weekly
Friday-as-of returns, estimates the model, saves posterior draws, checks
convergence using correlations P and log-volatility states h, and creates
figures.

## Options

    WarmupIterations       Initial iterations discarded per chain.
    RetainedDraws          Target saved draws per chain.
    NumChains              Number of independent chains.
    MaxHours               Maximum runtime.
    Thin                   Save every Thin-th post-warm-up iteration.
    ChunkSize              Draws stored in each posterior chunk file.
    ParameterSmoothing     "draws" or "mean".
    CorrelationBackend     "mex" or "matlab".
    CorrelationThreads     Threads used by correlation calculations.
    RunTests               Run acceptance tests before estimation.
    BuildReport            Build the publication report.
    VerifyReport           Verify the publication report.
    KeepNewestParameterFiles
                            Number of parameter snapshots to retain; default 3.
    KeepIntermediateArtifacts
                            Default false; remove regenerable handoff files.

Short test:

    result = run_weekly_model1_report( ...
        WarmupIterations=10, RetainedDraws=100, NumChains=1, ...
        MaxHours=1, BuildReport=false, VerifyReport=false);

## Outputs

The default output directory is outputs/weekly_research_model1/.
Each run is stored under outputs/weekly_research_model1/runs/<run-id>/.

Permanent results include:

    chains/                                  Posterior draw chunks.
    result.mat                               Aggregated run result.
    run_manifest.json                        Data and code provenance.
    convergence/convergence_diagnostics.json Convergence results.
    convergence/chain_mean_paths.pdf         Chain means for P and h entries.
    figures/                                 Bayesian correlation-path PDFs.
    latest_model_report.json                 Summary of the newest run.

Temporary panel-analysis files and regenerable run handoffs are removed after
the run by default. Set KeepIntermediateArtifacts=true only for debugging.

## Data And Parameters

The raw CSV is never modified. Levels are converted to weekly percent
log-returns using each market's latest observation on or before Friday.

Saved parameter files are kept under outputs/weekly_research_model1/parameters/.
The newest three are retained by default. To load an existing parameter file:

    result = run_weekly_model1_report( ...
        EstimateParameters=false, ...
        ParameterFile="/absolute/path/to/dsc_parameters_file.mat");

## Other Entry Points

    run_weekly_research('analysis')  Prepare data and priors only.
    run_weekly_acceptance()           Run MATLAB acceptance tests.

run_weekly_model1_report is the recommended entry point.
