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
the two matching files in `toolbox/`:

    dsc_corr_batch_mex.<mex extension>
    dsc_mean_mex.<mex extension>

Verify the build with:

    assert(exist("dsc_corr_batch_mex","file") == 3);
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
    MaxHours               Maximum runtime. DEFAULT: 12.
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
    EstimateParameters     Estimate parameters, then smooth P and h. DEFAULT: true.
    EstimationEndDate      Optional estimation cutoff date. DEFAULT: "".
    SourceFile             Input CSV. DEFAULT: configured Model 1 CSV.
    VariableNames          Tickers to estimate. DEFAULT: "all".
    ParameterSmoothing     "draws" or "mean". MODEL 1 DEFAULT: "draws".
    CorrelationBackend     "mex" or "matlab". DEFAULT: "mex".
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

When EstimateParameters=false, the model automatically loads the newest
compatible saved parameter file and performs only smoothing of P and h.
Short test:

    result = run_weekly_model_report( ...
        WarmupIterations=10, RetainedDraws=100, NumChains=1, ...
        MaxHours=1, BuildReport=false, VerifyReport=false, ...
        VariableNames=["NKYTR Index","SPXT Index","SX5T Index","TUKXG Index"]);

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
and removed after the run. EstimateParameters=false loads parameters from
outputs/parameters/ automatically.
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
