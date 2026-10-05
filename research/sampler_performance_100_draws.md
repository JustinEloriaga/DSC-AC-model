# 100-Draw Sampler Performance

Measured on 5 October 2026 with MATLAB R2024b Update 1 on this Mac.

## Result

| Worktree | Before (seconds) | After (seconds) | Saved (seconds) | Reduction | Speedup |
|---|---:|---:|---:|---:|---:|
| Model 1, `jfrrnotall-save` | 252.836 | 33.083 | 219.753 | 86.92% | 7.64x |
| Shared sampler, `jfrr` | 204.779 | 32.182 | 172.597 | 84.28% | 6.36x |

Each measurement is **100 retained draws after 10 warm-up sweeps** on the
same saved Model 1 input: 1,100 weekly dates, NKYTR, SPXT, SX5T and TUKXG,
seed 20260914, chain 1, thinning 1 and chunks of 10. Estimation jointly
samples parameters and historical states (`ParameterSmoothing="draws"`).
Original source trees were frozen before editing.

The timed interval includes sampling, every completed-sweep checkpoint,
all retained chunks and final sampler files. Warm-up, MATLAB startup,
data analysis, plotting and report generation are excluded. Both versions
retain the same number of draws and use the same panel, masks and priors.
These are measured times from separate runs, not projections; machine load
can change exact timings.

## Changes

The largest improvement comes from compiling and using the existing
batched correlation MEX kernel, which the preceding Model 1 run had not
used. Additional changes accelerate the remaining work:

- `dsc_mean_mex.cpp` implements the existing Kalman filter and backward
  state sampler with MATLAB's LAPACK. MATLAB supplies all random noise.
- `dsc_volatility_likelihood.m` evaluates the exact Gaussian likelihood
  change for one volatility coordinate from cached whitened residuals.
  It avoids solving every date's covariance for each slice proposal.
- Retained correlations are packed by indexing the matrix pages directly.
- Sampler signatures include the new helper and mean binary, preventing
  a resume with a different implementation.
- Model 1's dispersed chain initialization now uses explicit column
  parameters for gamma draws, supporting the Dynare function that shadows
  MATLAB's gamma function on this installation.

The statistical model, priors, likelihood tolerances, slice updates,
warm-up, thinning and completed-sweep checkpoint policy are preserved.

Model 1's measured stage costs across the 100 retained draws:

| Stage | Before (seconds) | After (seconds) |
|---|---:|---:|
| Historical means | 44.874 | 0.332 |
| Correlations | 201.542 | 29.917 |
| Volatility | 3.935 | 0.364 |
| Checkpoints | 1.808 | 1.915 |

The final profile is dominated by the compiled correlation kernel. Tests
with 8 and 16 threads yielded 32.0 and 29.6 seconds per 100 draws, versus
30.1 seconds in an earlier four-thread run. These small differences do
not justify changing the four-thread default.

## Verification

The final Model 1 comparison checks every retained chunk, including
`P_pairs`, `B`, `h`, `V`, `sig2h`, `sig2r` and their metadata.
All 100 retained draws were verified. RNG states and all proposal counts
match exactly. Maximum absolute differences were:

| Quantity | Difference |
|---|---:|
| Correlation paths | 4.29646e-9 |
| Historical means | 7.07631e-11 |
| Mean evolution covariance | 6.93127e-18 |
| Volatility paths and variance parameters | 0 |

The sampler, native kernels, serial chain coordinator and inference tests
report **52 passed, 0 failed and 1 optional test skipped**. The separate
14-variable CSV-to-estimation-to-smoothing-to-reload integration test also
passed. Current source and binary signatures match both final timed runs.
This performance benchmark does not establish posterior convergence.

## Running

From the Model 1 worktree, compile once for this MATLAB installation:

```matlab
build_dsc_mex();
result = run_weekly_model1_report( ...
    WarmupIterations=10, RetainedDraws=100, ...
    CorrelationBackend="mex", CorrelationThreads=4, RunTests=false);
```

Both kernels are already built locally. `CorrelationBackend="auto"` also
uses them. `"matlab"` selects the interpreted reference mean/correlation
implementations. Everything runs through MATLAB; Python is not needed.
Native binaries are excluded from Git and must be rebuilt on another
device. Historical checkpoints remain available but cannot be resumed
with the changed sampler implementation.

`benchmark_dsc_sampler` measures a sampler source tree on input from a
saved checkpoint, with optional profiling. `compare_dsc_sampler_benchmarks`
reconciles complete retained output and writes the timing comparison.

Raw reports, frozen baseline sources, profiles, tests, checkpoints and
draws remain local in the primary `DSC-AC-model` worktree under
`outputs/weekly_research/performance/sampler_optimization/`. Final evidence:

- `model1-optimized-final-100-20261005-214359-090/performance_comparison.json`
- `main-optimized-final-100-20261005-214435-010/performance_comparison.json`
- `model1_test_results.mat`

