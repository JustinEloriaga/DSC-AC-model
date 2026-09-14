# MATLAB sampler speed review and rerun

The optimized pilot completed 20 sweeps in 10.83 minutes. Its mean sweep cost was 32.46 seconds. On the first 13 completed sweeps, the old implementation averaged 132.58 seconds and the new implementation averaged 30.69 seconds: **4.32x faster**. These full-sweep measurements were taken in separate sessions on the same machine; the fixed-input microbenchmark below compares both kernels in one session.

## What changed

- A batched C++ MEX kernel calls MATLAB's LAPACK/BLAS for correlation reconstruction and likelihood evaluation.
- Fixed-point iterations calculate only the diagonal of the matrix exponential, and the full matrix is constructed after convergence.
- 4 native threads evaluate independent date blocks; likelihood terms are summed in date order. The sequential MCMC coordinate updates and random draws are unchanged.
- Standardized return innovations are calculated once per correlation sweep.
- Resolved backend, implementation hashes, invalid proposals and per-sweep timings are recorded. Checkpoints reject a changed implementation.

## Evidence

The fixed saved-state kernel benchmark improved from 0.3482 to 0.0807 seconds (**4.31x**), with absolute log-likelihood difference 1.08e-08. Common-sweep correlation and volatility proposal counts match the historical run: **True**. Matching counts are a reproducibility check, not proof of identical floating-point states.

The rerun retained all 1204 weeks, 14 variables and 91 pairs. All 1204 final-state correlation matrices passed symmetry, unit-diagonal, positive-definiteness and matrix-log roundtrip checks; maximum coordinate error was 1.58e-08. Invalid completed-sweep correlation proposals: 0; invalid volatility proposals: 0.

The final saved likelihood also agrees with the independent MATLAB calculation: absolute difference 1.33e-08, against a required tolerance of 1e-05. This validation consumes no random draws and does not advance the chain.

Peak MATLAB process resident memory: 1.149 GiB (macOS `time -l`, including preparation and validation).

| Stage | Historical total (13 sweeps), seconds | New total (20 sweeps), seconds |
|---|---:|---:|
| Mean | 10.38 | 3.16 |
| Mean variance | 0.05 | 0.05 |
| Correlation | 1711.12 | 642.77 |
| Correlation variance | 0.02 | 0.03 |
| Volatility | 1.99 | 3.09 |
| Volatility variance | 0.00 | 0.01 |
| Checkpoint | 0.50 | 0.70 |

## Production implications

At the measured new mean, 4,000 sweeps would take approximately **36.1 hours per chain**, or 144.3 hours for four serial chains. These are linear compute projections from warm-up: mean sweep timings exclude checkpoint writes, and this pilot does not measure production posterior-chunk storage. Later sampling and simultaneous chains can have different costs. Native date threads already share the CPU, so concurrent chains are not assumed to gain another fourfold speedup.

All 20 sweeps are warm-up. There are **0 retained posterior draws**, and convergence has not been established. The rerun therefore does not yet supply Bayesian correlation estimates or 68% posterior bands. Production still requires its own run budget and the agreed convergence and sensitivity checks.

## Saved artifacts

- New run: `outputs/weekly_research/runs/20260913-184612-504-chain1`
- Historical run: `outputs/weekly_research/runs/20260911-175022-407-chain1`
- Kernel benchmark: `kernel_benchmark.json` in the new run (copied from `outputs/weekly_research/performance/kernel_benchmark.json`).
- Acceptance results: `acceptance_tests.json` in the new run (33 passed).
- Numerical validation: `validation.json` in the new run.

Rebuild this summary without rerunning estimation:

```sh
python3 scripts/summarize_speed_pilot.py --run outputs/weekly_research/runs/20260913-184612-504-chain1
```
