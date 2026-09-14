#!/usr/bin/env python3
"""Summarize saved kernel benchmarks and a completed pilot; never starts MATLAB."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re
import statistics

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--benchmark", type=Path)
    parser.add_argument("--resources", type=Path)
    args = parser.parse_args()
    run = args.run.resolve()
    if args.benchmark is None:
        args.benchmark = run / "kernel_benchmark.json"
        if not args.benchmark.exists():
            args.benchmark = ROOT / "outputs/weekly_research/performance/kernel_benchmark.json"
    if args.resources is None and (run / "process_resources.txt").exists():
        args.resources = run / "process_resources.txt"
    latest = json.loads((run / "pilot_summary.json").read_text())
    benchmark = json.loads(args.benchmark.read_text())
    baseline_path = Path(benchmark["checkpoint"]).parent / "pilot_summary.json"
    baseline = json.loads(baseline_path.read_text())
    validation = json.loads((run / "validation.json").read_text())
    likelihood_validation = json.loads((run / "final_likelihood_validation.json").read_text())
    assert latest["status"] in ("iteration_limit", "time_limit")
    assert validation["status"] == "passed"
    assert likelihood_validation["status"] == "passed"
    assert latest["implementation"] == benchmark["implementation"], "Benchmark and pilot implementations differ"
    assert benchmark["source_hash"] == validation["source_hash"]
    for key in ("seed", "T", "m", "pairs", "first_date", "last_date"):
        assert baseline[key] == latest[key], f"Historical comparison mismatch: {key}"
    common = min(len(benchmark["baseline_sweep_seconds"]), len(latest["sweep_seconds"]))
    assert common > 0
    old_seconds = statistics.mean(benchmark["baseline_sweep_seconds"][:common])
    new_seconds = statistics.mean(latest["sweep_seconds"][:common])
    proposal_counts_match = all(
        benchmark["baseline_" + key][:common] == latest[key][:common]
        for key in ("r_proposals_by_sweep", "h_proposals_by_sweep")
    )
    result = {
        "baseline_summary": str(baseline_path), "new_summary": str(run / "pilot_summary.json"),
        "kernel_speedup": benchmark["kernel_speedup"],
        "common_completed_sweeps": common, "baseline_common_mean_seconds": old_seconds,
        "new_common_mean_seconds": new_seconds, "common_sweep_speedup": old_seconds / new_seconds,
        "common_proposal_counts_match": proposal_counts_match,
        "new_completed_sweeps": latest["completed_iterations"],
        "new_sampler_seconds": latest["total_seconds"],
        "new_mean_sweep_seconds": latest["mean_sweep_seconds"],
        "illustrative_4000_sweep_chain_hours": latest["mean_sweep_seconds"] * 4000 / 3600,
        "illustrative_four_serial_chains_hours": latest["mean_sweep_seconds"] * 16000 / 3600,
        "saved_posterior_draws": latest["saved_draws"],
        "convergence_established": latest["convergence_established"],
        "numerical_validation": validation,
        "final_likelihood_validation": likelihood_validation,
    }
    memory_line = "Peak process memory was not supplied."
    if args.resources:
        resources = args.resources.read_text()
        matched = re.search(r"(\d+)\s+maximum resident set size", resources)
        if matched:
            peak_bytes = int(matched.group(1))
            result["peak_process_rss_gib"] = peak_bytes / 2**30
            result["memory_method"] = "macOS /usr/bin/time -l; peak MATLAB process RSS in bytes, includes data preparation and validation"
            memory_line = f"Peak MATLAB process resident memory: {peak_bytes / 2**30:.3f} GiB (macOS `time -l`, including preparation and validation)."
        (run / "process_resources.txt").write_text(resources)
    (run / "performance_comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    (run / "kernel_benchmark.json").write_text(json.dumps(benchmark, indent=2) + "\n")
    test_path = run / "acceptance_tests.json"
    if not test_path.exists():
        test_path = ROOT / "outputs/weekly_research/matlab_test_results.json"
    test_evidence = json.loads(test_path.read_text())
    assert test_evidence["failed"] == 0 and test_evidence["incomplete"] == 0
    (run / "acceptance_tests.json").write_text(json.dumps(test_evidence, indent=2) + "\n")
    stages = "\n".join(
        f"| {name.replace('_', ' ').capitalize()} | {baseline['stage_seconds'][name]:.2f} | {latest['stage_seconds'][name]:.2f} |"
        for name in latest["stage_seconds"]
    )
    note = f"""# MATLAB sampler speed review and rerun

The optimized pilot completed {latest['completed_iterations']} sweeps in {latest['total_seconds']/60:.2f} minutes. Its mean sweep cost was {latest['mean_sweep_seconds']:.2f} seconds. On the first {common} completed sweeps, the old implementation averaged {old_seconds:.2f} seconds and the new implementation averaged {new_seconds:.2f} seconds: **{old_seconds/new_seconds:.2f}x faster**. These full-sweep measurements were taken in separate sessions on the same machine; the fixed-input microbenchmark below compares both kernels in one session.

## What changed

- A batched C++ MEX kernel calls MATLAB's LAPACK/BLAS for correlation reconstruction and likelihood evaluation.
- Fixed-point iterations calculate only the diagonal of the matrix exponential, and the full matrix is constructed after convergence.
- {latest['correlation_threads']} native threads evaluate independent date blocks; likelihood terms are summed in date order. The sequential MCMC coordinate updates and random draws are unchanged.
- Standardized return innovations are calculated once per correlation sweep.
- Resolved backend, implementation hashes, invalid proposals and per-sweep timings are recorded. Checkpoints reject a changed implementation.

## Evidence

The fixed saved-state kernel benchmark improved from {benchmark['median_matlab_seconds']:.4f} to {benchmark['median_native_seconds']:.4f} seconds (**{benchmark['kernel_speedup']:.2f}x**), with absolute log-likelihood difference {benchmark['absolute_likelihood_error']:.3g}. Common-sweep correlation and volatility proposal counts match the historical run: **{proposal_counts_match}**. Matching counts are a reproducibility check, not proof of identical floating-point states.

The rerun retained all {latest['T']} weeks, {latest['m']} variables and {latest['pairs']} pairs. All {validation['matrices_checked']} final-state correlation matrices passed symmetry, unit-diagonal, positive-definiteness and matrix-log roundtrip checks; maximum coordinate error was {validation['maximum_coordinate_roundtrip_error']:.3g}. Invalid completed-sweep correlation proposals: {latest['r_invalid_proposals']}; invalid volatility proposals: {latest['h_invalid_proposals']}.

The final saved likelihood also agrees with the independent MATLAB calculation: absolute difference {likelihood_validation['absolute_errors']['saved_vs_matlab']:.3g}, against a required tolerance of {likelihood_validation['absolute_tolerance']:.1g}. This validation consumes no random draws and does not advance the chain.

{memory_line}

| Stage | Historical total ({baseline['completed_iterations']} sweeps), seconds | New total ({latest['completed_iterations']} sweeps), seconds |
|---|---:|---:|
{stages}

## Production implications

At the measured new mean, 4,000 sweeps would take approximately **{result['illustrative_4000_sweep_chain_hours']:.1f} hours per chain**, or {result['illustrative_four_serial_chains_hours']:.1f} hours for four serial chains. These are linear compute projections from warm-up: mean sweep timings exclude checkpoint writes, and this pilot does not measure production posterior-chunk storage. Later sampling and simultaneous chains can have different costs. Native date threads already share the CPU, so concurrent chains are not assumed to gain another fourfold speedup.

All {latest['completed_iterations']} sweeps are warm-up. There are **{latest['saved_draws']} retained posterior draws**, and convergence has not been established. The rerun therefore does not yet supply Bayesian correlation estimates or 68% posterior bands. Production still requires its own run budget and the agreed convergence and sensitivity checks.

## Saved artifacts

- New run: `{run.relative_to(ROOT)}`
- Historical run: `{baseline_path.parent.relative_to(ROOT)}`
- Kernel benchmark: `kernel_benchmark.json` in the new run (copied from `{args.benchmark.resolve().relative_to(ROOT)}`).
- Acceptance results: `acceptance_tests.json` in the new run ({test_evidence['passed']} passed).
- Numerical validation: `validation.json` in the new run.

Rebuild this summary without rerunning estimation:

```sh
python3 scripts/summarize_speed_pilot.py --run {run.relative_to(ROOT)}
```
"""
    (ROOT / "research/sampler_performance.md").write_text(note)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
