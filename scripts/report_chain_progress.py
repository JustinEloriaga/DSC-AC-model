#!/usr/bin/env python3
"""Read this run's committed MATLAB checkpoints; never modify model outputs.

Times are sums of sampler sweep timers, not wall time: checkpoint writes,
posterior storage, diagnostics and report generation are excluded.  A cursor
is updated only explicitly. Use --ack-block STAGE:SLOT:END to acknowledge the
exact block endpoints delivered, even if newer iterations have since finished.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import sys
import tempfile

try:
    import h5py
except ImportError:
    h5py = None


def scalar(group, name):
    value = group[name][()]
    if value.size != 1:
        raise ValueError(f"{name} must be scalar")
    number = float(value.reshape(-1)[0])
    if not math.isfinite(number) or number < 0 or number != int(number):
        raise ValueError(f"{name} must be a nonnegative integer")
    return int(number)


def read_checkpoint(path, total_iterations):
    with h5py.File(path, "r") as saved:
        checkpoint = saved["checkpoint"]
        cfg = checkpoint["identity"]["cfg"]
        completed = scalar(checkpoint, "completed")
        retained = scalar(checkpoint, "saved")
        burnin = scalar(cfg, "burnin")
        thin = scalar(cfg, "thin")
        chain_id = scalar(cfg, "chain_id")
        timing = checkpoint["diagnostics"]["sweep_seconds"]
        # MATLAB v7.3 represents an empty array using its dimensions as data.
        times = [] if bool(timing.attrs.get("MATLAB_empty", False)) else [
            float(value) for value in timing[()].reshape(-1)
        ]
    if thin < 1 or chain_id < 1 or burnin > total_iterations:
        raise ValueError("Invalid chain identity, warm-up or thinning")
    if completed > total_iterations or len(times) != completed:
        raise ValueError("Completed count disagrees with target or sweep timers")
    if any(not math.isfinite(value) or value < 0 for value in times):
        raise ValueError("Sweep timers must be finite and nonnegative")
    if retained != max(0, (completed - burnin) // thin):
        raise ValueError("Saved count disagrees with completed/warm-up/thinning")
    return dict(chain_id=chain_id, completed_iterations=completed,
                remaining_iterations=total_iterations - completed,
                warmup_iterations=burnin, warmup_completed=min(burnin, completed),
                retained_draws=retained,
                retained_target=max(0, (total_iterations - burnin) // thin),
                thin=thin, times=times)


def read_cursor(cursor_path, identity):
    if cursor_path is None or not cursor_path.exists():
        return {}
    data = json.loads(cursor_path.read_text())
    if data.get("identity") != identity:
        raise ValueError("Cursor belongs to another run, target or block size")
    positions = data.get("positions")
    if not isinstance(positions, dict) or any(
        not isinstance(value, int) or value < 0
        or value % identity["block_size"] or value > identity["total_iterations"]
        for value in positions.values()
    ):
        raise ValueError("Invalid cursor positions")
    return positions


def snapshot(run_dir, total_iterations=520, block_size=10, num_chains=4,
             cursor_path=None):
    if h5py is None:
        raise RuntimeError("h5py is unavailable; install nothing automatically")
    if min(total_iterations, block_size, num_chains) < 1:
        raise ValueError("Iteration, block and chain counts must be positive")
    run_dir = Path(run_dir).resolve()
    cursor_path = Path(cursor_path) if cursor_path is not None else None
    identity = dict(run_dir=str(run_dir), total_iterations=total_iterations,
                    block_size=block_size, num_chains=num_chains)
    previous = read_cursor(cursor_path, identity)
    next_positions = dict(previous)
    result = dict(**identity, timing_basis="sum_of_committed_sampler_sweep_seconds",
                  timing_excludes=["posterior storage", "checkpoint writes",
                                   "convergence diagnostics", "report generation"],
                  blocks=[], chains=[], errors=[])
    for stage, base in [("estimation", run_dir / "estimation"),
                        ("smoothing", run_dir)]:
        for slot in range(1, num_chains + 1):
            path = base / "chains" / f"chain_{slot:03d}" / "checkpoint.mat"
            key = f"{stage}/chain_{slot:03d}"
            status = dict(stage=stage, chain_slot=slot, checkpoint=str(path),
                          status="waiting_for_checkpoint", completed_iterations=None,
                          remaining_iterations=None)
            result["chains"].append(status)
            if not path.is_file():
                continue
            try:
                current = read_checkpoint(path, total_iterations)
                times = current.pop("times")
                if previous.get(key, 0) > current["completed_iterations"]:
                    raise ValueError("Checkpoint is behind the acknowledged cursor")
                status.update(current)
                status["status"] = "iteration_target_reached" if not current[
                    "remaining_iterations"] else "in_progress"
                final_block = current["completed_iterations"] // block_size * block_size
                for end in range(previous.get(key, 0) + block_size,
                                 final_block + 1, block_size):
                    result["blocks"].append(dict(
                        stage=stage, chain_slot=slot, chain_id=current["chain_id"],
                        first_iteration=end - block_size + 1, last_iteration=end,
                        sampler_seconds=math.fsum(times[end - block_size:end]),
                        remaining_iterations=total_iterations - end,
                        warmup_completed=min(current["warmup_iterations"], end),
                        retained_draws=max(0, (end - current["warmup_iterations"]) // current["thin"]),
                        current_completed_iterations=current["completed_iterations"],
                        current_remaining_iterations=current["remaining_iterations"]))
                next_positions[key] = final_block
            except (OSError, KeyError, ValueError, TypeError) as error:
                # Atomic replacement may race a read. Do not invent progress or
                # acknowledge a block; the next invocation can retry safely.
                status["status"] = "checkpoint_unreadable"
                status["error"] = str(error)
                result["errors"].append(dict(stage=stage, chain_slot=slot,
                                             checkpoint=str(path), error=str(error)))
    return result, dict(identity=identity, positions=next_positions)


def acknowledge(path, cursor):
    """Atomically save the explicitly acknowledged blocks, not model files."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent,
                                     prefix=path.name + ".", delete=False) as handle:
        temporary = Path(handle.name)
        try:
            json.dump(cursor, handle, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def acknowledge_blocks(path, candidate, specifications):
    """Advance only supplied block endpoints, never newer unseen blocks."""
    identity = candidate["identity"]
    positions = read_cursor(Path(path), identity)
    for specification in specifications:
        try:
            stage, slot_text, end_text = specification.split(":")
            slot, end = int(slot_text), int(end_text)
        except ValueError as error:
            raise ValueError("--ack-block must be STAGE:SLOT:END, e.g. estimation:1:10") from error
        key = f"{stage}/chain_{slot:03d}"
        if (stage not in {"estimation", "smoothing"}
                or not 1 <= slot <= identity["num_chains"]
                or end < 0 or end % identity["block_size"]
                or end > candidate["positions"].get(key, 0)):
            raise ValueError(f"Cannot acknowledge unavailable block {specification}")
        positions[key] = max(positions.get(key, 0), end)
    acknowledge(path, dict(identity=identity, positions=positions))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--total-iterations", type=int, default=520)
    parser.add_argument("--block-size", type=int, default=10)
    parser.add_argument("--num-chains", type=int, default=4)
    parser.add_argument("--cursor", type=Path)
    acknowledgement = parser.add_mutually_exclusive_group()
    acknowledgement.add_argument("--ack", action="store_true",
                                 help="Acknowledge ALL blocks emitted by this call; prefer --ack-block after posting")
    acknowledgement.add_argument("--ack-block", action="append", default=[],
                                 metavar="STAGE:SLOT:END",
                                 help="Acknowledge only this delivered endpoint; repeat per chain/stage")
    args = parser.parse_args(argv)
    if (args.ack or args.ack_block) and args.cursor is None:
        parser.error("Acknowledgement requires --cursor")
    try:
        report, cursor = snapshot(args.run_dir, args.total_iterations,
                                  args.block_size, args.num_chains, args.cursor)
        print(json.dumps(report, indent=2, allow_nan=False), flush=True)
        if args.ack:
            acknowledge(args.cursor, cursor)
        elif args.ack_block:
            acknowledge_blocks(args.cursor, cursor, args.ack_block)
        return 0
    except (OSError, ValueError, RuntimeError) as error:
        print(json.dumps({"error": str(error)}), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
