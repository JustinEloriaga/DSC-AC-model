"""Small synthetic MATLAB v7.3 checkpoint fixtures for progress reporting."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import h5py
import numpy as np

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/report_chain_progress.py"
SPEC = importlib.util.spec_from_file_location("report_chain_progress", SCRIPT)
progress = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(progress)


class ChainProgressTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.run = Path(self.temporary.name) / "exact_run"
        self.cursor = Path(self.temporary.name) / "cursor.json"

    def checkpoint(self, count, stage="estimation", slot=1, seconds=None):
        base = self.run / "estimation" if stage == "estimation" else self.run
        path = base / "chains" / f"chain_{slot:03d}" / "checkpoint.mat"
        path.parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(path, "w") as saved:
            cp = saved.create_group("checkpoint")
            cp.attrs["MATLAB_class"] = np.bytes_("struct")
            cp["completed"] = [[float(count)]]
            cp["saved"] = [[float(max(0, count - 20))]]
            cfg = cp.create_group("identity").create_group("cfg")
            for key, value in dict(burnin=20, thin=1, chain_id=slot).items():
                cfg[key] = [[float(value)]]
            diagnostic = cp.create_group("diagnostics")
            if count:
                diagnostic["sweep_seconds"] = [seconds if seconds is not None else
                                                list(range(1, count + 1))]
            else:
                dataset = diagnostic.create_dataset("sweep_seconds", data=[0, 1])
                dataset.attrs["MATLAB_empty"] = 1
        return path

    def test_no_checkpoints_does_not_infer_another_run(self):
        report, _ = progress.snapshot(self.run)
        self.assertEqual(len(report["chains"]), 8)
        self.assertEqual(report["blocks"], [])
        self.assertTrue(all(c["completed_iterations"] is None for c in report["chains"]))
        self.assertFalse(self.run.exists())

    def test_exact_blocks_partial_tail_and_warmup(self):
        self.checkpoint(25)
        report, _ = progress.snapshot(self.run)
        self.assertEqual([b["sampler_seconds"] for b in report["blocks"]], [55, 155])
        self.assertEqual([b["remaining_iterations"] for b in report["blocks"]], [510, 500])
        self.assertEqual(report["chains"][0]["retained_draws"], 5)
        self.assertEqual(report["chains"][0]["remaining_iterations"], 495)
        self.assertEqual(report["blocks"][1]["current_completed_iterations"], 25)

    def test_empty_matlab_array_and_under_ten(self):
        self.checkpoint(0)
        self.checkpoint(9, slot=2)
        report, _ = progress.snapshot(self.run)
        self.assertEqual(report["blocks"], [])
        self.assertEqual(report["errors"], [])
        self.assertEqual(report["chains"][0]["completed_iterations"], 0)

    def test_no_duplicate_blocks_after_explicit_ack(self):
        self.checkpoint(20)
        report, cursor = progress.snapshot(self.run, cursor_path=self.cursor)
        self.assertEqual(len(report["blocks"]), 2)
        self.assertFalse(self.cursor.exists())
        again, _ = progress.snapshot(self.run, cursor_path=self.cursor)
        self.assertEqual(report, again)
        progress.acknowledge(self.cursor, cursor)
        report, _ = progress.snapshot(self.run, cursor_path=self.cursor)
        self.assertEqual(report["blocks"], [])
        self.checkpoint(31)
        report, _ = progress.snapshot(self.run, cursor_path=self.cursor)
        self.assertEqual([b["last_iteration"] for b in report["blocks"]], [30])
        self.assertEqual(report["blocks"][0]["sampler_seconds"], 255)
        self.assertEqual(report["blocks"][0]["retained_draws"], 10)

    def test_stage_and_chain_cursors_are_independent(self):
        self.checkpoint(20)
        _, cursor = progress.snapshot(self.run)
        progress.acknowledge(self.cursor, cursor)
        self.checkpoint(10, stage="smoothing")
        self.checkpoint(10, slot=2)
        report, _ = progress.snapshot(self.run, cursor_path=self.cursor)
        self.assertEqual([(b["stage"], b["chain_id"]) for b in report["blocks"]],
                         [("estimation", 2), ("smoothing", 1)])

    def test_twenty_warmup_five_hundred_retained_exact_target(self):
        for slot in range(1, 5):
            self.checkpoint(520, slot=slot, seconds=[2.5] * 520)
        report, _ = progress.snapshot(self.run)
        self.assertEqual(len(report["blocks"]), 208)
        self.assertTrue(all(c["retained_draws"] == 500 and c["remaining_iterations"] == 0
                            for c in report["chains"][:4]))
        self.assertTrue(all(b["sampler_seconds"] == 25 for b in report["blocks"]))
        self.assertEqual(report["blocks"][-1]["remaining_iterations"], 0)

    def test_corrupt_or_inconsistent_checkpoint_not_acknowledged(self):
        path = self.checkpoint(10)
        path.write_bytes(b"partial file")
        self.checkpoint(10, slot=2, seconds=[1] * 9)
        report, cursor = progress.snapshot(self.run)
        self.assertEqual(len(report["errors"]), 2)
        self.assertEqual(report["blocks"], [])
        self.assertEqual(cursor["positions"], {})

    def test_cursor_cannot_be_reused_for_another_run(self):
        _, cursor = progress.snapshot(self.run)
        progress.acknowledge(self.cursor, cursor)
        with self.assertRaisesRegex(ValueError, "another run"):
            progress.snapshot(self.run / "different", cursor_path=self.cursor)

    def test_cli_only_ack_writes_cursor(self):
        self.checkpoint(10)
        args = ["--run-dir", str(self.run), "--cursor", str(self.cursor)]
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(progress.main(args), 0)
        self.assertFalse(self.cursor.exists())
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(progress.main(args + ["--ack"]), 0)
        self.assertEqual(json.loads(self.cursor.read_text())["positions"]["estimation/chain_001"], 10)

    def test_precise_ack_keeps_blocks_completed_since_post_new(self):
        self.checkpoint(10)
        report, _ = progress.snapshot(self.run, cursor_path=self.cursor)
        self.assertEqual(report["blocks"][0]["last_iteration"], 10)
        # The user receives block 10, but block 20 finishes before acknowledgement.
        self.checkpoint(20)
        args = ["--run-dir", str(self.run), "--cursor", str(self.cursor),
                "--ack-block", "estimation:1:10"]
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(progress.main(args), 0)
        report, _ = progress.snapshot(self.run, cursor_path=self.cursor)
        self.assertEqual([b["last_iteration"] for b in report["blocks"]], [20])
        self.assertEqual(report["blocks"][0]["sampler_seconds"], 155)

    def test_precise_ack_rejects_unavailable_block(self):
        self.checkpoint(10)
        _, cursor = progress.snapshot(self.run)
        with self.assertRaisesRegex(ValueError, "unavailable block"):
            progress.acknowledge_blocks(self.cursor, cursor, ["estimation:1:20"])
        self.assertFalse(self.cursor.exists())


if __name__ == "__main__":
    unittest.main()
