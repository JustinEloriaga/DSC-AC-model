"""Publication provenance checks using isolated saved-artifact fixtures."""
from copy import deepcopy
from datetime import date, timedelta
import importlib.util
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/build_publication.py"
SPEC = importlib.util.spec_from_file_location("build_publication", SCRIPT)
publication = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(publication)


class InferencePublicationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.run = self.root / "run"
        self.figures = self.root / "generated"
        (self.run / "figures").mkdir(parents=True)
        self.figures.mkdir()
        self.args = SimpleNamespace(posterior_run=self.run, pilot=self.run / "pilot_summary.json")
        self.summary = {"source_hash": "synthetic", "first_return": "2020-01-03", "last_return": "2020-06-26"}
        self.sampler = {"saved_draws": 3, "burnin": 1, "convergence_established": False,
                        "first_date": "2020-01-03", "last_date": "2020-06-26",
                        "inference_mode": "fixed_parameter_smoothing"}
        self.metadata = {"inference_mode": "fixed_parameter_smoothing", "parameter_file": str(self.root / "parameters/model.mat"),
                         "parameters_estimated_at": "2026-09-26T11:22:33.000Z", "parameter_estimation_start": "2020-01-03",
                         "parameter_estimation_end": "2020-03-27", "parameter_estimation_draws": 5,
                         "parameter_estimation_run": str(self.root / "estimation"), "smoothing_end": "2020-06-26",
                         "parameter_smoothing": "mean", "parameter_uncertainty_in_bands": False}
        self.manifest = {"source_hash": "synthetic", "inference": deepcopy(self.metadata)}
        self.paths = {"retained_draws": 3, "warmup_removed": 1, "color_scale": [-1, 1],
                      "color_rule": "Continuous RdBu_r scale applied to the posterior median correlation.",
                      "run_dir": str(self.run), "first_date": "2020-01-03", "last_date": "2020-06-26",
                      "inference": deepcopy(self.metadata)}
        # PDF contents are opaque to make_bayesian_paths; it copies and hashes.
        # MATLAB plotting tests separately verify generated PDF figures.
        for name in ["selected"] + [f"all_{page:02d}" for page in range(1, 12)]:
            (self.run / f"figures/bayes_correlation_paths_{name}.pdf").write_bytes(b"%PDF-fixture\n")
        self.output = {}
        self.persist()

    def persist(self):
        for name, data in [("pilot_summary.json", self.sampler), ("inference_metadata.json", self.metadata),
                           ("run_manifest.json", self.manifest), ("figures/bayes_correlation_paths_manifest.json", self.paths)]:
            (self.run / name).write_text(json.dumps(data))

    def build_paths(self):
        return publication.make_bayesian_paths(self.args, self.output.__setitem__, self.figures, self.summary)

    def reserve_initial_weeks(self):
        self.summary.update(n_returns=26, n_series=3, n_pairs=3)
        self.metadata.update(calibration_weeks=4, calibration_start="2020-01-03",
                             calibration_end="2020-01-24", full_data_start="2020-01-03",
                             smoothing_start="2020-01-31", parameter_estimation_start="2020-01-31")
        self.sampler.update(T=22, m=3, pairs=3, first_date="2020-01-31", excluded_initial_weeks=4,
                            calibration_start="2020-01-03", calibration_end="2020-01-24",
                            completed_iterations=4, elapsed_seconds=8)
        self.paths["first_date"] = "2020-01-31"
        self.paths["inference"] = deepcopy(self.metadata)
        self.manifest["inference"] = deepcopy(self.metadata)
        self.persist()

    def pilot_fixture(self):
        folder = self.root / "data"
        folder.mkdir()
        (folder / "data_summary.json").write_text(json.dumps(self.summary))
        dates = [(date(2020, 1, 3) + timedelta(weeks=k)).isoformat() for k in range(26)]
        (folder / "weekly_returns.csv").write_text("Date\n" + "\n".join(dates) + "\n")
        prior = {"excluded_initial_weeks": 4,
                 "calibration": [{"window": 4, "sig2h_median": .001, "sig2r_median": .0001}]}
        (self.run / "prior_summary.json").write_text(json.dumps(prior))
        return folder

    def test_calibration_exclusion_is_published(self):
        self.reserve_initial_weeks()
        self.build_paths()
        text = self.output["bayesian_paths.tex"]
        for phrase in ["first 4 weekly returns", "2020-01-03", "2020-01-24", "excluded"]:
            self.assertIn(phrase, text)

    def test_calibration_count_mismatch_is_rejected(self):
        self.reserve_initial_weeks()
        self.sampler["T"] = 26
        self.persist()
        with self.assertRaisesRegex(ValueError, "return count"):
            self.build_paths()

    def test_calibration_boundary_mismatch_is_rejected(self):
        self.reserve_initial_weeks()
        self.metadata["calibration_end"] = "2020-01-17"
        self.paths["inference"] = deepcopy(self.metadata)
        self.manifest["inference"] = deepcopy(self.metadata)
        self.persist()
        with self.assertRaisesRegex(ValueError, "Calibration exclusion"):
            self.build_paths()

    def test_pilot_uses_modeled_suffix_and_discloses_prior_scope(self):
        self.reserve_initial_weeks()
        folder = self.pilot_fixture()
        publication.make_pilot(self.args, self.output.__setitem__, folder)
        self.assertIn("Weekly returns in the likelihood & 22", self.output["pilot.tex"])
        self.assertIn("2020-01-31 to 2020-06-26", self.output["pilot.tex"])
        self.assertIn("first 4 weekly returns", self.output["prior_sensitivity.tex"])
        self.assertIn("also enter estimation", self.output["prior_sensitivity.tex"])

    def test_pilot_rejects_reintroduced_calibration_rows(self):
        self.reserve_initial_weeks()
        folder = self.pilot_fixture()
        self.sampler.update(T=26, first_date="2020-01-03")
        self.persist()
        with self.assertRaisesRegex(ValueError, "Pilot T disagrees"):
            publication.make_pilot(self.args, self.output.__setitem__, folder)

    def test_fixed_scope_and_provenance_are_published(self):
        manifest = self.build_paths()
        self.assertEqual(manifest["inference"], self.metadata)
        self.assertIn("inference_metadata_sha256", manifest)
        self.assertEqual(len(manifest["figure_sha256"]), 12)
        text = self.output["bayesian_paths.tex"]
        for phrase in ["2020-03-27", "2020-06-26", "2026-09-26T11:22:33.000Z", "exclude parameter uncertainty",
                       "Historical state estimates can change", "dashed vertical line"]:
            self.assertIn(phrase, text)

    def test_draw_based_fixed_scope_includes_saved_parameter_variation(self):
        self.metadata.update(parameter_smoothing="draws", parameter_uncertainty_in_bands=True)
        self.paths["inference"] = deepcopy(self.metadata)
        self.manifest["inference"] = deepcopy(self.metadata)
        self.persist()
        self.build_paths()
        text = self.output["bayesian_paths.tex"]
        self.assertIn("saved retained parameter draws", text)
        self.assertIn("include variation across the saved parameter draws", text)
        self.assertNotIn("exclude parameter uncertainty", text)

    def test_joint_scope_includes_parameter_uncertainty(self):
        self.metadata.update(inference_mode="parameter_estimation", parameter_uncertainty_in_bands=True,
                             parameter_estimation_end=self.metadata["smoothing_end"])
        self.sampler["inference_mode"] = "parameter_estimation"
        self.paths["inference"] = deepcopy(self.metadata)
        self.manifest["inference"] = deepcopy(self.metadata)
        self.persist()
        self.build_paths()
        text = self.output["bayesian_paths.tex"]
        self.assertIn("joint parameter and state draws", text)
        self.assertNotIn("exclude parameter uncertainty", text)
        self.assertNotIn("dashed vertical line", text)

    def test_stale_figure_metadata_is_rejected_before_copying(self):
        self.paths["inference"]["parameter_estimation_end"] = "2020-02-28"
        self.persist()
        with self.assertRaisesRegex(ValueError, "disagrees"):
            self.build_paths()
        self.assertEqual(list(self.figures.iterdir()), [])

    def test_missing_conditional_metadata_is_rejected(self):
        (self.run / "inference_metadata.json").unlink()
        with self.assertRaisesRegex(ValueError, "metadata is missing"):
            self.build_paths()

    def test_mode_and_sample_mismatches_are_rejected(self):
        self.sampler["inference_mode"] = "parameter_estimation"
        self.persist()
        with self.assertRaisesRegex(ValueError, "different inference modes"):
            self.build_paths()
        self.sampler["inference_mode"] = "fixed_parameter_smoothing"
        self.paths["last_date"] = "2020-06-19"
        self.persist()
        with self.assertRaisesRegex(ValueError, "sample dates disagree"):
            self.build_paths()

    def test_legacy_no_metadata_run_still_builds(self):
        del self.sampler["inference_mode"]
        del self.paths["inference"]
        self.manifest["inference"] = "historical full-sample smoothing; empirical-Bayes priors"
        self.persist()
        (self.run / "inference_metadata.json").unlink()
        result = self.build_paths()
        self.assertNotIn("inference", result)
        self.assertIn("68\\% posterior bands", self.output["bayesian_paths.tex"])

    def test_prior_snapshot_is_bound_to_selected_run(self):
        shared = self.root / "data"
        shared.mkdir()
        (shared / "prior_summary.json").write_text('{"calibration_window":52}')
        with self.assertRaisesRegex(FileNotFoundError, "snapshot"):
            publication.selected_prior_summary(self.args.pilot, shared)
        snapshot = self.run / "prior_summary.json"
        snapshot.write_text('{"calibration_window":104}')
        self.assertEqual(publication.selected_prior_summary(self.args.pilot, shared), snapshot)
        snapshot.unlink()
        (self.run / "inference_metadata.json").unlink()
        self.assertEqual(publication.selected_prior_summary(self.args.pilot, shared), shared / "prior_summary.json")

    def test_posterior_run_alone_selects_its_own_pilot_and_priors(self):
        pilot = publication.selected_pilot_summary(None, self.run)
        self.assertEqual(pilot, self.run / "pilot_summary.json")
        snapshot = self.run / "prior_summary.json"
        snapshot.write_text('{"calibration_window":104}')
        self.assertEqual(publication.selected_prior_summary(pilot, self.root / "data"), snapshot)
        self.assertIsNone(publication.selected_pilot_summary(None, None))
        self.assertEqual(publication.selected_pilot_summary(self.args.pilot, self.run), self.args.pilot)


if __name__ == "__main__":
    unittest.main()
