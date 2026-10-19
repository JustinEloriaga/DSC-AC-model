"""Pure display-limit checks and six-axis in-memory renderer regression tests."""
from datetime import date, timedelta
import os
from pathlib import Path
import types
import unittest

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/replot_correlation_axes.py"
renderer = types.ModuleType("replot_correlation_axes")
renderer.__file__ = str(SCRIPT)
# Avoid importlib's optional __pycache__ writes, including for the test subject.
exec(compile(SCRIPT.read_bytes(), str(SCRIPT), "exec"), renderer.__dict__)


class CorrelationPanelLimitTests(unittest.TestCase):
    def test_all_bands_and_fits_enter_the_range(self):
        baseline = np.array([[.10, .11], [.20, .21], [.30, .31]])
        loose = np.array([[-.40, -.39], [.25, .26], [.70, .71]])
        np.testing.assert_allclose(renderer.panel_limits(baseline, loose), [-.41, .72], rtol=0, atol=1e-16)

    def test_nonfinite_values_are_ignored(self):
        self.assertEqual(renderer.panel_limits([np.nan, .2, np.inf], [-np.inf, .5]), (.19, .51))

    def test_constant_panel_gets_exact_padding(self):
        np.testing.assert_allclose(renderer.panel_limits(np.full((3, 4), .4)), [.39, .41], rtol=0, atol=1e-16)

    def test_negative_and_cross_zero_ranges(self):
        np.testing.assert_allclose(renderer.panel_limits([-.8, -.3]), [-.81, -.29], rtol=0, atol=1e-16)
        self.assertEqual(renderer.panel_limits([-.1, .3]), (-.11, .31))

    def test_no_zero_inclusion_or_correlation_clamping(self):
        np.testing.assert_allclose(renderer.panel_limits([.8, 1]), [.79, 1.01], rtol=0, atol=1e-16)
        np.testing.assert_allclose(renderer.panel_limits([-1, -.8]), [-1.01, -.79], rtol=0, atol=1e-16)

    def test_empty_or_all_nonfinite_panel_fails(self):
        for arrays in ((), ([],), ([np.nan, np.inf, -np.inf],), ([], [np.nan])):
            with self.subTest(arrays=arrays), self.assertRaisesRegex(ValueError, "finite displayed"):
                renderer.panel_limits(*arrays)

    def test_invalid_padding_fails(self):
        for padding in (0, -1, np.nan, np.inf):
            with self.subTest(padding=padding), self.assertRaisesRegex(ValueError, "padding"):
                renderer.panel_limits([.2], padding=padding)

    def test_six_actual_axes_keep_independent_limits_after_draw(self):
        # Caller supplies an existing cache location; no fixture directories,
        # images, PDFs, or other test outputs are authored.
        launch = Path(os.environ.get("MPLCONFIGDIR", ""))
        plotting = renderer.load_plotting(launch)
        dates = [date(2020, 1, 3) + timedelta(weeks=i) for i in range(4)]
        baseline = np.empty((3, 6, 4))
        loose = np.empty_like(baseline)
        for k in range(6):
            center = -.6 + k * .2
            baseline[:, k, :] = np.array([center - .04, center, center + .05])[:, None]
            loose[:, k, :] = np.array([center - .08, center + .01, center + .07])[:, None]

        def stage(label, bands):
            return {"label": label, "bands": bands, "dates": dates, "count": 1000,
                    "estimation_count": 1000, "provenance": {"target_complete": True},
                    "inference": {"smoothing_start": dates[0].isoformat(), "smoothing_end": dates[-1].isoformat()}}

        stages = (stage("baseline", baseline), stage("loose", loose))
        original = [s["bands"].copy() for s in stages]
        for selected in ((stages[0],), (stages[1],), stages):
            with self.subTest(fits=[s["label"] for s in selected]):
                fig, records = renderer.make_figure(selected, "Test correlations", plotting)
                try:
                    fig.canvas.draw()
                    self.assertEqual(len(fig.axes), 6)
                    for k, ax in enumerate(fig.axes):
                        expected = renderer.panel_limits(*(s["bands"][:, k, :] for s in selected))
                        np.testing.assert_array_equal(ax.get_ylim(), expected)
                        np.testing.assert_array_equal(records[k]["ylim"], expected)
                        for other in fig.axes:
                            if other is not ax:
                                self.assertFalse(ax.get_shared_y_axes().joined(ax, other))
                    self.assertEqual(len({ax.get_ylim() for ax in fig.axes}), 6)
                finally:
                    plotting.plt.close(fig)
        for saved, unchanged in zip(stages, original):
            np.testing.assert_array_equal(saved["bands"], unchanged)


if __name__ == "__main__":
    unittest.main()
