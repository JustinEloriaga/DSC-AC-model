#!/usr/bin/env python3
"""Write new axis-rescaled revisions using validated saved correlations only.

Usage: python -B scripts/replot_correlation_axes.py LAUNCH REPORT --check-only
Then set MPLCONFIGDIR to the same existing LAUNCH and omit --check-only.
No estimation, smoothing, MATLAB, directory creation, or original-file update
occurs. The original SHA-bound comparison helper is executed in memory, and
its strict validate() supplies the unchanged MATLAB-compatible Hazen bands.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import re
import types

import numpy as np


PADDING = 0.01
QUANTILES = [16, 50, 84]
QUANTILE_METHOD = "Hazen/type 5; MATLAB R2024b prctile"
TICKERS = ["NKYTR Index", "SPXT Index", "SX5T Index", "TUKXG Index"]
PAIR_I = [2, 3, 4, 3, 4, 4]
PAIR_J = [1, 1, 1, 2, 2, 3]
COLORS = {"baseline": "#1763A6", "loose": "#D47724"}
LABELS = {"baseline": "Baseline posterior", "loose": "Looser-prior posterior"}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def existing_directory(path, label):
    require(path.is_absolute() and path.is_dir(), f"Use an existing absolute {label}: {path}")
    require(not any(part.is_symlink() for part in (path, *path.parents)),
            f"Linked directory/ancestor is not allowed for {label}: {path}")


def finite_extrema(*displayed):
    """Extrema of finite displayed values; NaN/Inf never enlarge the range."""
    lows, highs = [], []
    for values in displayed:
        array = np.asarray(values, dtype=float)
        finite = array[np.isfinite(array)]
        if finite.size:
            lows.append(float(np.min(finite)))
            highs.append(float(np.max(finite)))
    require(bool(lows), "A correlation panel needs at least one finite displayed value")
    return min(lows), max(highs)


def panel_limits(*displayed, padding=PADDING):
    """Return min(displayed)-padding, max(displayed)+padding, without clamping.

    Pass every displayed median/band array from every fit in a panel. Raw
    posterior extrema, invisible zero lines, and other panels are not inputs.
    """
    require(np.isscalar(padding) and np.isfinite(padding) and padding > 0,
            "Panel padding must be finite and positive")
    low, high = finite_extrema(*displayed)
    limits = (low - float(padding), high + float(padding))
    require(all(np.isfinite(limits)) and limits[0] < limits[1], "Invalid padded panel limits")
    return limits


def pair_label(k):
    return (f"{TICKERS[PAIR_I[k] - 1].replace(' Index', '')} / "
            f"{TICKERS[PAIR_J[k] - 1].replace(' Index', '')}")


def panel_records(stages):
    records = []
    for k in range(6):
        displayed = [stage["bands"][:, k, :] for stage in stages]
        low, high = finite_extrema(*displayed)
        records.append({"pair_i": PAIR_I[k], "pair_j": PAIR_J[k], "pair": pair_label(k),
                        "displayed_min": low, "displayed_max": high,
                        "ylim": list(panel_limits(*displayed)),
                        "fits": [stage["label"] for stage in stages]})
    return records


def figure_specs(launch, baseline, loose):
    match = re.fullmatch(r"loose-prior-mean\.(\d{14})", launch.name)
    require(match is not None, "Expected an approved loose-prior-mean.TIMESTAMP launch")
    suffix = match[1]
    return [
        (f"baseline_posterior_correlations_zoomed_{suffix}",
         "Baseline posterior correlations", (baseline,)),
        (f"loose_prior_posterior_correlations_zoomed_{suffix}",
         "Looser-prior posterior correlations", (loose,)),
        (f"four_ticker_loose_prior_sensitivity_zoomed_{suffix}",
         "Prior-concentration sensitivity", (baseline, loose)),
    ]


def outputs_for(launch, report, specs):
    paths = [report / f"{stem}.{extension}" for stem, _, _ in specs for extension in ("pdf", "png")]
    paths.append(launch / "axis_rescale_provenance.json")
    require(not any(path.exists() or path.is_symlink() for path in paths),
            "Axis-rescale output exists/is linked; saved artifacts must not be overwritten")
    return paths


def load_saved(launch, report):
    existing_directory(launch, "launch directory")
    existing_directory(report, "report directory")
    original_proof_path = launch / "comparison_provenance.json"
    helper_path = launch / "compare_loose_saved.py"
    for path in (original_proof_path, helper_path):
        require(path.is_file() and not path.is_symlink(), f"Missing/linked saved input: {path}")
    original = json.loads(original_proof_path.read_text(encoding="utf-8"))
    require(isinstance(original, dict), "Original comparison provenance must be an object")
    helper_sha = sha256(helper_path)
    require(helper_sha == original.get("helper_sha256"),
            "Original helper differs from its saved comparison provenance")
    require(original.get("quantiles_percent") == QUANTILES
            and original.get("quantile_method") == QUANTILE_METHOD
            and original.get("static_parameter_uncertainty_included") is False
            and original.get("convergence_checked") is False,
            "Unexpected original comparison quantiles/inference scope")
    require(Path(original.get("launch_dir", "")).resolve() == launch.resolve(),
            "Original comparison provenance belongs to another launch")
    module = types.ModuleType("_frozen_loose_saved_reader")
    module.__file__ = str(helper_path)
    # No importlib, __pycache__, or change to the frozen helper/provenance.
    exec(compile(helper_path.read_bytes(), str(helper_path), "exec"), module.__dict__)
    require(module.TICKERS == TICKERS and list(module.PAIR_I) == PAIR_I
            and list(module.PAIR_J) == PAIR_J, "Frozen helper ticker/pair map differs")
    state, state_path, baseline, loose, proof, rows = module.validate(launch, report)
    require(sha256(state_path) == original.get("launch_status_sha256"),
            "Launch status differs from original saved comparison")
    require(original.get("stages") == [baseline["provenance"], loose["provenance"]],
            "Saved chunks/metadata differ from original comparison provenance")
    require(original.get("validation", {}).get("parameter_artifacts") == proof["parameter_artifacts"],
            "Saved parameter artifacts differ from original comparison provenance")
    specs = figure_specs(launch, baseline, loose)
    outputs_for(launch, report, specs)
    return {"state": state, "state_path": state_path, "baseline": baseline, "loose": loose,
            "validation": proof, "parameter_scalar_count": len(rows), "specs": specs,
            "original_comparison": {"helper_path": str(helper_path), "helper_sha256": helper_sha,
                                    "provenance_path": str(original_proof_path),
                                    "provenance_sha256": sha256(original_proof_path)}}


def load_plotting(launch):
    """Validate the existing cache location before importing matplotlib."""
    existing_directory(launch, "launch directory")
    mpldir = Path(os.environ.get("MPLCONFIGDIR", ""))
    existing_directory(mpldir, "MPLCONFIGDIR")
    require(mpldir.resolve() == launch.resolve(), "Set MPLCONFIGDIR to the existing launch directory")
    require(os.access(mpldir, os.W_OK), "MPLCONFIGDIR must be writable; refuse matplotlib cache fallback")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator
    return types.SimpleNamespace(plt=plt, dates=mdates, Line2D=Line2D, MaxNLocator=MaxNLocator)


def make_figure(stages, title, plotting):
    """Build six independent axes in memory; do not save or alter band data."""
    records = panel_records(stages)
    fig, axes = plotting.plt.subplots(3, 2, figsize=(12, 10.4), sharex=True, sharey=False)
    fig.subplots_adjust(left=.075, right=.98, top=.85, bottom=.215, hspace=.30, wspace=.21)
    for k, ax in enumerate(axes.flat):
        for stage in stages:
            q, color = stage["bands"], COLORS[stage["label"]]
            ax.fill_between(stage["dates"], q[0, k], q[2, k], color=color, alpha=.16, linewidth=0)
            ax.plot(stage["dates"], q[1, k], color=color, linewidth=1.05)
        ax.set_title(pair_label(k), loc="left", fontsize=11, fontweight="medium")
        ax.set_xlim(stages[0]["dates"][0], stages[0]["dates"][-1])
        ax.yaxis.set_major_locator(plotting.MaxNLocator(nbins=5))
        ax.grid(axis="y", color="#DCE1E6", linewidth=.55)
        ax.xaxis.set_major_locator(plotting.dates.YearLocator(5))
        ax.xaxis.set_major_formatter(plotting.dates.DateFormatter("%Y"))
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=9)
        if k % 2 == 0:
            ax.set_ylabel("Correlation", fontsize=10)
        # Set limits last: neither tick placement nor zero lines may expand them.
        ax.set_ylim(*records[k]["ylim"])
    complete = all(stage["provenance"]["target_complete"] for stage in stages)
    if not complete:
        title += " - PARTIAL TARGETS"
    fig.suptitle(title, x=.075, y=.969, ha="left", fontsize=20, fontweight="bold")
    fig.text(.075, .931, "Four-ticker historical correlations, conditional on each fitted posterior parameter mean",
             ha="left", fontsize=11.1, color="#455565")
    labels = [f"{LABELS[s['label']]}: {s['count']:,} paths"
              + (" (partial)" if not s["provenance"]["target_complete"] else "") for s in stages]
    fig.legend([plotting.Line2D([0], [0], color=COLORS[s["label"]], lw=2) for s in stages], labels,
               loc="upper left", bbox_to_anchor=(.067, .914), ncols=len(stages), frameon=False, fontsize=10)
    counts = "; ".join(f"{s['label']} {s['estimation_count']:,} parameter draws / {s['count']:,} paths" for s in stages)
    notes = [
        "Lines: pointwise medians. Shading: 16th-84th percentiles; MATLAB-compatible Hazen quantiles.",
        "Panel y-axis bounds differ: displayed band/median minimum - 0.01 to maximum + 0.01; no clamping.",
        f"Modeled weeks: {stages[0]['inference']['smoothing_start']} to {stages[0]['inference']['smoothing_end']}; initial 104 calibration weeks excluded.",
        f"Actual retained counts: {counts}.",
        "Bands are conditional on posterior parameter means and exclude static parameter uncertainty.",
        "Convergence mode: off; one chain; convergence not checked. Exploratory historical state estimates.",
    ]
    if len(stages) == 2:
        notes.append("Data and static prior means unchanged; IG shape 10 -> 3 and IW df 104 -> 20 change prior concentration.")
    if not complete:
        notes.append("Partial results use actual saved draws/paths only; requested targets were not fully completed.")
    fig.text(.075, .163, "\n".join(notes), ha="left", va="top", fontsize=9, color="#455565", linespacing=1.5)
    return fig, records


def counts(stage):
    return {"fit": stage["label"], "retained_parameter_draws": stage["estimation_count"],
            "retained_paths": stage["count"], "modeled_weeks": len(stage["dates"]),
            "target_complete": stage["provenance"]["target_complete"]}


def check_summary(saved):
    return {"validated": True, "saved_results_only": True, "padding": PADDING,
            "quantiles_percent": QUANTILES, "quantile_method": QUANTILE_METHOD,
            "counts": [counts(saved[key]) for key in ("baseline", "loose")],
            "figures": [{"stem": stem, "panels": panel_records(stages)} for stem, _, stages in saved["specs"]]}


def author(launch, report, saved):
    outputs = outputs_for(launch, report, saved["specs"])
    plotting = load_plotting(launch)
    figures, artifacts = [], []
    for stem, title, stages in saved["specs"]:
        fig, panels = make_figure(stages, title, plotting)
        try:
            for extension in ("pdf", "png"):
                path = report / f"{stem}.{extension}"
                options = ({"metadata": {"Title": title,
                            "Subject": "Saved posterior-mean smoothing; convergence not checked",
                            "Creator": "replot_correlation_axes.py (saved results only)"}}
                           if extension == "pdf" else {"dpi": 160})
                # Exclusive file creation also protects against an output appearing after preflight.
                with path.open("xb") as stream:
                    fig.savefig(stream, format=extension, **options)
                artifacts.append({"path": str(path), "sha256": sha256(path)})
            figures.append({"stem": stem, "panels": panels})
        finally:
            plotting.plt.close(fig)
    stages = [saved[key] for key in ("baseline", "loose")]
    provenance = {"schema_version": 1, "created_at_utc": datetime.now(timezone.utc).isoformat(),
                  "operation": "Rescale figures from validated direct stage-2 retained chunks only; no estimation or smoothing",
                  "launch_dir": str(launch), "launch_status_sha256": sha256(saved["state_path"]),
                  "renderer_path": str(Path(__file__).resolve()), "renderer_sha256": sha256(Path(__file__)),
                  "original_comparison": saved["original_comparison"],
                  "source_hash": saved["baseline"]["manifest"]["source_hash"],
                  "quantiles_percent": QUANTILES, "quantile_method": QUANTILE_METHOD,
                  "padding": PADDING, "limits_rule": "finite displayed minimum - 0.01, maximum + 0.01, independently per panel",
                  "tickers": TICKERS, "pair_i": PAIR_I, "pair_j": PAIR_J,
                  "modeled_dates": [d.isoformat() for d in stages[0]["dates"]],
                  "counts": [counts(stage) for stage in stages], "figures": figures,
                  "static_parameter_uncertainty_included": False, "convergence_checked": False,
                  "convergence_mode": "off", "conditional_on": "each fit's posterior static parameter mean",
                  "validation": saved["validation"], "stages": [stage["provenance"] for stage in stages],
                  "displayed_band_arrays": [{"fit": stage["label"], "shape": list(stage["bands"].shape),
                                            "sha256_float64_little_endian_c_order": hashlib.sha256(
                                                np.asarray(stage["bands"], dtype="<f8", order="C").tobytes()).hexdigest()}
                                           for stage in stages],
                  "artifacts": artifacts}
    with outputs[-1].open("x", encoding="utf-8") as stream:
        json.dump(provenance, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return {"figures": artifacts, "provenance": str(outputs[-1])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("launch_dir", type=Path)
    parser.add_argument("report_dir", type=Path)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    saved = load_saved(args.launch_dir, args.report_dir)
    result = check_summary(saved) if args.check_only else author(args.launch_dir, args.report_dir, saved)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
