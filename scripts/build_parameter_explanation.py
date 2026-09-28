#!/usr/bin/env python3
"""Build a self-contained LaTeX explanation from saved DSC parameter draws.

No estimation is performed. The template is research/parameter_explanation.template.tex.
Run with --parameter-file PATH to pin a saved run; otherwise use the newest saved file.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.stats import gaussian_kde, invgamma


ROOT = Path(__file__).resolve().parents[1]


def tex(value):
    return str(value).replace("_", r"\_").replace("%", r"\%").replace("&", r"\&")


def number(value):
    if value == 0:
        return "0"
    mantissa, exponent = f"{value:.3e}".split("e")
    return rf"{mantissa}\times10^{{{int(exponent)}}}"


def coordinates(x, y):
    return " ".join(f"({a:.6g},{b:.6g})" for a, b in zip(x, y))


def density_plot(title, values, prior):
    mean = prior.mean()
    normalized = values / mean
    upper = max(prior.ppf(0.999) / mean, normalized.max() * 1.05)
    lower = min(prior.ppf(0.0001) / mean, normalized.min() * 0.8)
    x = np.linspace(lower, upper, 150)
    density = gaussian_kde(np.log(normalized))(np.log(x)) / x
    return rf"""\begin{{tikzpicture}}
\begin{{axis}}[width=.475\linewidth,height=4.85cm,title={{{title}}},
 title style={{font=\small}},xlabel={{Parameter / prior mean}},ylabel={{Density}},
 xmin=0,xmax={upper:.6g},ymin=0,axis lines=left,grid=major,
 grid style={{gray!15}},tick label style={{font=\scriptsize}},
 label style={{font=\scriptsize}},scaled ticks=false]
\addplot[prior,line width=1.2pt,no marks] coordinates {{{coordinates(x, mean * prior.pdf(x * mean))}}};
\addplot[posterior,line width=1.2pt,no marks] coordinates {{{coordinates(x, density)}}};
\end{{axis}}\end{{tikzpicture}}"""


def trace_plot(title, values, prior_mean, ids, iterations):
    plots = []
    for chain, color in zip(np.unique(ids), ["chainone", "chaintwo", "chainthree", "chainfour"]):
        keep = ids == chain
        plots.append(rf"\addplot[{color},line width=.45pt,no marks] coordinates "
                     + "{" + coordinates(iterations[keep], values[keep] / prior_mean) + "};")
    return rf"""\begin{{tikzpicture}}
\begin{{axis}}[width=.475\linewidth,height=4.7cm,title={{{title}}},
 title style={{font=\small}},xlabel={{Completed iteration}},ylabel={{Parameter / prior mean}},
 axis lines=left,grid=major,grid style={{gray!15}},
 tick label style={{font=\scriptsize}},label style={{font=\scriptsize}},scaled ticks=false]
{chr(10).join(plots)}
\end{{axis}}\end{{tikzpicture}}"""


def panel_grid(panels):
    return "\n\\par\\smallskip\n".join(
        panels[i] + r"\hfill" + panels[i + 1] for i in range(0, len(panels), 2))


def make_record(family, label, values, mean, sd):
    quantiles = np.quantile(values, [0.05, 0.5, 0.95])
    return dict(family=family, parameter=label, prior_mean=float(mean), prior_sd=float(sd),
                draw_mean=float(np.mean(values)), draw_q05=float(quantiles[0]),
                draw_median=float(quantiles[1]), draw_q95=float(quantiles[2]))


def appendix(records):
    pages = []
    groups = [
        ("Mean-evolution variances and volatility-evolution variances", records[:28]),
        ("Correlation-coordinate evolution variances: 1--31", records[28:59]),
        ("Correlation-coordinate evolution variances: 32--62", records[59:90]),
        ("Correlation-coordinate evolution variances: 63--91", records[90:119]),
        ("Mean-evolution covariances: 1--31", records[119:150]),
        ("Mean-evolution covariances: 32--62", records[150:181]),
        ("Mean-evolution covariances: 63--91", records[181:210]),
    ]
    for title, rows in groups:
        data = []
        for row in rows:
            data.append(tex(row["parameter"]) + " & " + " & ".join(
                "$" + number(row[key]) + "$" for key in
                ["prior_mean", "prior_sd", "draw_median", "draw_q05", "draw_q95"]) + r" \\")
        pages.append(r"\clearpage\section*{" + title + "}\n" + r"""
\small All rows use the same saved run. Prior moments are analytical.
The last three columns are empirical retained-draw summaries; convergence is unverified.
Intervals are marginal 5th--95th percentiles, not simultaneous intervals.
Ticker labels omit only the common ``Index'' suffix.
\par\medskip
\begingroup\fontsize{7.8}{10.5}\selectfont\setlength{\tabcolsep}{3pt}
\renewcommand{\arraystretch}{1.45}
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}lrrrrr@{}}
\toprule
Parameter & Prior mean & Prior SD & Draw median & Draw 5\% & Draw 95\%\\
\midrule
""" + "\n".join(data) + r"""
\bottomrule
\end{tabular*}\endgroup
\par\medskip\footnotesize
V(i,j) describes covariance of changes in expected returns, in squared percentage-point units.
h(j) denotes $\sigma_{h,j}^{2}$; r(i,j) denotes $\sigma_{r,k}^{2}$ for the matrix-log coordinate
associated with that pair. Neither r(i,j) nor V(i,j) is a return correlation.
""")
    return "\n".join(pages)


def build(parameter_file, output):
    params = loadmat(parameter_file, simplify_cells=True)["parameters"]
    run = Path(params["run_dir"])
    summary_path = run / "pilot_summary.json"
    manifest_path = run.parent / "run_manifest.json"
    summary = json.loads(summary_path.read_text())
    manifest = json.loads(manifest_path.read_text())
    priors = params["priors"]
    draws = params["parameter_draws"]
    labels = [name.removesuffix(" Index") for name in manifest["tickers"]]
    m, q, count = len(labels), len(params["pair_i"]), int(params["retained_draws"])
    assert m == 14 and q == 91
    assert params["model"] == "joint_dsc_p0"
    assert manifest["source_hash"] == params["source_hash"] == priors["source_hash"]
    assert summary["T"] == params["training_returns"].shape[0]
    assert summary["saved_draws"] == count
    assert not bool(params["convergence_established"]), "Update prose for a validated run."
    assert summary["status"] == "time_limit"
    assert params["convergence"]["status"] == "not_checked"
    assert priors["ig_shape"] == 10 and priors["nuB"] == 104
    assert params["estimation_config"]["kB"] == 0.01
    assert priors["h0_scale"] == priors["r0_scale"] == 10
    assert summary["burnin"] == 20 and params["estimation_config"]["thin"] == 1
    assert params["calibration_weeks"] == 104
    expected_pairs = [(i, j) for j in range(m) for i in range(j + 1, m)]
    pairs = list(zip(params["pair_i"].astype(int) - 1, params["pair_j"].astype(int) - 1))
    assert pairs == expected_pairs
    V, h, r = draws["V"], draws["sig2h"], draws["sig2r"]
    ids, iterations = draws["chain_id"].astype(int), draws["iteration"].astype(int)
    assert V.shape == (m, m, count) and h.shape == (count, m) and r.shape == (count, q)
    assert np.all(np.isfinite(V)) and np.all(np.isfinite(h)) and np.all(np.isfinite(r))
    assert np.all(h > 0) and np.all(r > 0)
    assert np.allclose(V, V.swapaxes(0, 1), rtol=0, atol=1e-12)
    assert np.min(np.linalg.eigvalsh(np.moveaxis(V, 2, 0))) > 0
    assert np.all(iterations > summary["burnin"])
    for chain, retained in zip(params["chain_ids"], params["retained_draws_per_chain"]):
        keep = ids == chain
        assert np.sum(keep) == retained and np.all(np.diff(iterations[keep]) > 0)
    for key, values, axis in [("V", V, 2), ("sig2h", h, 0), ("sig2r", r, 0)]:
        np.testing.assert_allclose(np.mean(values, axis=axis), params["parameter_estimate"][key], rtol=1e-10)

    shape, nu, scale = float(priors["ig_shape"]), float(priors["nuB"]), priors["V0B"]
    hp = invgamma(shape, scale=float(priors["ig_scale_h"]))
    rp = invgamma(shape, scale=float(priors["ig_scale_r"]))
    vp = [invgamma((nu - m + 1) / 2, scale=scale[j, j] / 2) for j in range(m)]
    records = [make_record("V_diagonal", f"V({label},{label})", V[j, j], vp[j].mean(), vp[j].std())
               for j, label in enumerate(labels)]
    records += [make_record("sig2h", f"h({label})", h[:, j], hp.mean(), hp.std())
                for j, label in enumerate(labels)]
    records += [make_record("sig2r", f"r({labels[i]},{labels[j]})", r[:, k], rp.mean(), rp.std())
                for k, (i, j) in enumerate(pairs)]
    for i, j in pairs:
        mean = scale[i, j] / (nu - m - 1)
        variance = ((nu - m + 1) * scale[i, j] ** 2 +
                    (nu - m - 1) * scale[i, i] * scale[j, j]) / (
                        (nu - m) * (nu - m - 1) ** 2 * (nu - m - 3))
        records.append(make_record("V_offdiagonal", f"V({labels[i]},{labels[j]})", V[i, j], mean, np.sqrt(variance)))
    assert len(records) == 210 and len({row["parameter"] for row in records}) == 210

    aud, eur, sp, sx = [labels.index(label) for label in ["AUDUSD", "EURUSD", "SPXT", "SX5T"]]
    fxpair, eqpair = pairs.index((eur, aud)), pairs.index((sx, sp))
    selected = [
        (r"$\sigma^2_h$: AUDUSD", h[:, aud], hp),
        (r"$\sigma^2_h$: SPXT", h[:, sp], hp),
        (r"$\sigma^2_r$: EURUSD / AUDUSD", r[:, fxpair], rp),
        (r"$\sigma^2_r$: SX5T / SPXT", r[:, eqpair], rp),
        (r"$V_{jj}$: AUDUSD", V[aud, aud], vp[aud]),
        (r"$V_{jj}$: SPXT", V[sp, sp], vp[sp]),
    ]
    initial_rows = []
    for j, label in enumerate(labels):
        initial_rows.append(f"{label} & {priors['Bbar'][j]:.4f} & "
                            f"{2*np.sqrt(priors['VBbar'][j,j]):.4f} & {priors['mh0'][j]:.4f}" + r" \\")
    calibration_rows = []
    for row in priors["calibration"]:
        calibration_rows.append(f"{int(row['window'])} & ${number(row['sig2h_median'])}$ & "
                                f"${number(row['sig2r_median'])}$" + r" \\")
    table_stats = "\n".join(
        tex(label) + " & $" + number(prior.mean()) + "$ & $" + number(np.mean(values)) +
        "$ & " + f"{np.mean(values)/prior.mean():.2f}" + r" \\"
        for label, values, prior in [
            ("Volatility: AUDUSD", h[:, aud], hp), ("Volatility: SPXT", h[:, sp], hp),
            ("Coordinate: EURUSD / AUDUSD", r[:, fxpair], rp),
            ("Coordinate: SX5T / SPXT", r[:, eqpair], rp)])
    replacements = {
        "DRAWCOUNT": str(count), "COUNTS": ", ".join(map(str, params["retained_draws_per_chain"].astype(int))),
        "MODELSTART": summary["first_date"], "MODELEND": summary["last_date"], "T": str(summary["T"]),
        "CALSTART": summary["calibration_start"], "CALEND": summary["calibration_end"],
        "ESTIMATEDAT": tex(params["estimated_at"]), "SOURCEHASH": params["source_hash"],
        "PARAMETERFILE": parameter_file.name, "RUNDIR": run.parent.name,
        "HMEAN": number(hp.mean()), "RMEAN": number(rp.mean()),
        "HSCALE": number(float(priors["ig_scale_h"])), "RSCALE": number(float(priors["ig_scale_r"])),
        "CALIBRATION": "\n".join(calibration_rows), "INITIALROWS": "\n".join(initial_rows),
        "DENSITIES": panel_grid([density_plot(*entry) for entry in selected]),
        "TRACES": panel_grid([trace_plot(title, values, prior.mean(), ids, iterations)
                                for title, values, prior in selected[:4]]),
        "COMPARISONS": table_stats, "APPENDIX": appendix(records),
        "POSTSHAPE": str(int(shape + summary["T"] / 2)),
        "POSTDF": str(int(nu + summary["T"] - 1)),
        "HMINRATIO": f"{np.min(h.mean(axis=0))/hp.mean():.2f}",
        "HMAXRATIO": f"{np.max(h.mean(axis=0))/hp.mean():.2f}",
        "RMINRATIO": f"{np.min(r.mean(axis=0))/rp.mean():.2f}",
        "RMAXRATIO": f"{np.max(r.mean(axis=0))/rp.mean():.2f}",
    }
    source = (ROOT / "research/parameter_explanation.template.tex").read_text()
    for key, value in replacements.items():
        assert "@@" + key + "@@" in source, key
        source = source.replace("@@" + key + "@@", value)
    assert "@@" not in source
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(source)
    generated = ROOT / "research/generated/parameter_explanation"
    generated.mkdir(parents=True, exist_ok=True)
    with (generated / "parameter_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    evidence = {
        "parameter_file": str(parameter_file), "parameter_file_sha256": hashlib.sha256(parameter_file.read_bytes()).hexdigest(),
        "source_hash": params["source_hash"], "run_dir": str(run), "retained_draws": count,
        "chain_counts": params["retained_draws_per_chain"].astype(int).tolist(),
        "convergence_established": False, "fixed_parameters_summarized": len(records),
        "positive_definite_V_draws": count, "generated_source": str(output),
        "summary_scope": "All retained post-burn-in draws pooled; chain IDs preserved for traces.",
        "prior_moments": "Analytical inverse-gamma and inverse-Wishart moments.",
        "density_method": "Gaussian KDE in log(parameter / prior mean), transformed back with Jacobian.",
        "sources": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in [summary_path, manifest_path, ROOT / "toolbox/dsc_sample.m",
                              ROOT / "toolbox/dsc_calibrate_priors.m", ROOT / "toolbox/dsc_prepare_inference.m"]},
    }
    (generated / "manifest.json").write_text(json.dumps(evidence, indent=2) + "\n")
    print(json.dumps(evidence, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameter-file", type=Path)
    parser.add_argument("--output", type=Path, default=ROOT / "research/parameter_priors_posteriors.tex")
    args = parser.parse_args()
    parameter_file = args.parameter_file
    if parameter_file is None:
        parameter_file = max((ROOT / "outputs/weekly_research/parameters").glob("dsc_parameters_*.mat"))
    build(parameter_file.resolve(), args.output.resolve())
