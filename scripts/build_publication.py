#!/usr/bin/env python3
"""Rebuild the LaTeX report and Beamer PDF from saved results; never runs MCMC.

Requires Python 3, numpy, pandas, matplotlib, and a TeX distribution with latexmk.
Generated tables and vector evidence charts share one numerical source of truth.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / "research/generated/.matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(__file__).resolve().parents[1] / "research/generated/.cache"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
NAVY, TEAL = "#17324D", "#148C8A"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.labelcolor": NAVY, "text.color": NAVY,
                     "axes.edgecolor": "#A9B3BF", "pdf.fonttype": 42,
                     "savefig.bbox": "tight"})


def tex(value):
    return "".join({"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%",
                    "$": r"\$", "#": r"\#", "_": r"\_", "{": r"\{",
                    "}": r"\}", "~": r"\textasciitilde{}",
                    "^": r"\textasciicircum{}"}.get(c, c) for c in str(value))


def short(ticker):
    return str(ticker).removesuffix(" Index")


def number(x, places=2):
    return "--" if pd.isna(x) else f"{float(x):.{places}f}"


def pvalue(row):
    bound = str(getattr(row, "PValueBound", ""))
    for prefix, symbol in [("<=", r"\leq"), (">=", r"\geq"), ("<", "<"), (">", ">")]:
        if bound.startswith(prefix):
            return "$" + symbol + " " + bound[len(prefix):] + "$"
    return number(row.PValue, 3)


def yes(series):
    return series.astype(str).str.lower().isin(["1", "true", "1.0"])


def table(headers, rows, align=None, small=True, long=False):
    align = align or ("l" + "r" * (len(headers) - 1))
    env = "longtable" if long else "tabular"
    result = "\\par\\medskip\\noindent\n" + ("{\\footnotesize\n" if long else ("{\\small\n" if small else "")) + f"\\begin{{{env}}}{{@{{}}{align}@{{}}}}\n\\toprule\n"
    head = " & ".join(headers) + r"\\" + "\n\\midrule\n"
    result += head
    if long:
        result += "\\endfirsthead\n\\toprule\n" + head + "\\endhead\n"
    result += "\n".join(" & ".join(str(c) for c in row) + r"\\" for row in rows)
    result += f"\n\\bottomrule\n\\end{{{env}}}\n" + ("}\n" if small or long else "") + "\\par\\medskip\n"
    return result


def read_csv(folder, name, required):
    df = pd.read_csv(folder / name)
    absent = set(required) - set(df.columns)
    if absent:
        raise ValueError(f"{name}: missing columns {sorted(absent)}")
    return df


def build(args):
    folder = args.data.resolve()
    generated = ROOT / "research/generated"
    figures = generated / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    def write(name, value):
        (generated / name).write_text(value, encoding="utf-8")

    summary = json.loads((folder / "data_summary.json").read_text())
    if summary.get("closure_policy") != "asof":
        raise ValueError("The baseline publication uses asof marks. Supply baseline outputs and describe mask results as a sensitivity.")
    returns = read_csv(folder, "weekly_returns.csv", ["Date"])
    levels = read_csv(folder, "weekly_levels.csv", ["Date"])
    quality = read_csv(folder, "weekly_quality.csv", ["AgeDays", "Ticker"])
    full = read_csv(folder, "full_sample_correlations.csv", ["PairIndex", "TickerI", "TickerJ", "Correlation", "N"])
    rolling = read_csv(folder, "rolling_correlations.csv", ["Date", "Window", "PairIndex", "TickerI", "TickerJ", "Correlation", "N"])
    tests = read_csv(folder, "stationarity.csv", ["Ticker", "Transform", "Test", "Lags", "Statistic", "PValue", "Reject5pct", "N"])
    diagnostic = pd.read_csv(folder / "diagnostics.csv")
    tickers = list(returns.columns[1:])
    if len(tickers) != 14 or len(full) != 91 or any("USDCNH" in x for x in tickers):
        raise ValueError("Expected 14 non-USDCNH series and 91 distinct correlations")
    dates = pd.to_datetime(returns.Date)
    if (len(returns) != summary["n_returns"] or len(levels) != summary["n_levels"] or
            returns.iloc[0, 0] != summary["first_return"] or returns.iloc[-1, 0] != summary["last_return"]):
        raise ValueError("Weekly sample dimensions or dates disagree with data_summary.json")
    if not dates.is_unique or not dates.is_monotonic_increasing or not dates.dt.dayofweek.eq(4).all():
        raise ValueError("Weekly dates must be unique increasing Fridays")
    if len(dates) > 1 and not dates.diff().dropna().eq(pd.Timedelta(days=7)).all():
        raise ValueError("Weekly calendar has missing dates")
    if list(levels.columns[1:]) != tickers or len(levels) != len(returns) + 1:
        raise ValueError("Level/return panel dimensions or ticker labels disagree")
    reconstructed = 100 * np.diff(np.log(levels[tickers].to_numpy(dtype=float)), axis=0)
    if not np.allclose(reconstructed, returns[tickers].to_numpy(dtype=float), rtol=1e-9, atol=1e-10, equal_nan=True):
        raise ValueError("Published returns differ from percentage log changes of the saved levels")
    sourcefile = Path(summary.get("source_file", ""))
    if sourcefile.is_file() and hashlib.sha256(sourcefile.read_bytes()).hexdigest() != summary["source_hash"]:
        raise ValueError("Source file changed after data analysis")
    rolling["Date"] = pd.to_datetime(rolling.Date)
    if rolling.duplicated(["Date", "Window", "PairIndex"]).any():
        raise ValueError("Duplicate rolling-correlation output keys")
    if set(rolling.Window.unique()) != {52, 104}:
        raise ValueError("Expected exactly 52- and 104-week rolling windows")
    for window in [52, 104]:
        sample = rolling[rolling.Window.eq(window)]
        counts = sample.groupby("Date").size()
        expected_dates = pd.DatetimeIndex(dates.iloc[window-1:])
        if not counts.index.equals(expected_dates) or not counts.eq(91).all():
            raise ValueError(f"Incomplete date/pair coverage in {window}-week correlations")
        if not sample.N.between(0, window).all():
            raise ValueError(f"Invalid sample counts in {window}-week correlations")
    full = full.sort_values("PairIndex")
    # Verify the canonical column-major lower triangle, independently of labels.
    expected = [(tickers[i], tickers[j]) for j in range(14) for i in range(j + 1, 14)]
    if list(zip(full.TickerI, full.TickerJ)) != expected:
        raise ValueError("Pair ordering differs from MATLAB column-major lower triangle")
    first, last = returns.Date.iloc[0], returns.Date.iloc[-1]
    macros = {"FirstReturn": first, "LastReturn": last, "NumReturns": f"{len(returns):,}",
              "NumLevels": f"{len(levels):,}", "InputHash": summary["source_hash"]}
    write("metrics.tex", "\n".join(f"\\newcommand{{\\{key}}}{{{tex(value)}}}" for key, value in macros.items()) + "\n")

    corr = np.eye(14)
    for row in full.itertuples():
        i, j = tickers.index(row.TickerI), tickers.index(row.TickerJ)
        corr[i, j] = corr[j, i] = row.Correlation
    fig, ax = plt.subplots(figsize=(8.2, 7))
    im = ax.imshow(corr, vmin=-1, vmax=1, cmap="RdBu_r")
    labels = [short(t) for t in tickers]
    ax.set(xticks=np.arange(14), yticks=np.arange(14), xticklabels=labels, yticklabels=labels)
    ax.tick_params(axis="both", length=0, labelsize=8)
    plt.setp(ax.get_xticklabels(), rotation=60, ha="right", rotation_mode="anchor")
    for i in range(14):
        for j in range(14):
            ax.text(j, i, f"{corr[i,j]:.2f}", ha="center", va="center", fontsize=6.3,
                    color="white" if abs(corr[i,j]) > .64 else NAVY)
    fig.colorbar(im, ax=ax, shrink=.75, label="Pearson correlation")
    ax.set_title(f"{first} to {last}  |  pairwise available weekly returns", fontsize=10, pad=14)
    fig.tight_layout()
    fig.savefig(figures / "full_correlation.pdf")
    plt.close(fig)

    ages = quality.AgeDays.value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(6.3, 3.7))
    ax.bar(ages.index.astype(int), ages.values, color=[TEAL if x else NAVY for x in ages.index], width=.7)
    ax.set(xlabel="Calendar days before Friday", ylabel="Endpoint observations", xticks=range(8))
    ax.set_yscale("log")
    ax.set_ylim(.7, max(ages) * 4)
    for x, count in ages.items():
        ax.text(x, count * 1.2, f"{int(count):,}", ha="center", fontsize=8)
    ax.grid(axis="y", alpha=.18)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(figures / "endpoint_age.pdf")
    plt.close(fig)

    def pair(a, b):
        chosen = full[(full.TickerI.map(short).eq(a) & full.TickerJ.map(short).eq(b)) |
                      (full.TickerI.map(short).eq(b) & full.TickerJ.map(short).eq(a))]
        if len(chosen) != 1:
            raise ValueError(f"No unique requested pair: {a}/{b}")
        return chosen.iloc[0]

    pairs = [pair("SPXT", "SX5T"), pair("SPXT", "NKYTR"), pair("SPXT", "LUATTRUU"),
             pair("SX5T", "LEATTREU"), pair("EURUSD", "GBPUSD"), pair("AUDUSD", "USDJPY")]

    def pairplot(rows, filename, ncols=2, height=2.4, numbered=False, width=10):
        nrows = (len(rows) + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(width, height * nrows), squeeze=False)
        for ax, row in zip(axes.flat, rows):
            data = rolling[rolling.PairIndex.eq(row.PairIndex)]
            for window, color in [(104, TEAL), (52, NAVY)]:
                d = data[data.Window.eq(window)].sort_values("Date")
                ax.plot(d.Date, d.Correlation, color=color, lw=.85, label=f"{window} weeks")
            title = f"{short(row.TickerI)} / {short(row.TickerJ)}"
            ax.set_title((f"{int(row.PairIndex):02d}. " if numbered else "") + title, fontsize=9, loc="left")
            ax.set(ylim=(-1, 1), yticks=[-1, 0, 1])
            ax.axhline(0, color="#A9B3BF", lw=.6)
            ax.grid(axis="y", alpha=.12)
            ax.xaxis.set_major_locator(mdates.YearLocator(5))
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
            ax.tick_params(labelsize=8)
        for ax in list(axes.flat)[len(rows):]:
            ax.axis("off")
        handles, labs = axes.flat[0].get_legend_handles_labels()
        fig.legend(handles[::-1], labs[::-1], loc="lower center", ncol=2, frameon=False, fontsize=9)
        fig.tight_layout(rect=[0, .05 / nrows, 1, 1], h_pad=1.3)
        fig.savefig(figures / filename)
        plt.close(fig)

    pairplot(pairs, "selected_pairs.pdf", height=2.4)
    for name, rows in [("equity", pairs[:2]), ("cross", pairs[2:4]), ("currency", pairs[4:])]:
        pairplot(rows, name + "_pairs.pdf", height=3.2)
    appendix = []
    all_rows = list(full.itertuples())
    start = 0
    for page, size in enumerate([9, 9, 9] + [8] * 8):
        rows = all_rows[start:start + size]
        filename = f"all_pairs_{page+1:02d}.pdf"
        pairplot(rows, filename, ncols=3, height=3.2, numbered=True, width=8.2)
        appendix += ["\\clearpage\n" if page else "", f"\\subsection*{{Pairs {start+1}--{start+size}}}\n",
                     f"\\includegraphics[width=\\linewidth]{{generated/figures/{filename}}}\n"]
        start += size
    write("all_pairs.tex", "".join(appendix))

    rawmiss = int(summary["raw_missing_cells"])
    carried = int(summary["carried_endpoint_cells"])
    closure = int(summary["closure_return_cells"])
    qrows = [["Raw daily rows", f"{summary['raw_rows']:,}"],
             ["Raw daily missing cells", f"{rawmiss:,} ({100*summary['raw_missing_fraction']:.2f}\\%)"],
             ["Weekly values not from a Friday close", f"{carried:,} ({100*summary['carried_endpoint_fraction']:.2f}\\%)"],
             ["Fridays with at least one non-Friday value", f"{summary['fridays_with_carry']:,}"],
             ["Oldest value used in a Friday row", f"{summary['max_age_days']} days"],
             ["Flagged NKYTR returns around 2019 closure", str(closure)]]
    write("data_quality.tex", table(["Quality measure", "Count / rate"], qrows))

    descrows = []
    for ticker in tickers:
        x = pd.to_numeric(returns[ticker], errors="coerce").dropna()
        descrows.append([tex(short(ticker)), str(len(x)), number(x.mean(), 3), number(x.std(), 3),
                         number(x.quantile(.01), 2), number(x.quantile(.99), 2)])
    write("descriptive_table.tex", table(["Ticker", "$N$", "Mean (\\%)", "SD (\\%)", "1st pct.", "99th pct."], descrows))

    r52 = rolling[(rolling.Window.eq(52)) & rolling.Correlation.notna()]
    latest = r52[r52.Date.eq(r52.Date.max())].set_index("PairIndex")
    ranges = r52.groupby("PairIndex").Correlation.agg(["min", "max"])
    ranges["range"] = ranges["max"] - ranges["min"]
    greatest = int(ranges["range"].idxmax())
    widest = full[full.PairIndex.eq(greatest)].iloc[0]
    strong = full.iloc[np.abs(full.Correlation.to_numpy()).argmax()]
    stronglabel = f"{short(strong.TickerI)} / {short(strong.TickerJ)}"
    widelabel = f"{short(widest.TickerI)} / {short(widest.TickerJ)}"
    executive = (
        f"The weekly panel contains \\textbf{{{len(returns):,} dates, 14 variables, and 91 pairs}}. "
        "Each weekly row is a Friday row: for each series, we use its last published level on or before that Friday. "
        "This keeps a local holiday in one market from deleting the whole week for every market. "
        "Japan's NKYTR series has one unusual 2019 holiday week; the affected returns are flagged so they can be checked separately.\n\n"
    )
    write("executive.tex", executive)
    windowrows = []
    for row in pairs:
        current104 = rolling[(rolling.Window.eq(104)) & rolling.PairIndex.eq(row.PairIndex)].sort_values("Date").iloc[-1]
        windowrows.append([tex(f"{short(row.TickerI)} / {short(row.TickerJ)}"), number(row.Correlation),
                           number(latest.loc[row.PairIndex, "Correlation"]), number(current104.Correlation),
                           number(ranges.loc[row.PairIndex, "min"]), number(ranges.loc[row.PairIndex, "max"])])
    write("rolling_findings.tex", table(["Pair", "Full sample", "Latest 52", "Latest 104", "52 min", "52 max"], windowrows) +
          f"Latest estimates end {last}. Minimum and maximum refer to the entire available 52-week path.\n")
    write("slides_windows.tex", table(["Pair", "Latest 52", "Latest 104"], [[r[0], r[2], r[3]] for r in windowrows], small=False) +
          f"\\par\\medskip\\small End date: {last}. Pearson correlations on fixed calendar windows.\n")

    transform = tests.Transform.astype(str).str.lower()
    rt = tests[transform.str.contains("return")].copy()
    if rt.empty:
        raise ValueError(f"Cannot identify return diagnostics among transforms {tests.Transform.unique()}")
    adf = rt[rt.Test.astype(str).str.upper().eq("ADF")]
    kpss = rt[rt.Test.astype(str).str.upper().eq("KPSS")]
    adfn = int(yes(adf.Reject5pct).sum())
    k13 = kpss[kpss.Lags.eq(13)]
    flagged = [short(t) for t in k13[yes(k13.Reject5pct)].Ticker]
    flaggedtext = ", ".join(tex(t) for t in flagged) if flagged else "none"
    stationarity_findings = (f"Table~\\ref{{tab:return-stationarity}} reports the tests on percentage weekly log returns. "
                            f"ADF rejects a unit root for {adfn} of {len(adf)} return series at 5\\%. ")
    kpss_groups = {}
    for bandwidth in [4, 13, 26]:
        kb = kpss[kpss.Lags.eq(bandwidth)]
        names = ", ".join(tex(short(t)) for t in kb[yes(kb.Reject5pct)].Ticker) or "none"
        kpss_groups.setdefault(names, []).append(str(bandwidth))
    for names, bandwidths in kpss_groups.items():
        stationarity_findings += (f"At bandwidths {', '.join(bandwidths)}, KPSS rejects stationarity "
                                 f"around a constant mean for: {names}. ")
    stationarity_findings += "These are return results; the KPSS null's name does not describe the input transformation.\n"
    la = tests[tests.Transform.eq("log_level") & tests.Test.eq("ADF")]
    drift = la[la.Model.eq("ARD")]
    trend = la[la.Model.eq("TS")]
    trend_names = ", ".join(tex(short(t)) for t in trend[yes(trend.Reject5pct)].Ticker) or "none"
    lk = tests[tests.Transform.eq("log_level") & tests.Test.eq("KPSS") & tests.Model.eq("level")]
    level_findings = (f"Table~\\ref{{tab:level-stationarity}} reports the tests on weekly log levels. "
                      f"The constant-only ADF rejects a unit root for {int(yes(drift.Reject5pct).sum())} of {len(drift)} series. ")
    level_groups = {}
    for bandwidth in [4, 13, 26]:
        kb = lk[lk.Lags.eq(bandwidth)]
        counts = (int(yes(kb.Reject5pct).sum()), len(kb))
        level_groups.setdefault(counts, []).append(str(bandwidth))
    for (rejections, total), bandwidths in level_groups.items():
        level_findings += (f"KPSS rejects constant-mean stationarity for {rejections} of {total} "
                           f"log-level series at bandwidths {', '.join(bandwidths)}. ")
    level_findings += (f"With a deterministic trend, ADF rejects for {trend_names}. "
                       "Thus the broad evidence against stationarity concerns log levels. Failure to reject an ADF unit root alone does not prove that a unit root exists.\n")
    write("level_stationarity_findings.tex", level_findings)
    write("stationarity_findings.tex", stationarity_findings)
    levelrows = []
    for ticker in tickers:
        a = drift[drift.Ticker.eq(ticker)].iloc[0]
        at = trend[trend.Ticker.eq(ticker)].iloc[0]
        row = [tex(short(ticker)), pvalue(a), pvalue(at)]
        for lag in [4, 13, 26]:
            row.append(pvalue(lk[lk.Ticker.eq(ticker) & lk.Lags.eq(lag)].iloc[0]))
        levelrows.append(row)
    write("level_stationarity_table.tex", table(
        ["Ticker", "ADF const.", "ADF trend", "KPSS(4)", "KPSS(13)", "KPSS(26)"],
        levelrows, align="lrrrrr"))
    strows = []
    for ticker in tickers:
        a = adf[adf.Ticker.eq(ticker)].iloc[0]
        row = [tex(short(ticker)), str(int(a.Lags)), number(a.Statistic, 3), pvalue(a)]
        for lag in [4, 13, 26]:
            k = kpss[kpss.Ticker.eq(ticker) & kpss.Lags.eq(lag)].iloc[0]
            row.append(pvalue(k))
        strows.append(row)
    write("stationarity_table.tex", table(["Ticker", "ADF lag", "ADF stat.", "ADF $p$", "KPSS(4)", "KPSS(13)", "KPSS(26)"], strows, align="lrrrrrr"))
    write("return_adf_table.tex", table(["Ticker", "ADF lag", "ADF statistic", "$p$-value"],
                                       [row[:4] for row in strows], align="lrrr"))
    write("return_adf_conclusion.tex",
          f"ADF rejects a unit root for {adfn} of {len(adf)} return series at the 5\\% level. "
          f"The largest ADF $p$-value shown in the table is {adf.PValue.max():.3f} (subject to the tabulated bounds in the table).\n")
    fullst = []
    for row in tests.itertuples():
        bound = str(getattr(row, "PValueBound", ""))
        transform_label = "return" if "return" in row.Transform else row.Transform.replace("log_level", "log level")
        fullst.append([tex(short(row.Ticker)), tex(transform_label), tex(row.Test), tex(getattr(row, "Model", "")),
                       str(int(row.Lags)), number(row.Statistic, 3), pvalue(row),
                       "Y" if str(row.Reject5pct).lower() in ["true", "1", "1.0"] else "N", str(int(row.N))])
    write("stationarity_full.tex", table(["Ticker", "Transform", "Test", "Model", "Lag", "Statistic", "$p$", "Reject", "$N$"], fullst,
                                        align="llllrrrrr", long=True))

    if {"Test", "Reject5pct"}.issubset(diagnostic.columns):
        dep = diagnostic.assign(_reject=yes(diagnostic.Reject5pct)).groupby("Test")._reject.agg(["sum", "count"])
        deprows = [[tex(name), f"{int(row['sum'])} / {int(row['count'])}"] for name, row in dep.iterrows()]
    else:
        raise ValueError("diagnostics.csv must include Test and Reject5pct")
    raw_ljung_box = diagnostic[diagnostic.Test.eq("LjungBoxReturns")]
    rejected_raw = raw_ljung_box[yes(raw_ljung_box.Reject5pct)]
    rejected_names = ", ".join(tex(short(t)) for t in rejected_raw.Ticker) or "none"
    arch = dep.loc["ARCH"]
    squared = dep.loc["LjungBoxSquaredReturns"]
    write("dependence_findings.tex",
          "\\par\\textbf{Do earlier returns help explain this week's return?} "
          "The Ljung--Box test checks whether a series is correlated with its own returns over the previous 13 weeks, roughly three months. "
          f"It finds evidence of a relationship in {len(rejected_raw)} of {len(raw_ljung_box)} series: {rejected_names}. "
          "This suggests that past returns may contain some information about subsequent returns. "
          "The test alone does not tell us how large that relationship is, whether it reflects continuation or reversal, or whether it is useful for forecasting.\n\n"
          "\\par\\textbf{Does the size of market moves vary in persistent spells?} "
          "The ARCH test and the Ljung--Box test on squared, demeaned returns focus on the size of moves, regardless of whether prices rise or fall. "
          f"They find evidence of dependence in {int(arch['sum'])} of {int(arch['count'])} and "
          f"{int(squared['sum'])} of {int(squared['count'])} series, respectively, using 13 weekly lags. "
          "This is consistent with volatility clustering: turbulent weeks tend to occur in spells, as do calmer weeks. "
          "It supports allowing the model's estimate of risk to change over time. It does not, by itself, show that correlations between different markets change.\n\n"
          "\\par\\textbf{What does this mean for our model?} "
          "We allow volatility to vary and will assess whether adding a simple dependence on last week's return improves the mean equation. "
          "These checks describe the observed data before model fitting. After estimation, we must repeat them on the part of returns the model leaves unexplained (the innovations) "
          "to see whether patterns remain. Dependence across weeks is not the same as a unit root and does not by itself call for differencing returns again.\n\n"
          "\\par{\\footnotesize Reading the statistics: a rejection means the test found evidence against its no-dependence assumption at the chosen 5\\% threshold. "
          "Testing many series creates more opportunities for chance flags; these counts have not been adjusted for that. "
          "Source: saved \\texttt{diagnostics.csv}, percentage weekly log returns, "
          f"{first}--{last}; individual statistics and $p$-values remain in that file.}}\n")
    write("slides_diagnostics.tex", table(["Dependence diagnostic", "Rejections / tests"], deprows, small=False) +
          "\\par\\medskip\\small These are nominal 5\\% decisions without multiplicity correction. The full report records stationarity specifications and all 91 correlation paths.\n")
    write("slides_stationarity.tex", f"\\begin{{itemize}}\\item Return ADF rejects a unit root for {adfn} of 14 variables."
          f"\\item KPSS at 13 weeks flags: {flaggedtext}."
          "\\item Lag and bandwidth sensitivity matter for mean specification."
          "\\item Serial-dependence tests motivate an AR-mean sensitivity before production.\\end{itemize}\n")
    write("slides_findings.tex", "\\begin{itemize}"
          f"\\item {len(returns):,} weekly observations support 91 pairwise paths."
          f"\\item Largest absolute full-sample correlation: {tex(stronglabel)} ({strong.Correlation:.2f})."
          f"\\item The 52-week {tex(widelabel)} path ranges from {ranges.loc[greatest,'min']:.2f} to {ranges.loc[greatest,'max']:.2f}."
          f"\\item ADF rejects a unit root for {adfn} return series. KPSS sensitivity qualifies the stationarity conclusion."
          "\\end{itemize}\n")

    make_pilot(args, write, folder)
    posterior = make_bayesian_paths(args, write, figures, summary)
    manifest = {"data_dir": str(folder), "pilot": str(args.pilot.resolve()) if args.pilot else None,
                "posterior_run": str(args.posterior_run.resolve()) if args.posterior_run else None,
                "n_returns": len(returns), "n_pairs": len(full), "pair_order": expected,
                "source_hash": summary["source_hash"],
                "csv_hashes": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(folder.glob("*.csv"))}}
    if args.pilot:
        manifest["pilot_sha256"] = hashlib.sha256(args.pilot.read_bytes()).hexdigest()
        resource_path = args.pilot.parent / "resource_observations.json"
        if resource_path.exists():
            manifest["resource_observations_sha256"] = hashlib.sha256(resource_path.read_bytes()).hexdigest()
        for label, evidence_path in pilot_evidence_paths(args.pilot, folder).items():
            if evidence_path.exists():
                manifest[label + "_sha256"] = hashlib.sha256(evidence_path.read_bytes()).hexdigest()
                manifest[label + "_path"] = str(evidence_path.resolve())
    if posterior:
        manifest["posterior"] = posterior
    write("publication_manifest.json", json.dumps(manifest, indent=2) + "\n")
    if not args.no_compile:
        output = ROOT / "output/pdf"
        output.mkdir(parents=True, exist_ok=True)
        latexmk = shutil.which("latexmk") or "/Library/TeX/texbin/latexmk"
        env = dict(os.environ)
        env["PATH"] = str(Path(latexmk).parent) + os.pathsep + env.get("PATH", "")
        stems = ["weekly_correlation_report"] if args.report_only else ["weekly_correlation_report", "weekly_correlation_slides"]
        for stem in stems:
            subprocess.run([latexmk, "-pdf", "-interaction=nonstopmode", "-halt-on-error", "-file-line-error",
                            "-outdir=../output/pdf", stem + ".tex"], cwd=ROOT / "research", env=env, check=True)
            print(output / (stem + ".pdf"))


def make_bayesian_paths(args, write, figures, data_summary):
    """Copy and describe correlation paths from the explicitly selected run."""
    if args.posterior_run is None:
        write("bayesian_paths.tex", "")
        return None
    run = args.posterior_run.resolve()
    summary_path = run / "pilot_summary.json"
    manifest_path = run / "figures/bayes_correlation_paths_manifest.json"
    run_manifest_path = run / "run_manifest.json"
    required = [summary_path, manifest_path, run_manifest_path,
                run / "figures/bayes_correlation_paths_selected.pdf"]
    required += [run / f"figures/bayes_correlation_paths_all_{page:02d}.pdf" for page in range(1, 12)]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Posterior run is missing generated path files: " + ", ".join(missing))
    if args.pilot and args.pilot.resolve() != summary_path:
        raise ValueError("--pilot and --posterior-run must refer to the same sampler run")
    sampler = json.loads(summary_path.read_text())
    paths = json.loads(manifest_path.read_text())
    run_manifest = json.loads(run_manifest_path.read_text())
    if run_manifest.get("source_hash") != data_summary["source_hash"]:
        raise ValueError("Posterior paths and publication data have different source hashes")
    draws = int(paths["retained_draws"])
    warmup = int(paths["warmup_removed"])
    if draws < 2 or draws != int(sampler.get("saved_draws", -1)):
        raise ValueError("Posterior-path draw count disagrees with the sampler summary")
    if paths.get("color_scale") != [-1, 1] or "Continuous RdBu_r" not in paths.get("color_rule", ""):
        raise ValueError("Posterior paths must use the Figure 1 continuous correlation scale")
    selected_source = run / "figures/bayes_correlation_paths_selected.pdf"
    selected_target = figures / "bayesian_paths_selected.pdf"
    shutil.copy2(selected_source, selected_target)
    pages = []
    copied = [selected_target]
    for page in range(1, 12):
        source = run / f"figures/bayes_correlation_paths_all_{page:02d}.pdf"
        target = figures / f"bayesian_paths_all_{page:02d}.pdf"
        shutil.copy2(source, target)
        copied.append(target)
        prefix = "" if page == 1 else "\\clearpage\n"
        height = ".80" if page == 1 else ".90"
        pages.append(prefix +
                     f"\\includegraphics[width=\\linewidth,height={height}\\textheight,keepaspectratio]"
                     f"{{generated/figures/{target.name}}}\n")
    convergence = bool(sampler.get("convergence_established", False))
    status_text = ("Convergence diagnostics passed." if convergence else
                   "Convergence has not been established, so the bands are preliminary.")
    write("bayesian_paths.tex",
          "\\section{Bayesian correlation paths}\n"
          f"The selected sampler run retained {draws:,} draws after {warmup:,} warm-up iterations. "
          "The dark line is the posterior median of the actual correlation $P_{ij,t}$ and the gray envelope "
          "contains the pointwise 16th and 84th percentiles. Vertical colors use the same continuous scale "
          "as Figure~1: blue for negative values, white near zero, and red for positive values, with intensity "
          "set by the posterior median. " + status_text + "\n"
          "\\begin{figure}[p]\\centering\n"
          "\\includegraphics[width=\\linewidth,height=.80\\textheight,keepaspectratio]"
          "{generated/figures/bayesian_paths_selected.pdf}\n"
          "\\caption{Selected smoothed conditional innovation-correlation paths. The intervals are pointwise "
          "68\\% posterior bands from the saved run; they are not simultaneous bands.}\n"
          "\\end{figure}\n"
          "\\clearpage\n"
          "\\subsection{All 91 Bayesian correlation paths}\n"
          "The following pages use the same vertical scale, color scale, and posterior summaries for every pair.\n" +
          "".join(pages) + "\\clearpage\n")
    return {
        "run_dir": str(run),
        "retained_draws": draws,
        "warmup_iterations": warmup,
        "convergence_established": convergence,
        "summary_sha256": hashlib.sha256(summary_path.read_bytes()).hexdigest(),
        "paths_manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "figure_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in copied},
    }


def pilot_evidence_paths(pilot_path, folder):
    run = pilot_path.parent
    acceptance = run / "acceptance_tests.json"
    return {
        "warmup_validation": run / "validation.json",
        "matlab_test_results": acceptance if acceptance.exists() else folder.parent / "matlab_test_results.json",
        "performance_comparison": run / "performance_comparison.json",
        "kernel_benchmark": run / "kernel_benchmark.json",
        "final_likelihood_validation": run / "final_likelihood_validation.json",
    }


def make_pilot(args, write, folder):
    """Only report recorded pilot values; absent results stay explicitly unavailable."""
    priorfile = folder / "prior_summary.json"
    prior = json.loads(priorfile.read_text()) if priorfile.exists() else {}
    sensitivity = prior.get("calibration", [])
    if sensitivity:
        rows = [[str(item["window"]), number(item["sig2h_median"], 7), number(item["sig2r_median"], 7)] for item in sensitivity]
        write("prior_sensitivity.tex", table(["Window (weeks)", "Log-variance prior mean", "Correlation-coordinate prior mean"], rows) +
              "These are prior means for innovation variances, not return correlations. A longer calibration window lowers the scale substantially.\n")
    else:
        write("prior_sensitivity.tex", "No saved prior-sensitivity results were supplied to this build.\n")
    if args.pilot is None:
        write("pilot_executive.tex", "No pilot summary was supplied, so the production cost remains unassessed.\n")
        write("pilot.tex", "No pilot summary was supplied to this build. Bayesian timing, memory and posterior results are unavailable.\n")
        write("slides_pilot.tex", "No pilot summary supplied. Numerical results and computational timing are unavailable in this build.\n")
        write("production.tex", "A production cost recommendation requires a completed, measured pilot.\n")
        write("slides_production.tex", "The production decision requires a measured pilot and a review of model diagnostics.\n")
        return
    pilot = json.loads(args.pilot.read_text())
    data_summary = json.loads((folder / "data_summary.json").read_text())
    for pilot_key, data_key in [("T", "n_returns"), ("m", "n_series"), ("pairs", "n_pairs"),
                                ("first_date", "first_return"), ("last_date", "last_return")]:
        if pilot_key in pilot and pilot[pilot_key] != data_summary[data_key]:
            raise ValueError(f"Pilot {pilot_key} disagrees with the published panel")
    evidence_paths = pilot_evidence_paths(args.pilot, folder)
    evidence = {label: json.loads(path.read_text()) for label, path in evidence_paths.items() if path.exists()}
    performance = evidence.get("performance_comparison", {})
    benchmark = evidence.get("kernel_benchmark", {})
    implementation = pilot.get("implementation", {})
    accelerated = pilot.get("correlation_backend", implementation.get("correlation_backend")) == "mex"
    if performance:
        if Path(performance["new_summary"]).resolve() != args.pilot.resolve():
            raise ValueError("Performance comparison belongs to a different pilot")
        if (performance["new_completed_sweeps"] != pilot.get("completed_iterations") or
                performance["numerical_validation"]["source_hash"] != data_summary["source_hash"]):
            raise ValueError("Performance comparison disagrees with the pilot or data source")
    if benchmark:
        if benchmark.get("source_hash") != data_summary["source_hash"] or benchmark.get("implementation") != implementation:
            raise ValueError("Kernel benchmark belongs to a different input or sampler implementation")
    write("pilot_raw.json", json.dumps(pilot, indent=2) + "\n")
    # Explicit aliases support the sampler's versioned JSON and older pilot summaries.
    def get(*keys, default=None):
        for key in keys:
            if key in pilot:
                return pilot[key]
        return default
    n = get("completed_iterations", "iterations_completed", "completedIterations")
    elapsed = get("elapsed_seconds", "wall_seconds", "elapsedSeconds", "session_seconds", "total_seconds")
    per = get("mean_iteration_seconds", "seconds_per_iteration", "meanIterationSeconds", "mean_sweep_seconds")
    if n is None or elapsed is None:
        raise ValueError("Pilot summary needs completed_iterations and elapsed_seconds")
    if per is None:
        per = float(elapsed) / n if n else None
    reason = get("stop_reason", "reason", "status", default="bounded pilot")
    status = get("status", default="bounded_pilot")
    compact_status = {"time_limit": "Time cap reached", "iteration_limit": "Iteration cap reached",
                      "failed": "Numerical failure", "bounded_pilot": "Bounded pilot"}.get(status, "See recorded reason")
    rows = [["Completed iterations", str(n)], ["Retained posterior draws", str(get("saved_draws", default=0))],
            ["Elapsed wall time", f"{float(elapsed)/60:.2f} minutes"],
            ["Mean time per completed iteration", f"{float(per):.2f} seconds" if n and per is not None else "Unavailable"], ["Stop status", compact_status]]
    if get("r_proposals") is not None:
        rows.append(["Correlation / volatility proposals", f"{get('r_proposals'):,} / {get('h_proposals', default=0):,}"])
    resource_text = "Peak resident process memory was not measured."
    peak_rss = performance.get("peak_process_rss_gib")
    if peak_rss is not None:
        if not np.isfinite(float(peak_rss)) or float(peak_rss) <= 0:
            raise ValueError("Peak process RSS must be positive finite GiB")
        rows.append(["Peak MATLAB process RSS", f"{float(peak_rss):.3f} GiB"])
        resource_text = ("Peak MATLAB process RSS was measured by macOS \\texttt{time -l}; "
                         "it includes data preparation and validation. ")
    resource_path = args.pilot.parent / "resource_observations.json"
    if resource_path.exists():
        resources = json.loads(resource_path.read_text())
        observations = resources.get("observations", [])
        for observation in observations:
            rss_kib = float(observation["rss_kib"])
            if not np.isfinite(rss_kib) or rss_kib <= 0:
                raise ValueError("Resource observation needs a positive finite RSS in KiB")
            rss_gib = rss_kib / 2**20
            process_elapsed = tex(observation.get("process_elapsed", "unavailable"))
            rows.append([f"Sampled RSS at process elapsed {process_elapsed}", f"{rss_gib:.3f} GiB"])
        if observations:
            stamp = tex(observations[-1].get("observed_at_utc", "unrecorded"))
            sample_text = (
                "The resident-set observation uses macOS \\texttt{ps}, converting its KiB output to GiB. "
                f"The latest recorded observation time is {stamp}. "
                "This is point-in-time MATLAB process memory, not a measured peak. "
                "It includes the MATLAB runtime and data-preparation overhead and excludes helper processes. "
                "Process elapsed time starts with MATLAB, which can precede the sampler timer."
            )
            resource_text = resource_text + sample_text if peak_rss is not None else sample_text
    for label, keys, fmt in [("Likelihood proposals", ["proposal_count", "total_proposals", "likelihood_evaluations"], ",.0f"),
                             ("Estimated working-array memory", ["estimated_memory_mb", "working_memory_mb"], ".1f"),
                             ("Peak resident process memory (MB)", ["peak_rss_mb", "peak_memory_mb"], ".1f")]:
        value = get(*keys)
        if value is not None:
            rows.append([label, format(value, fmt)])
    storage = []
    for key, label in [("old_dense_brownian_GiB", "Original dense random-walk factors"),
                       ("old_4000_draw_h_r_P_GiB", "Original 4,000-draw h/r/P storage"),
                       ("packed_P_2000_draw_GiB", "Distinct correlations, 2,000 retained draws")]:
        if get(key) is not None:
            storage.append([label, f"{float(get(key)):.2f} GiB"])
    storagetex = table(["Analytical array/storage estimate", "Size"], storage) if storage else ""
    stage = get("stage_seconds", default={})
    stage_text = ""
    if stage and sum(stage.values()) > 0:
        correlation_seconds = float(stage.get("correlation", 0))
        correlation_fraction = 100 * correlation_seconds / sum(stage.values())
        stage_text = (
            "The slow part of the pilot is the Bayesian correlation step. "
            "At each iteration, the sampler turns proposed correlation states into weekly "
            "correlation matrices and evaluates the likelihood of the observed returns. "
            f"This step used {correlation_fraction:.1f}\\% of the recorded sampler time "
            f"({correlation_seconds:.2f} seconds). "
        )
        if not accelerated:
            stage_text += "This makes the correlation step the main target for speed improvements. "
        if status == "time_limit":
            stage_text += (f"The unfinished sweep used {float(get('unfinished_sweep_seconds', default=0)):.2f} seconds "
                           "before the time guard discarded it. ")
        stage_text += "\\par\n"
    speed_text = ""
    slide_speed = ""
    if accelerated:
        threads = get("correlation_threads")
        thread_text = f"{int(threads)} native threads" if threads is not None else "native threads"
        speed_text = (
            f"We rewrote that calculation in compiled code using {thread_text}. "
            "The statistical model is unchanged. "
        )
    if benchmark:
        speed_text += (
            f"On the same saved input, the compiled calculation was {benchmark['kernel_speedup']:.2f} times faster "
            "and matched the MATLAB likelihood to numerical precision. "
        )
    if performance:
        speed_text += (f"Across {performance['common_completed_sweeps']} matched pilot iterations, "
                       f"the full sampler was {performance['common_sweep_speedup']:.2f} times faster. ")
        if performance.get("common_proposal_counts_match"):
            speed_text += "The matched iterations used the same proposal counts, so this is a like-for-like timing comparison. "
        slide_speed = (f"\\small Whole sweeps: {performance['common_sweep_speedup']:.2f} times faster "
                       f"over {performance['common_completed_sweeps']} common sweeps")
        if benchmark:
            slide_speed += f"; fixed-input kernel: {benchmark['kernel_speedup']:.2f} times faster"
        slide_speed += ".\\par\n"
    speed_text = ""
    validation_text = ""
    if "matlab_test_results" in evidence:
        test_results = evidence["matlab_test_results"]
        validation_text += (f"The saved MATLAB test run records {test_results['passed']} passing tests, "
                            f"with {test_results['failed']} failures. ")
    if "warmup_validation" in evidence:
        validation = evidence["warmup_validation"]
        if validation.get("source_hash") != data_summary["source_hash"]:
            raise ValueError("Warm-up validation belongs to a different data source")
        validation_text += (
            f"The last warm-up state was also checked across all {validation['matrices_checked']:,} weekly correlation matrices. "
            "These checks confirm that the saved pilot is numerically consistent; they are not convergence evidence.\\par\n"
        )
    if "final_likelihood_validation" in evidence:
        likelihood = evidence["final_likelihood_validation"]
        if (likelihood.get("status") != "passed" or likelihood.get("source_hash") != data_summary["source_hash"] or
                likelihood.get("completed_iterations") != n or likelihood.get("implementation") != implementation):
            raise ValueError("Final likelihood validation disagrees with the completed pilot")
        maximum_error = max(likelihood["absolute_errors"].values())
        validation_text += ("A separate likelihood check gives the same result from MATLAB and the compiled backend, "
                            "up to numerical rounding.\\par\n")
    write("pilot.tex", table(["Measured pilot quantity", "Value"], rows,
                             align=r"p{.58\linewidth}p{.34\linewidth}") + storagetex +
          "\\textbf{Recorded stop reason:} " + tex(reason) + "\\par\n" +
          stage_text)
    slide_rows = [row for row in rows if not row[0].startswith(("Sampled RSS", "Peak MATLAB"))][:6]
    sampled = [row for row in rows if row[0].startswith("Sampled RSS")]
    if peak_rss is not None:
        slide_rows.append(["Peak MATLAB process RSS", f"{float(peak_rss):.3f} GiB"])
    elif sampled:
        slide_rows.append(["Sampled MATLAB RSS (not peak)", sampled[-1][1]])
    slide_status_note = "\\small The time guard discarded the unfinished sweep and retained the last completed checkpoint.\n" if status == "time_limit" else ""
    write("slides_pilot.tex", table(["Pilot quantity", "Recorded value"], slide_rows,
                                    align=r"p{.58\linewidth}p{.34\linewidth}") + slide_status_note + slide_speed)
    if not n or per is None or not np.isfinite(float(per)):
        write("pilot_executive.tex", "The supplied sampler run completed no full iteration, so a per-iteration runtime is unavailable.\n")
        write("production.tex", "The pilot completed no full iteration. A reliable per-iteration cost and production projection are unavailable. Investigate the recorded unfinished stage before authorizing production.\n")
        write("slides_production.tex", "The pilot completed no full iteration. The next step is to investigate the unfinished stage. A production cost projection is unavailable.\n")
        return
    hours = float(per) * 4000 / 3600
    next_step = ("Acceleration is implemented; mean/prior sensitivity and convergence assessment remain. " if accelerated else
                 "Mean/prior sensitivity and profiling or acceleration should precede a production decision. ")
    cost_scope = "Projections exclude checkpoint writes and production posterior-chunk storage. "
    write("pilot_executive.tex", f"The supplied sampler run completed {n} iterations in {float(elapsed)/60:.1f} minutes "
          f"and retained {get('saved_draws', default=0)} posterior draws. Its measured mean time was "
          f"{float(per):.2f} seconds per completed iteration.\n")
    write("production.tex", "")
    write("slides_production.tex", "\\textbf{Recommendation: hold production.}\\par\\medskip\\small " + next_step + "\\par\\medskip\n" +
          f"Illustrative 4,000-iteration chain: \\textbf{{{hours:.1f} hours}} at pilot speed. "
          f"Four serial chains: \\textbf{{{4*hours:.1f} hours}}. " + cost_scope + "No convergence guarantee.\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "outputs/weekly_research/data")
    parser.add_argument("--results", type=Path, help="Alternative: research output root containing data/")
    parser.add_argument("--pilot", type=Path, help="Saved pilot_summary.json; no estimation is launched")
    parser.add_argument("--posterior-run", type=Path,
                        help="Sampler run containing generated Bayesian correlation-path PDFs")
    parser.add_argument("--report-only", action="store_true", help="Compile the report but not the Beamer slides")
    parser.add_argument("--no-compile", action="store_true", help="Regenerate tables/charts only")
    args = parser.parse_args()
    if args.results:
        args.data = args.results / "data"
    build(args)


if __name__ == "__main__":
    main()
