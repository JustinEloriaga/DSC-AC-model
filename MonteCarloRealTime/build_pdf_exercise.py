#!/usr/bin/env python3
"""Insert checked saved Monte Carlo results into the standalone DSC LaTeX note."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]
BEGIN = "% BEGIN MONTE_CARLO_VALUES"
END = "% END MONTE_CARLO_VALUES"
LABELS = [r"\nu_V", r"\Psi_{11}", r"\Psi_{12}", r"\Psi_{22}",
          r"a_h", r"b_{h1}", r"b_{h2}", r"a_r", r"b_r"]
KEYS = ["nu_V", "Psi_V11", "Psi_V12", "Psi_V22", "a_h", "b_h1", "b_h2", "a_r", "b_r"]


def csv_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def matrix(values, digits=6):
    values = np.asarray(values)
    return r"\begin{pmatrix}" + r"\\".join(
        "&".join(f"{v:.{digits}f}" for v in row) for row in values
    ) + r"\end{pmatrix}"


def scientific(value):
    mantissa, exponent = f"{value:.2e}".split("e")
    return rf"{mantissa}\times10^{{{int(exponent)}}}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()
    folder = args.run_dir.resolve()
    summary = json.loads((folder / "summary.json").read_text())
    prior = json.loads((folder / "prior.json").read_text())
    old = loadmat(folder / "posterior_T3.mat", simplify_cells=True)["previous"]
    saved = loadmat(folder / "posterior_T4.mat", simplify_cells=True)
    batch, sequential = saved["batch"], saved["sequential"]
    np.testing.assert_array_equal(batch["hyper"], sequential["hyper"])
    np.testing.assert_allclose(batch["weights"], sequential["weights"], atol=1e-14, rtol=1e-12)
    for field in ["B", "h", "r", "V", "sig2h", "sig2r"]:
        np.testing.assert_array_equal(batch["particles"][field], sequential["particles"][field])
    hyper_rows = csv_rows(folder / "posterior_hyperparameters.csv")
    np.testing.assert_array_equal([r["Hyperparameter"] for r in hyper_rows], KEYS)
    for column, posterior in [("After_1_3", old), ("Batch_1_4", batch), ("Sequential_3_then_4", sequential)]:
        np.testing.assert_allclose(
            [float(row[column]) for row in hyper_rows],
            posterior["weights"] @ posterior["hyper"], atol=1e-13, rtol=1e-13,
        )

    # Check the selected component against its latent path, independently of the CSV.
    index = 0
    path = batch["particles"]
    B, h, r = path["B"][index].T, path["h"][index].T, path["r"][index]
    for T, posterior in [(3, old), (4, batch)]:
        delta = np.diff(B[:T], axis=0)
        psi = np.asarray(prior["Psi"]) + delta.T @ delta
        bh = prior["bh"] + ((h[0] - prior["mh"])**2 / (prior["ch"] + 1)
                            + (np.diff(h[:T], axis=0)**2).sum(axis=0)) / 2
        br = prior["br"] + ((r[0] - prior["mr"])**2 / (prior["cr"] + 1)
                            + (np.diff(r[:T])**2).sum()) / 2
        expected = [prior["nu"] + T - 1, psi[0, 0], psi[0, 1], psi[1, 1],
                    prior["a"] + T / 2, *bh, prior["a"] + T / 2, br]
        np.testing.assert_allclose(posterior["hyper"][index], expected, atol=1e-14, rtol=1e-13)

    macros = []

    def define(name, body):
        macros.append("\\newcommand{\\" + name + "}{" + body + "}")

    define("MCMeanPrior", matrix(np.asarray(prior["mB"])[:, None], 2))
    define("MCCovPrior", matrix(prior["CB"], 2))
    define("MCPsiPrior", matrix(prior["Psi"], 2))
    define("MCHMeanPrior", matrix(np.asarray(prior["mh"])[:, None], 9))
    define("MCRMeanPrior", f"{prior['mr']:.9f}")
    for name, value in [("MCNuPrior", prior["nu"]), ("MCAprior", prior["a"]),
                        ("MCBhPrior", prior["bh"]), ("MCBrPrior", prior["br"]),
                        ("MCChPrior", prior["ch"]), ("MCCrPrior", prior["cr"]),
                        ("MCInitialHScale", prior["ch"] + 1), ("MCInitialRScale", prior["cr"] + 1),
                        ("MCNuThree", prior["nu"] + 2), ("MCNuFour", prior["nu"] + 3),
                        ("MCShapeThree", prior["a"] + 1.5), ("MCShapeFour", prior["a"] + 2),
                        ("MCParticles", summary["particles"]), ("MCRepetitions", summary["repetitions"])]:
        define(name, f"{value:g}")
    define("MCObservations", "\n".join(
        f"{t} & {y[0]:.9f} & {y[1]:.9f}" + r" \\" for t, y in enumerate(summary["observations"], 1)
    ))
    define("MCComponentPaths", "\n".join(
        str(t + 1) + " & " + " & ".join(f"{v:.6f}" for v in [*B[t], *h[t], r[t]]) + r" \\"
        for t in range(4)
    ))
    old_hyper = old["hyper"][index]
    full_hyper = batch["hyper"][index]
    define("MCComponentHyperRows", "\n".join(
        f"${label}$ & {full_hyper[k]:.9f} & {old_hyper[k]:.9f} & {sequential['hyper'][index, k]:.9f}" + r" \\"
        for k, label in enumerate(LABELS)
    ))
    define("MCWeightedHyperRows", "\n".join(
        f"${label}$ & {float(row['Batch_1_4']):.9f} & {float(row['After_1_3']):.9f} & {float(row['Sequential_3_then_4']):.9f}" + r" \\"
        for label, row in zip(LABELS, hyper_rows)
    ))
    define("MCComponentBrThree", f"{old_hyper[8]:.12f}")
    define("MCComponentBrFour", f"{full_hyper[8]:.12f}")
    define("MCComponentBrIncrement", f"{(r[3] - r[2])**2 / 2:.12f}")
    define("MCComponentWeightThree", scientific(old["weights"][index]))
    define("MCComponentWeightFour", scientific(batch["weights"][index]))
    define("MCMaxWeightError", scientific(summary["max_weight_difference"]))
    define("MCMaxMeanError", scientific(summary["max_parameter_mean_difference"]))
    define("MCRunName", folder.name.replace("_", r"\_"))
    text = (ROOT / "research/sequential_bayes_dsc_n2.tex").read_text()
    if text.count(BEGIN) != 1 or text.count(END) != 1:
        raise ValueError("Expected exactly one generated-value region in the LaTeX source")
    before, remainder = text.split(BEGIN)
    _, after = remainder.split(END)
    generated = BEGIN + "\n% Generated from saved artifacts; rebuild with build_pdf_exercise.py.\n"
    generated += "\n".join(macros) + "\n" + END
    (ROOT / "research/sequential_bayes_dsc_n2.tex").write_text(before + generated + after)
    print(json.dumps({"status": "passed", "run": str(folder), "particles": len(batch["weights"]),
                      "selected_component": index + 1, "exact_component_hyperparameter_match": True}, indent=2))


if __name__ == "__main__":
    main()
