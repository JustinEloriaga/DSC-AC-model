# MonteCarloRealTime

A reproducible two-series (`n=2`), four-date (`T=4`) DSC experiment comparing:

1. Original prior plus observations 1:4.
2. Original prior plus observations 1:3, followed by the saved joint posterior as the prior and observation 4 only.

Run in MATLAB, from the repository root:

```matlab
addpath('MonteCarloRealTime');
result = run_monte_carlo();             % 20,000 particles, 5 paired repetitions
disp(result.output_dir);
```

Requires MATLAB and Statistics and Machine Learning Toolbox. No production estimation is launched. The experiment uses the same distributions as `toolbox/dsc_sample.m`, with fixed synthetic hyperparameters and the exact bivariate matrix-log map `rho = tanh(r)`.

The two routes use matching random numbers and retain weights without resampling. This makes both finite-particle posterior representations equal to machine precision. Independent Monte Carlo runs target the same posterior but need not return identical numerical estimates.

The observed-data posterior is not described by a single conjugate hyperparameter vector because the paths are latent. The report gives the conditional hyperparameters of a posterior mixture, their weighted summaries, parameter posterior means/SDs, and old-state smoothing. It checks every mixture component, not just the reported averages. The generated true states are not passed to inference.

## Results

Each run creates a new folder under `results/` and preserves earlier runs:

- `REPORT.md`: interpretation and side-by-side numerical tables.
- `observations.csv`, `prior.json`, `synthetic_data.mat`: data and generation inputs.
- `posterior_T3.mat`: reusable joint date-3 posterior with weights, states, parameter draws, conditional hyperparameters, and RNG state.
- `posterior_T4.mat`: both full posterior representations.
- `posterior_hyperparameters.csv`: shapes/degrees of freedom and weighted scales.
- `example_mixture_components.csv`: the first twelve component comparisons.
- `posterior_parameters.csv`: static-parameter posterior means and SDs.
- `smoothed_states.csv`: state means/SDs and pre-update means.
- `replications.csv`, `summary.json`: numerical discrepancies and importance ESS.

Only `dsc_rt_update(previous,y_new)` performs the sequential update. It receives the old posterior and the new observation; its signature does not accept earlier observations. Each append uses the same static parameter draw associated with its old path. Conditional hyperparameters are recalculated for batch and incremented for sequential.

The existing PDF at `output/pdf/sequential_bayes_dsc_n2.pdf` includes a worked exercise
using these saved outputs. To refresh the numerical values in its standalone source,
without rerunning inference:

```sh
python3 MonteCarloRealTime/build_pdf_exercise.py MonteCarloRealTime/results/RUN_DIRECTORY
```

The builder uses NumPy/SciPy to verify the saved MATLAB posterior objects and CSV
tables, then updates only the marked generated-value block in
`research/sequential_bayes_dsc_n2.tex`. Recompile that source to refresh the PDF.

For fewer particles or another output directory:

```matlab
result = run_monte_carlo(4000, 2);
```

Run the focused regression tests:

```matlab
addpath('MonteCarloRealTime');
results = runtests('MonteCarloRealTime/test_monte_carlo.m');
assertSuccess(results);
```
