# Two-Series DSC: Batch Versus Sequential Bayes

Both routes agree to floating-point precision. This run used 20000 particles in each of 5 paired repetitions, with the same four observations.

Maximum component-hyperparameter difference: `0`. Maximum normalized-weight difference: `8.67e-18`. Maximum parameter-mean difference: `1.39e-17`.

## Data and Starting Prior

Data were simulated from the DSC prior with seed 20260929. Inference receives only the observations; the simulated latent truth is saved separately for inspection and is never treated as observed.

| Date | y1 | y2 |
|---|---:|---:|
| 1 | -0.4837400063 | -0.1248558243 |
| 2 | -0.2245600949 | -0.5459561184 |
| 3 | -0.3221661259 | -0.5848109623 |
| 4 | -0.2900639822 | -1.0977376731 |

`B1 ~ N(mB, CB)`, independently of V. `V ~ IW(Psi, nu)`.

`sig2h_j ~ IG(a, bh)` and `sig2r ~ IG(a, br)`, with density proportional to `q^(-a-1) exp(-b/q)`.

`h1_j | sig2h_j ~ N(mh_j, (ch+1)*sig2h_j)`; `r1 | sig2r ~ N(mr, (cr+1)*sig2r)`. Later states follow random walks. `rho_t = tanh(r_t)` and `Sigma_t = diag(exp(h_t/2))*R_t*diag(exp(h_t/2))`.

All given hyperparameters, identical in both routes:

```json
{
  "mB": [
    0,
    0
  ],
  "CB": [
    [
      0.25,
      0.05
    ],
    [
      0.05,
      0.25
    ]
  ],
  "Psi": [
    [
      0.15,
      0.03
    ],
    [
      0.03,
      0.2
    ]
  ],
  "nu": 8,
  "mh": [
    -0.35667494393873245,
    0.18232155679395459
  ],
  "mr": 0.20273255405408219,
  "ch": 10,
  "cr": 10,
  "a": 6,
  "bh": 0.15,
  "br": 0.05
}
```

## What Is Compared

1. **Batch:** draw complete candidate paths from the original prior and weight them by `L1*L2*L3*L4`.
2. **Sequential:** draw paths through date 3, weight by `L1*L2*L3`, save and reload that joint posterior, append date-4 states from their transition distributions, and multiply the old weights by **L4 only**. `dsc_rt_update` accepts no old observations.
3. Use matched random numbers and no intermediate resampling. Thus the two routes have the same candidate paths and the same final weights, up to floating-point arithmetic. Old path weights change, so dates 1--3 are smoothed using observation 4.

The common prior particles are a controlled comparison, not independent Monte Carlo runs. The original continuous posterior is approximated. Exact equality of the paired approximations does not prove that finite-sample Monte Carlo error is zero. Different repetition seeds illustrate that uncertainty.

## Posterior Hyperparameters

Because states are latent, the observed-data posterior has **no single IW/IG hyperparameter vector**. We represent it by a weighted mixture: for each candidate path s,

- `V | path_s, y ~ IW(Psi_s, nu_s)`;
- `sig2h_j | path_s, y ~ IG(a_h, b_hj_s)`;
- `sig2r | path_s, y ~ IG(a_r, b_r_s)`.

Every component hyperparameter is checked between routes. The table reports the shared degrees of freedom/shapes and the **posterior-weighted mean of each scale**, using repetition 1. Scale averages do not turn the mixture into one IW/IG distribution.

| Hyperparameter | Prior | After 1:3 | Batch 1:4 | Sequential 3+1 | Max component difference |
|---|---:|---:|---:|---:|---:|
| nu_V | 8.000000000000 | 10.000000000000 | 11.000000000000 | 11.000000000000 | 0 |
| Psi_V11 | 0.150000000000 | 0.203390538022 | 0.226104931188 | 0.226104931188 | 0 |
| Psi_V12 | 0.030000000000 | 0.041736647358 | 0.046961012944 | 0.046961012944 | 0 |
| Psi_V22 | 0.200000000000 | 0.275065970367 | 0.315644134807 | 0.315644134807 | 0 |
| a_h | 6.000000000000 | 7.500000000000 | 8.000000000000 | 8.000000000000 | 0 |
| b_h1 | 0.150000000000 | 0.201370500453 | 0.229821458635 | 0.229821458635 | 0 |
| b_h2 | 0.150000000000 | 0.203066867568 | 0.216380891678 | 0.216380891678 | 0 |
| a_r | 6.000000000000 | 7.500000000000 | 8.000000000000 | 8.000000000000 | 0 |
| b_r | 0.050000000000 | 0.066647998441 | 0.072056132465 | 0.072056132465 | 0 |

In particular, `nu_V: 8 -> 10 -> 11`, and `a_h = a_r: 6 -> 7.5 -> 8`. Each path-specific scale is updated by its new squared increment (or outer product for V). The h/r initial quadratic term is included.

Complete mixture components and weights are in `posterior_T4.mat`; twelve illustrative components are in `example_mixture_components.csv`. These are inferred candidate paths, not the simulated true path.

## Static Parameter Posterior Summaries

Means and standard deviations integrate each conditional IW/IG component analytically and then average over its posterior path weight. This reduces Monte Carlo noise and retains parameter uncertainty.

| Parameter | Mean after 1:3 | Batch mean | Sequential mean | Batch SD | Sequential SD |
|---|---:|---:|---:|---:|---:|
| V11 | 0.029055791146 | 0.028263116399 | 0.028263116399 | 0.020566657375 | 0.020566657375 |
| V12 | 0.005962378194 | 0.005870126618 | 0.005870126618 | 0.016634405151 | 0.016634405151 |
| V22 | 0.039295138624 | 0.039455516851 | 0.039455516851 | 0.030385850748 | 0.030385850748 |
| sig2h1 | 0.030980076993 | 0.032831636948 | 0.032831636948 | 0.018632919860 | 0.018632919860 |
| sig2h2 | 0.031241056549 | 0.030911555954 | 0.030911555954 | 0.015682525810 | 0.015682525810 |
| sig2r | 0.010253538222 | 0.010293733209 | 0.010293733209 | 0.005223681306 | 0.005223681306 |

## The Old History Also Updates

The first-date states below change after observation 4. Both routes give the same new smoothed means. All dates, including SDs, are in `smoothed_states.csv`.

| First-date state | Before y4 | Batch after y4 | Sequential after y4 |
|---|---:|---:|---:|
| B1 | -0.206780107967 | -0.230870145186 | -0.230870145186 |
| B2 | -0.192512714132 | -0.311032149600 | -0.311032149600 |
| h1 | -0.695878293170 | -0.870896371046 | -0.870896371046 |
| h2 | -0.187041590945 | -0.200373043701 | -0.200373043701 |
| r | 0.241091134504 | 0.215741108362 | 0.215741108362 |
| rho | 0.213671738368 | 0.191997113657 | 0.191997113657 |

For date 4, the `Mean_before_y4` column in the CSV is a predictive mean; for dates 1--3 it is the old smoothed mean. State labels B1/B2 and h1/h2 identify series; the Time column identifies the date.

## Monte Carlo Diagnostics

ESS here means importance-weight ESS, `1/sum(w.^2)`, not MCMC autocorrelation ESS.

| Repeat | ESS after 1:3 | ESS after 1:4 | Max hyperparameter difference | Max weight difference |
|---|---:|---:|---:|---:|
| 1 | 3600.7 | 1030.2 | 0 | 6.94e-18 |
| 2 | 3505.6 | 1010.6 | 0 | 6.07e-18 |
| 3 | 3857.7 | 1364.4 | 0 | 4.34e-18 |
| 4 | 3143.8 | 1040.7 | 0 | 8.67e-18 |
| 5 | 2731.6 | 1316.3 | 0 | 6.94e-18 |

Across-repetition SD of the estimated means (V11, V12, V22, sig2h1, sig2h2, sig2r):

`[0.000361651 0.000180675 0.000351813 0.000596614 0.000406502 0.000210164]`

This is empirical Monte Carlo variation for one fixed data set, not a posterior SD and not a study of repeated data sets. Each pair of routes shares this numerical error.

Elapsed time: 2.40 seconds. Both evidence calculations agree: direct `mean(L1*L2*L3*L4)` and old evidence times the sequential predictive density.

The final weighted posterior is the same, with no refitting of hyperpriors, no use of latent truth, and no repeated old-data likelihood in the sequential update. This is a small teaching exercise, not an alteration of the production sampler.
