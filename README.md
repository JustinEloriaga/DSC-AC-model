# RANDOMCORR — TVP-VAR with Dynamic Stochastic Correlation

TVP-VAR with Dynamic Stochastic Correlation applied to MARS portfolio return data. The model is the random coefficient VAR from Arias, Rubio-Ramirez, and Shin (2023, *Journal of Econometrics*), adapted for financial portfolio returns.

## Model

The model is a **Time-Varying Parameter VAR with Multivariate Stochastic Volatility** (TVP-VAR-MSV):

- **Random coefficients**: VAR coefficients B(t) follow a random walk
- **Stochastic volatility**: Log-volatilities h(t) follow random walks with Cogley-Sargent priors
- **Dynamic correlations**: Hansen's parameterization — each correlation element follows a random walk
- **Estimation**: MCMC via Elliptical Slice Sampling

## Quick Start

```matlab
main
```

That is all. Everything is configured at the top of `main.m`.

## Data

**`data_JRR.csv`** — daily portfolio log-returns for three MARS strategies:

| Column | Description |
|--------|-------------|
| `obs_date` | Trading date |
| `mars_equities_portfolio` | Equities |
| `mars_bonds_portfolio` | Bonds |
| `mars_commodities_portfolio` | Commodities |

Daily returns are summed within each period to produce the model frequency. The end date is always set automatically to the last available observation in the file.

## Configuration

All settings live at the top of **`main.m`** — nothing else needs to be edited.

### Frequency

```matlab
info.freq = 'monthly';   % ~500 obs — default, fast
info.freq = 'weekly';    % ~2200 obs — slower; consider increasing info.T0
```

Switching frequency changes:
- How daily returns are aggregated (monthly sum vs. weekly sum)
- The output filename: `RANDOMCORR_monthly_YYYY-MM-DD.mat` or `RANDOMCORR_weekly_YYYY-MM-DD.mat`
- The output PDF: `correlations_monthly.pdf` or `correlations_weekly.pdf`

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `info.freq` | `'monthly'` | Data frequency: `'monthly'` or `'weekly'` |
| `info.eval_T1` | `'9999-12-31'` | End of estimation window (sentinel = use all data) |
| `info.p` | `0` | VAR lags |
| `info.nex` | `1` | Include constant |
| `info.T0` | `40` | Training sample length (periods, used for prior) |
| `info.kB` | `0.01` | Prior scaling for VAR coefficients |
| `info.ndraws` | `50` | MCMC draws (use 10000+ for production) |
| `info.nburn` | `10` | Burn-in draws |
| `info.nthin` | `1` | Thinning interval |
| `info.nreport` | `2` | Report progress every N draws |

## Outputs

Each run produces two files in the project root:

| File | Contents |
|------|----------|
| `RANDOMCORR_monthly_YYYY-MM-DD.mat` | MCMC draws: `r.B`, `r.V`, `r.r`, `r.P`, `r.h`, `r.sig2h`, `r.sig2r` |
| `correlations_monthly.pdf` | Time-series plots of the three pairwise correlations (mean ± 1 std) |

For weekly runs the filenames are `RANDOMCORR_weekly_...` and `correlations_weekly.pdf`. Old `.mat` files are deleted automatically before each run.

## File Structure

```
DSC-AC-model/
├── main.m                          ← run this
├── data_JRR.csv
├── core/                           ← MCMC sampler and VAR utilities
│   ├── tvsvar_modified_msv2_gam2_gen.m   MCMC sampler
│   ├── kfilter.m, kback.m                Kalman filter / smoother
│   ├── make_varXY.m, lag.m               VAR helpers
│   ├── olsblock.m                        OLS
│   ├── mvnrnd_modified.m                 MVN draws
│   ├── LogAbsDet.m, lognormpdf.m         Math utilities
├── toolbox/                        ← likelihoods, slice sampler, estimation
│   ├── estimation_RANDOMCORR.m           Data loading, aggregation, MCMC call, save
│   ├── loglike_yt_given_ht_Pt.m          Log-likelihood for volatility
│   ├── loglike_yt_given_rt_ht.m          Log-likelihood for correlation
│   ├── slice_sampling_v02.m              Elliptical slice sampler
│   ├── veclAtoC.m                        Correlation vector → matrix (Hansen)
│   ├── log_mvnpdf.m                      Log MVN PDF
└── archive_original/               ← original code before refactor (for comparison)
```

## Reference

Arias, J. E., Rubio-Ramirez, J. F., & Shin, M. (2023). Macroeconomic forecasting and variable ordering in multivariate stochastic volatility models. *Journal of Econometrics*, 235(2), 1054–1086.
