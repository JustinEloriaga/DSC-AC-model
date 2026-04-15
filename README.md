# DSC-AC Model for MARS Portfolio Returns

TVP-VAR with Dynamic Stochastic Correlation (DSC-AC) applied to MARS portfolio return data. The model is the random coefficient VAR from Arias, Rubio-Ramirez, and Shin (2023, *Journal of Econometrics*), adapted for financial portfolio returns.

## Model

The DSC-AC model is a **Time-Varying Parameter VAR with Multivariate Stochastic Volatility** (TVP-VAR-MSV):

- **Random coefficients**: VAR coefficients B(t) follow a random walk: B(t) = B(t-1) + e(t)
- **Stochastic volatility**: Log-volatilities h(t) follow random walks with Cogley-Sargent priors
- **Dynamic correlations**: Hansen's parameterization of the correlation matrix, with each element following a random walk
- **Estimation**: MCMC via Elliptical Slice Sampling (generalized to work with any number of variables)

## Data

**`data_JRR.csv`** contains daily portfolio returns:

| Column | Description |
|--------|-------------|
| `obs_date` | Trading date |
| `mars_equities_portfolio` | Daily log-return, equities |
| `mars_bonds_portfolio` | Daily log-return, bonds |
| `mars_commodities_portfolio` | Daily log-return, commodities |
| `mars_inflation_portfolio` | Daily log-return, inflation (available from ~2006 only) |

The 3-variable model uses equities, bonds, and commodities (1984-05-04 to 2026-02-27, ~10,700 daily obs). Daily returns are aggregated to **monthly** frequency (~500 obs) before estimation, since the TVP-VAR-MSV model builds T x T covariance matrices internally.

## Quick Start

```matlab
cd fcst_dcc_msv_var
main_v3c_msv2_gam2_mars3
```

This runs the 3-variable model (equities, bonds, commodities) with default settings. Output is saved to `pred_v4_msv2_gam2_mars3/`.

## Configuration

Edit `get_default_info_v3_msv2_gam2_mars3.m` to change:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ndraws` | 1000 | MCMC draws (use 50000+ for production) |
| `nburn` | 100 | Burn-in draws |
| `nthin` | 5 | Thinning interval |
| `p` | 2 | VAR lags |
| `hmax` | 8 | Forecast horizon (months) |
| `eval_T0` | `'1990-01-31'` | Start of evaluation period |
| `eval_T1` | `'2025-12-31'` | End of evaluation period |
| `primiceri` | 3 | Model type (3 = Hansen's MSV2, the JE paper model) |
| `T0` | 40 | Training sample size for prior |

## File Structure

```
DSC-AC-model/
+-- data_JRR.csv                          # MARS portfolio data
+-- README.md
+-- fcst_dcc_msv_var/
    +-- main_v3c_msv2_gam2_mars3.m        # Main runfile (3-var)
    +-- main_estimation_and_forecast_v3c_msv_gam2_mars3.m  # Core estimation & forecast
    +-- get_default_info_v3_msv2_gam2_mars3.m              # Hyperparameters & config
    +-- compute_crps2.m                   # CRPS scoring
    +-- fcst_Primiceri_v00/               # Estimation internals
    |   +-- tvsvar_modified_msv2_gam2_gen.m   # MCMC sampler (generalized)
    |   +-- fcst_var_primiceri_msv2_gen.m      # Forecast evaluation (generalized)
    |   +-- kfilter.m, kback.m                # Kalman filter/smoother
    |   +-- lag.m, trimr.m, make_varXY.m      # VAR utilities
    |   +-- olsblock.m, ols1.m                # OLS helpers
    |   +-- mvnrnd_modified.m                 # MVN random draws
    |   +-- LogAbsDet.m, lognormpdf.m         # Math utilities
    |   +-- tvsvar_modified.m                 # Original Primiceri (alt. model)
    |   +-- fcst_var_primiceri.m              # Original Primiceri forecast
    |   +-- step_sv*.m, tria*.m, cols.m       # Supporting functions
    +-- toolbox_msv/                      # MSV toolbox
    |   +-- loglike_yt_given_ht_Pt.m          # Likelihood (volatility)
    |   +-- slice_sampling_v02.m              # Elliptical slice sampler
    |   +-- LogAbsDet.m                       # Log absolute determinant
    +-- toolbox_msv_corr/                 # Dynamic correlation toolbox
        +-- loglike_yt_given_rt_ht.m          # Likelihood (correlation)
        +-- veclAtoC.m                        # Vector to correlation matrix
        +-- log_mvnpdf.m                      # Log MVN PDF
```

## References

Arias, J. E., Rubio-Ramirez, J. F., & Shin, M. (2023). Macroeconomic forecasting and variable ordering in multivariate stochastic volatility models. *Journal of Econometrics*, 235(2), 1054-1086.
