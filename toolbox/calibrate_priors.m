function calibrate_priors(projroot, include_inflation)
% Empirical calibration of sig2h / sig2r prior means from MONTHLY data.
%
% Method:
%   1. Aggregate daily MARS returns to monthly (matches estimation_RANDOMCORR).
%   2. For each rolling window of `win` months, compute
%        h_t = log(diag(Var(returns_window_t)))
%        r_t = vech_lower(matrix-log(Corr(returns_window_t)))
%   3. Take first differences; their sample variance estimates the random-walk
%      innovation variance.
%   4. Write calibrated_priors.mat (read by tvsvar_modified_msv2_gam2_gen).
%
% Usage: calibrate_priors            % defaults: pwd, include_inflation=true
%        calibrate_priors(projroot)
%        calibrate_priors(projroot, include_inflation)

if nargin < 1 || isempty(projroot);          projroot          = pwd;  end
if nargin < 2 || isempty(include_inflation); include_inflation = true; end

addpath(fullfile(projroot, 'core'));
addpath(fullfile(projroot, 'toolbox'));

windows = [12 24 36];   % rolling windows in months; recommendation uses middle

%% Load daily returns, aggregate to monthly
tmp = readtable(fullfile(projroot, 'data.csv'));
dates_raw = datetime(tmp.obs_date);
if include_inflation
    data_raw = [tmp.mars_equities_portfolio, ...
                tmp.mars_bonds_portfolio, ...
                tmp.mars_commodities_portfolio, ...
                tmp.mars_inflation_portfolio];
    var_names = {'eq','bd','cm','infl'};
else
    data_raw = [tmp.mars_equities_portfolio, ...
                tmp.mars_bonds_portfolio, ...
                tmp.mars_commodities_portfolio];
    var_names = {'eq','bd','cm'};
end
valid_idx   = all(~isnan(data_raw), 2);
data_daily  = data_raw(valid_idx, :);
dates_daily = dates_raw(valid_idx);

grp = year(dates_daily)*100 + month(dates_daily);
[~, ~, ic] = unique(grp);
n_periods  = max(ic);
m          = size(data_daily, 2);
data       = zeros(n_periods, m);
for pp = 1:n_periods
    data(pp, :) = sum(data_daily(ic == pp, :), 1);
end
data = data * 100;            % match estimation_RANDOMCORR scaling
T    = size(data, 1);
n_r  = m*(m-1)/2;

fprintf('Loaded %d monthly observations, m=%d variables, n_r=%d pairs\n', T, m, n_r);

%% Pair labels matching veclAtoC's column-major tril(...,-1) ordering
pair_labels = cell(n_r, 1);
k = 0;
for j = 1:m-1
    for i = j+1:m
        k = k + 1;
        pair_labels{k} = sprintf('%s-%s', var_names{i}, var_names{j});
    end
end

%% Loop over windows
results = struct();
for w_idx = 1:numel(windows)
    win = windows(w_idx);
    if win >= T; warning('window %d >= T=%d, skipping', win, T); continue; end

    log_var_t = nan(T, m);
    r_t       = nan(T, n_r);

    for t = win:T
        Xw      = data(t-win+1:t, :);
        Vw      = var(Xw);
        Cw      = corr(Xw);
        log_var_t(t, :) = log(Vw);

        [Q, L] = eig((Cw + Cw')/2, 'vector');
        L      = max(real(L), 1e-10);          % defensive: enforce PD
        A      = Q * diag(log(L)) * Q';
        A      = real((A + A')/2);
        r_t(t, :) = A(tril(true(m), -1));      % column-major lower-tri
    end

    dh = diff(log_var_t);
    dh = dh(all(~isnan(dh), 2), :);
    dr = diff(r_t);
    dr = dr(all(~isnan(dr), 2), :);

    sig2h_emp = var(dh, 0, 1);
    sig2r_emp = var(dr, 0, 1);

    results(w_idx).win = win;
    results(w_idx).sig2h_emp = sig2h_emp;
    results(w_idx).sig2r_emp = sig2r_emp;
end

%% Report
for w_idx = 1:numel(windows)
    win = windows(w_idx);
    sig2h_emp = results(w_idx).sig2h_emp;
    sig2r_emp = results(w_idx).sig2r_emp;

    fprintf('\n=================== WINDOW = %d months ===================\n', win);
    fprintf('sig2h (Var of diff log-var) per variable:\n');
    for i = 1:m
        fprintf('  sig2h_%-4s = %.4g\n', var_names{i}, sig2h_emp(i));
    end
    fprintf('  mean across variables = %.4g\n', mean(sig2h_emp));
    fprintf('  median across vars    = %.4g\n', median(sig2h_emp));

    fprintf('\nsig2r (Var of diff matrix-log-corr) per pair:\n');
    for i = 1:n_r
        fprintf('  sig2r_{%s} = %.4g\n', pair_labels{i}, sig2r_emp(i));
    end
    fprintf('  mean across pairs = %.4g\n', mean(sig2r_emp));
    fprintf('  median across pairs   = %.4g\n', median(sig2r_emp));
end

%% Choose recommendation: use middle window's MEDIAN (robust)
choose_idx = 2;
sig2h_prior_mean = median(results(choose_idx).sig2h_emp);
sig2r_prior_mean = median(results(choose_idx).sig2r_emp);
chosen_window    = results(choose_idx).win;

fprintf('\n=================== RECOMMENDATION ========================\n');
fprintf('Using window = %d, median across series:\n', chosen_window);
fprintf('  sig2h_prior_mean = %.4g\n', sig2h_prior_mean);
fprintf('  sig2r_prior_mean = %.4g\n', sig2r_prior_mean);

fprintf('\nWith T0 = 10 pseudo-observations, set inside MCMC function:\n');
fprintf('  prior_sig2h.T0 = ones(m,1)   * 10;\n');
fprintf('  prior_sig2h.V0 = ones(m,1)   * %.4g;     %% mean = V0/(T0-1) = %.4g\n', ...
        sig2h_prior_mean*9, sig2h_prior_mean);
fprintf('  prior_sig2r.T0 = ones(n_r,1) * 10;\n');
fprintf('  prior_sig2r.V0 = ones(n_r,1) * %.4g;     %% mean = %.4g\n', ...
        sig2r_prior_mean*9, sig2r_prior_mean);
fprintf('\nAlso initialize:\n');
fprintf('  sig2h_old = %.4g * ones(1, m);\n', sig2h_prior_mean);
fprintf('  sig2r_old = %.4g * ones(1, n_r);\n', sig2r_prior_mean);

% Save inside toolbox/ so it sits on the MATLAB path and is auto-loaded
% by tvsvar_modified_msv2_gam2_gen via which('calibrated_priors.mat').
out_path = fullfile(projroot, 'toolbox', 'calibrated_priors.mat');
save(out_path, ...
     'sig2h_prior_mean', 'sig2r_prior_mean', ...
     'results', 'chosen_window', 'include_inflation', 'var_names', 'pair_labels');
fprintf('\nSaved %s\n', out_path);
end
