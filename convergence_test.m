% Convergence diagnostics for the TVP-VAR + Hansen MSV2 sampler.
%
% Runs a SINGLE long chain with NO thinning, then computes:
%   - Autocorrelation function (ACF) for a panel of scalar series
%   - Effective Sample Size (ESS, Geyer initial monotone sequence)
%   - Geweke z-score (first 10% vs last 50% of post-burn-in samples)
% and recommends per-chain `ndraws` and `nthin` for the multi-chain config.

clc; clear variables; close all;
workpath = pwd;
rng('default'); rng(83);

addpath([workpath, filesep, 'core']);
addpath([workpath, filesep, 'toolbox']);

%% Configuration — single long chain, no thinning
info = [];
info.workpath          = workpath;
info.datapath          = workpath;
info.freq              = 'monthly';
info.include_inflation = true;
info.eval_T1           = '9999-12-31';
info.p                 = 0;
info.nex               = 1;
info.kB                = 0.01;

info.ndraws  = 2000;   % long enough to expose autocorrelation structure
info.nburn   = 500;    % discarded inside this script as well
info.nthin   = 1;      % NO thinning -- we want raw ACF
info.nreport = 100;
info.nchains = 1;      % single chain for this diagnostic

%% Run
delete(fullfile(workpath, 'RANDOMCORR_monthly_*.mat'));
disp('Starting convergence-test chain ...');
tStart = tic;
estimation_RANDOMCORR(info);
disp(['Elapsed: ', num2str(toc(tStart),'%.1f'), ' sec']);

%% Load result
mat_files = dir(fullfile(workpath, 'RANDOMCORR_monthly_*.mat'));
[~, latest_idx] = max([mat_files.datenum]);
matfile = fullfile(workpath, mat_files(latest_idx).name);
S = load(matfile, 'r'); r = S.r;

% Strip burn-in (mat_* arrays contain all M = ndraws + nburn draws since nthin=1)
burn      = info.nburn;
mat_h     = r.h(:,:,burn+1:end);
mat_sig2h = r.sig2h(burn+1:end, :);
mat_sig2r = r.sig2r(burn+1:end, :);
mat_P     = r.P(:,:,:,burn+1:end);

[T, m, N] = size(mat_h);
n_r = m*(m-1)/2;

%% Build panel of scalar series
series_names = {};
series_data  = [];

% Innovation-variance hyperparameters
for i = 1:m
    series_names{end+1} = sprintf('sig2h_%d', i);
    series_data(:, end+1) = mat_sig2h(:, i);
end
for i = 1:n_r
    series_names{end+1} = sprintf('sig2r_%d', i);
    series_data(:, end+1) = mat_sig2r(:, i);
end

% Representative correlations P(i,j,t) at 3 t's
t_check = round([T/4, T/2, 3*T/4]);
for ipair = 1:m
    for jpair = ipair+1:m
        for tt = t_check
            series_names{end+1} = sprintf('P_%d%d_t%d', ipair, jpair, tt);
            series_data(:, end+1) = squeeze(mat_P(ipair, jpair, tt, :));
        end
    end
end

S_n = numel(series_names);
fprintf('Diagnosing %d scalar series over %d post-burn-in draws.\n', S_n, N);

%% ACF
max_lag = min(200, floor(N/4));
acf     = zeros(max_lag+1, S_n);
for s = 1:S_n
    x = series_data(:, s) - mean(series_data(:, s));
    var_x = mean(x.^2);
    if var_x < 1e-14
        acf(:, s) = 0; acf(1, s) = 1;
        continue
    end
    for L = 0:max_lag
        acf(L+1, s) = mean(x(1:N-L) .* x(L+1:N)) / var_x;
    end
end

acf_threshold = 0.1;
lag_below = zeros(S_n, 1);
for s = 1:S_n
    idx = find(abs(acf(:, s)) < acf_threshold, 1, 'first');
    if isempty(idx); lag_below(s) = max_lag; else; lag_below(s) = idx - 1; end
end

%% ESS via Geyer initial monotone (positive) sequence
ess = zeros(S_n, 1);
for s = 1:S_n
    rho = acf(2:end, s);                       % rho(k) at lag k
    pair_sum = rho(1:2:end-1) + rho(2:2:end);
    K = find(pair_sum < 0, 1, 'first');
    if isempty(K)
        rho_sum = sum(rho);
    else
        rho_sum = sum(rho(1:2*K-2));           % stop at last positive pair
    end
    ess(s) = N / (1 + 2*max(rho_sum, 0));
end

%% Geweke (first 10% vs last 50%)
geweke_z = zeros(S_n, 1);
n1 = floor(0.1 * N);
n2 = floor(0.5 * N);
for s = 1:S_n
    x1 = series_data(1:n1, s);
    x2 = series_data(end-n2+1:end, s);
    se1 = std(x1) / sqrt(n1);
    se2 = std(x2) / sqrt(n2);
    if (se1^2 + se2^2) == 0
        geweke_z(s) = 0;
    else
        geweke_z(s) = (mean(x1) - mean(x2)) / sqrt(se1^2 + se2^2);
    end
end

%% Per-series table
fprintf('\n========================= DIAGNOSTICS =========================\n');
fprintf('Post burn-in N = %d (info.ndraws = %d, info.nburn = %d)\n', N, info.ndraws, info.nburn);
fprintf('%-20s  ACF<%.2g@lag    ESS     |Geweke z|\n', 'Series', acf_threshold);
fprintf('---------------------------------------------------------------\n');
for s = 1:S_n
    fprintf('%-20s   %5d        %6.1f     %5.2f\n', ...
            series_names{s}, lag_below(s), ess(s), abs(geweke_z(s)));
end

%% Recommendations
worst_lag   = max(lag_below);
median_ess  = median(ess);
min_ess     = min(ess);
max_geweke  = max(abs(geweke_z));

% nthin chosen as the worst ACF-decay lag (each thinned sample roughly independent)
nthin_rec = max(1, worst_lag);

% Target ~400 effective draws TOTAL across multi-chain run; this is a standard
% rule-of-thumb for stable posterior summaries
target_eff_total = 400;
% Estimate ESS-per-iteration for the worst series (proxy for sampling efficiency)
ess_per_iter = min_ess / N;
% Per-chain ndraws to hit target with 4 chains
nchains_target = 4;
ndraws_per_chain_rec = ceil(target_eff_total / (ess_per_iter * nchains_target));

fprintf('\n========================= RECOMMENDATIONS =====================\n');
fprintf('Worst ACF-decay lag across all series : %d\n', worst_lag);
fprintf('Median ESS / Min ESS                  : %.0f / %.0f\n', median_ess, min_ess);
fprintf('Max |Geweke z|                        : %.2f  (>2 => non-converged)\n', max_geweke);
fprintf('Sampling efficiency (min ESS / N)     : %.3f\n', ess_per_iter);
fprintf('\n');
fprintf('Recommended info.nthin                : %d\n', nthin_rec);
fprintf('Recommended info.ndraws per chain     : %d   (for %d chains, ~%d eff. draws total)\n', ...
        ndraws_per_chain_rec, nchains_target, target_eff_total);

if max_geweke > 2
    fprintf('\n*** WARNING: Geweke z > 2 on at least one series.\n');
    fprintf('    Increase info.nburn before trusting these recommendations.\n');
end

%% Plot trace + ACF for a representative subset
nplot = min(6, S_n);
sel = round(linspace(1, S_n, nplot));
fig = figure('Color','w','Position',[100 100 1400 900]);
for idx = 1:nplot
    s = sel(idx);
    subplot(nplot, 2, 2*idx-1);
    plot(series_data(:, s));
    title(['Trace: ', series_names{s}], 'Interpreter','none');
    xlabel('Draw'); grid on;

    subplot(nplot, 2, 2*idx);
    stem(0:max_lag, acf(:, s), 'Marker','none');
    yline(acf_threshold, 'r--'); yline(-acf_threshold, 'r--');
    title(['ACF: ', series_names{s}], 'Interpreter','none');
    xlabel('Lag'); ylim([-0.3 1.05]); grid on;
end
exportgraphics(fig, fullfile(workpath, 'convergence_diagnostics.pdf'));
fprintf('\nTrace+ACF plots saved to: convergence_diagnostics.pdf\n');
