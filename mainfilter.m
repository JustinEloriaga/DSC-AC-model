% TVP-VAR with Hansen's MSV2 model — MARS portfolio data
% FILTER-ONLY VARIANT (no backward smoother for B)
% 3-variable VAR: [equities, bonds, commodities]

%% Housekeeping
clc; clear variables; close all;
workpath = pwd;
rng('default');
rng(83);

addpath([workpath, filesep, 'core'])
addpath([workpath, filesep, 'toolbox'])

%% Configuration
info = [];

% Paths
info.workpath = workpath;
info.datapath = workpath;

% Data frequency: 'monthly', 'weekly', or 'daily'
info.freq = 'monthly';

% Include inflation swaps as 4th variable (true/false)
info.include_inflation = true;

% Evaluation window (eval_T1 = '9999-12-31' always uses last available date)
info.eval_T1 = '9999-12-31';

% VAR
info.p   = 0; % number of lags
info.nex = 1; % include constant

% Prior
% info.T0 is set automatically to 10% of sample inside estimation_RANDOMCORR_filter
info.kB  = 0.01; % prior scaling for VAR coefficients

% MCMC
info.ndraws  = 10000;
info.nburn   = 1000;
info.nthin   = 3;
info.nreport = 100;

%% Estimation
% Only delete prior FILTER outputs so the smoother run's mat is preserved
delete(fullfile(workpath, 'RANDOMCORR_FILTER_*.mat'));

disp('Starting filter-only estimation ...');
tStart = tic;
estimation_RANDOMCORR_filter(info);
tElapsed = toc(tStart);
disp(['Done. Total time: ', num2str(tElapsed), ' sec']);

%% Plot correlations and export PDF
pdfname = ['correlations_filter_', info.freq, '.pdf'];
if exist(fullfile(workpath, pdfname), 'file')
    delete(fullfile(workpath, pdfname));
end

mat_files = dir(fullfile(workpath, 'RANDOMCORR_FILTER_*.mat'));
if isempty(mat_files)
    warning('No .mat result file found — skipping plots.');
else
    [~, latest_idx] = max([mat_files.datenum]);
    matfile = fullfile(workpath, mat_files(latest_idx).name);
    load(matfile, 'r');

    X       = r.P;
    X_mean  = mean(X, 4);
    X_std   = std(X, 0, 4);
    X_plus  = X_mean + X_std;
    X_minus = X_mean - X_std;

    dates   = build_correlation_dates(fullfile(workpath, 'data.csv'), info.p, info.freq, info.include_inflation, matfile, size(X, 3));
    pdfpath = fullfile(workpath, pdfname);

    pairs = {[1,2], 'Equity - Bonds (filter-only)'; ...
             [1,3], 'Equity - Commodities (filter-only)'; ...
             [2,3], 'Bonds - Commodities (filter-only)'};
    if isfield(info, 'include_inflation') && info.include_inflation
        pairs = [pairs; {[1,4], 'Equity - Inflation (filter-only)'; ...
                         [2,4], 'Bonds - Inflation (filter-only)'; ...
                         [3,4], 'Commodities - Inflation (filter-only)'}];
    end

    n_pairs = size(pairs, 1);
    if n_pairs <= 3
        nrows = 1; ncols = 3;
    else
        nrows = 2; ncols = 3;
    end

    fig = figure('Color', 'w', 'Position', [100 100 1400 800]);
    for kk = 1:n_pairs
        ij = pairs{kk, 1};
        ti = pairs{kk, 2};
        subplot(nrows, ncols, kk);
        plot_correlation_axes(dates, ...
            squeeze(X_mean(ij(1), ij(2), :)), ...
            squeeze(X_plus(ij(1), ij(2), :)), ...
            squeeze(X_minus(ij(1), ij(2), :)), ti);
    end
    exportgraphics(fig, pdfpath);

    disp(['Plots saved to: ', pdfpath]);
    %delete(fullfile(workpath, 'RANDOMCORR_FILTER_*.mat'));
end

%% ---- Local functions ----

function corr_dates = build_correlation_dates(datafile, p, freq, include_inflation, matfile, n_obs)
tmp         = readtable(datafile);
dates_raw   = datetime(tmp.obs_date);
if include_inflation
    data_raw = [tmp.mars_equities_portfolio, tmp.mars_bonds_portfolio, tmp.mars_commodities_portfolio, tmp.mars_inflation_portfolio];
else
    data_raw = [tmp.mars_equities_portfolio, tmp.mars_bonds_portfolio, tmp.mars_commodities_portfolio];
end
valid_idx   = all(~isnan(data_raw), 2);
dates_daily = dates_raw(valid_idx);

switch lower(freq)
    case 'daily'
        Xcal = dates_daily(p+1:end);
    case 'weekly'
        grp = year(dates_daily)*100 + week(dates_daily, 'weekofyear');
        [unique_grp, ~, ic] = unique(grp);
        period_end_dates = NaT(numel(unique_grp), 1);
        for mm = 1:numel(unique_grp)
            pd = dates_daily(ic == mm);
            period_end_dates(mm) = pd(end);
        end
        Xcal = period_end_dates(p+1:end);
    otherwise % monthly
        grp = year(dates_daily)*100 + month(dates_daily);
        [unique_grp, ~, ic] = unique(grp);
        period_end_dates = NaT(numel(unique_grp), 1);
        for mm = 1:numel(unique_grp)
            pd = dates_daily(ic == mm);
            period_end_dates(mm) = pd(end);
        end
        Xcal = period_end_dates(p+1:end);
end

tokens = regexp(matfile, '_(\d{4}-\d{2}-\d{2})\.mat$', 'tokens', 'once');
if isempty(tokens)
    error('Could not parse evaluation date from filename: %s', matfile);
end
eval_date = datetime(tokens{1}, 'InputFormat', 'yyyy-MM-dd');
eval_idx  = find(Xcal == eval_date, 1);
if isempty(eval_idx)
    [~, eval_idx] = min(abs(days(Xcal - eval_date)));
    warning('Evaluation date not found exactly; using nearest: %s', datestr(Xcal(eval_idx), 'yyyy-mm-dd'));
end

start_idx = eval_idx - n_obs + 1;
if start_idx < 1
    error('Not enough dates to plot %d observations ending at %s.', n_obs, datestr(Xcal(eval_idx), 'yyyy-mm-dd'));
end
corr_dates = Xcal(start_idx:eval_idx);
end

function plot_correlation_axes(dates, mean_vals, upper_vals, lower_vals, plot_title)
dates      = dates(:);
mean_vals  = mean_vals(:);
upper_vals = upper_vals(:);
lower_vals = lower_vals(:);

plot(dates, mean_vals,  'k',   'LineWidth', 2);
hold on
plot(dates, upper_vals, 'r--', 'LineWidth', 1.25);
plot(dates, lower_vals, 'r--', 'LineWidth', 1.25);
yline(0, ':', 'Color', [0.4 0.4 0.4]);
title(plot_title);
xlabel('Date');
ylabel('Correlation');
legend('Mean', '+1 std', '-1 std', 'Location', 'best');
grid on; box on;
xlim([dates(1), dates(end)]);
ylim([-1, 1]);
xtickformat('yyyy');
end
