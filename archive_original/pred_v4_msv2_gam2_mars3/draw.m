clear
close all
clc

if exist('correlations.pdf','file'); delete('correlations.pdf'); end

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
repo_root = fileparts(script_dir);

addpath(repo_root)

mat_files = dir(fullfile(script_dir, 'msv2_mars3_long_1_pred_vMARS3_t_*.mat'));
if isempty(mat_files)
    error('No correlation result files found in %s', script_dir);
end

[~, latest_idx] = max([mat_files.datenum]);
matfile = fullfile(script_dir, mat_files(latest_idx).name);
load(matfile, 'r');

X = r.P;

% Compute mean and standard deviation across saved draws (4th dimension).
X_mean = mean(X, 4);
X_std = std(X, 0, 4);
X_plus = X_mean + X_std;
X_minus = X_mean - X_std;

info = get_default_info_v3_msv2_gam2_mars3();
dates = build_correlation_dates(fullfile(repo_root, 'data_JRR.csv'), info.p, matfile, size(X, 3));

plot_correlation_series(dates, squeeze(X_mean(1,2,:)), squeeze(X_plus(1,2,:)), squeeze(X_minus(1,2,:)), 'Equity - Bonds');
exportgraphics(gcf, 'correlations.pdf', 'Append', true);
plot_correlation_series(dates, squeeze(X_mean(1,3,:)), squeeze(X_plus(1,3,:)), squeeze(X_minus(1,3,:)), 'Equity - Commodities');
exportgraphics(gcf, 'correlations.pdf', 'Append', true);
plot_correlation_series(dates, squeeze(X_mean(2,3,:)), squeeze(X_plus(2,3,:)), squeeze(X_minus(2,3,:)), 'Bonds - Commodities');
exportgraphics(gcf, 'correlations.pdf', 'Append', true);

function corr_dates = build_correlation_dates(datafile, p, matfile, n_obs)
tmp = readtable(datafile);

dates_raw = datetime(tmp.obs_date);
data_raw = [tmp.mars_equities_portfolio, ...
            tmp.mars_bonds_portfolio, ...
            tmp.mars_commodities_portfolio];

valid_idx = all(~isnan(data_raw), 2);
dates_daily = dates_raw(valid_idx);

ym = year(dates_daily) * 100 + month(dates_daily);
[unique_ym, ~, ic] = unique(ym);

month_end_dates = NaT(numel(unique_ym), 1);
for mm = 1:numel(unique_ym)
    month_dates = dates_daily(ic == mm);
    month_end_dates(mm) = month_dates(end);
end

Xcal = month_end_dates(p+1:end);

tokens = regexp(matfile, '_t_(\d{4}-\d{2}-\d{2})\.mat$', 'tokens', 'once');
if isempty(tokens)
    error('Could not parse evaluation date from %s', matfile);
end

eval_date = datetime(tokens{1}, 'InputFormat', 'yyyy-MM-dd');
eval_idx = find(Xcal == eval_date, 1);
if isempty(eval_idx)
    [~, eval_idx] = min(abs(days(Xcal - eval_date)));
    warning('Evaluation date %s not found exactly. Using nearest plotted date %s.', ...
        datestr(eval_date, 'yyyy-mm-dd'), datestr(Xcal(eval_idx), 'yyyy-mm-dd'));
end

start_idx = eval_idx - n_obs + 1;
if start_idx < 1
    error('Not enough monthly dates to plot %d correlation observations ending at %s.', ...
        n_obs, datestr(Xcal(eval_idx), 'yyyy-mm-dd'));
end

corr_dates = Xcal(start_idx:eval_idx);
end

function plot_correlation_series(dates, mean_vals, upper_vals, lower_vals, plot_title)
dates = dates(:);
mean_vals = mean_vals(:);
upper_vals = upper_vals(:);
lower_vals = lower_vals(:);

figure('Color', 'w');
plot(dates, mean_vals, 'k', 'LineWidth', 2);
hold on
plot(dates, upper_vals, 'r--', 'LineWidth', 1.25);
plot(dates, lower_vals, 'r--', 'LineWidth', 1.25);
yline(0, ':', 'Color', [0.4 0.4 0.4]);

title(plot_title);
xlabel('Date');
ylabel('Correlation');
legend('Mean', '+1 std', '-1 std', 'Location', 'best');
grid on
box on
xlim([dates(1), dates(end)]);
ylim([-1, 1]);
xtickformat('yyyy');
end
