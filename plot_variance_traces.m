% Plot trace + ACF for every sig2h and sig2r entry from the most recent
% RANDOMCORR_monthly_*.mat (assumes it was produced by convergence_test.m
% with nthin = 1 and a leading info.nburn portion).

clear variables; close all;
workpath = pwd;

% Burn-in used inside convergence_test.m
nburn = 2500;

mat_files = dir(fullfile(workpath, 'RANDOMCORR_monthly_*.mat'));
[~, latest_idx] = max([mat_files.datenum]);
matfile = fullfile(workpath, mat_files(latest_idx).name);
S = load(matfile, 'r'); r = S.r;

mat_sig2h = r.sig2h(nburn+1:end, :);   % N x m
mat_sig2r = r.sig2r(nburn+1:end, :);   % N x n_r

m   = size(mat_sig2h, 2);
n_r = size(mat_sig2r, 2);
N   = size(mat_sig2h, 1);
fprintf('Loaded %s\n', matfile);
fprintf('  m = %d, n_r = %d, post-burn-in N = %d\n', m, n_r, N);

%% ACF helper
max_lag = min(200, floor(N/4));
acf_of = @(x) local_acf(x, max_lag);

%% Variable labels — assumes [equities, bonds, commodities, inflation]
var_names = {'eq', 'bd', 'cm', 'infl'};
sig2h_labels = arrayfun(@(i) sprintf('sig2h_{%s}', var_names{i}), 1:m, 'UniformOutput', false);

% sig2r ordering follows tril(...,-1) column-major:
% for m=4 -> pairs in order (2,1),(3,1),(4,1),(3,2),(4,2),(4,3)
pair_labels = cell(n_r, 1);
k = 0;
for j = 1:m-1
    for i = j+1:m
        k = k + 1;
        pair_labels{k} = sprintf('sig2r_{%s-%s}', var_names{i}, var_names{j});
    end
end

%% ===== Figure 1: sig2h (m rows, trace + ACF columns) =====
fig1 = figure('Color','w','Position',[100 100 1300 200*m+80]);
for i = 1:m
    series = mat_sig2h(:, i);
    subplot(m, 2, 2*i-1);
    plot(series, 'LineWidth', 0.8); grid on;
    title(['Trace: ', sig2h_labels{i}]); xlabel('Draw'); xlim([1 N]);

    subplot(m, 2, 2*i);
    stem(0:max_lag, acf_of(series), 'Marker','none'); grid on;
    yline(0.1, 'r--'); yline(-0.1, 'r--');
    title(['ACF: ', sig2h_labels{i}]); xlabel('Lag'); ylim([-0.3 1.05]);
end
sgtitle('Variance-of-log-volatility innovations: sig2h');
exportgraphics(fig1, fullfile(workpath, 'traces_sig2h.pdf'));
fprintf('Saved traces_sig2h.pdf\n');

%% ===== Figure 2: sig2r (n_r rows) =====
fig2 = figure('Color','w','Position',[100 100 1300 200*n_r+80]);
for i = 1:n_r
    series = mat_sig2r(:, i);
    subplot(n_r, 2, 2*i-1);
    plot(series, 'LineWidth', 0.8); grid on;
    title(['Trace: ', pair_labels{i}]); xlabel('Draw'); xlim([1 N]);

    subplot(n_r, 2, 2*i);
    stem(0:max_lag, acf_of(series), 'Marker','none'); grid on;
    yline(0.1, 'r--'); yline(-0.1, 'r--');
    title(['ACF: ', pair_labels{i}]); xlabel('Lag'); ylim([-0.3 1.05]);
end
sgtitle('Variance-of-correlation-walk innovations: sig2r');
exportgraphics(fig2, fullfile(workpath, 'traces_sig2r.pdf'));
fprintf('Saved traces_sig2r.pdf\n');

%% ===== Combined: all on one figure =====
total = m + n_r;
fig3 = figure('Color','w','Position',[100 100 1300 160*total+80]);
row = 0;
for i = 1:m
    row = row + 1;
    series = mat_sig2h(:, i);
    subplot(total, 2, 2*row-1);
    plot(series, 'LineWidth', 0.8); grid on;
    title(['Trace: ', sig2h_labels{i}]); xlabel('Draw'); xlim([1 N]);

    subplot(total, 2, 2*row);
    stem(0:max_lag, acf_of(series), 'Marker','none'); grid on;
    yline(0.1, 'r--'); yline(-0.1, 'r--');
    title(['ACF: ', sig2h_labels{i}]); xlabel('Lag'); ylim([-0.3 1.05]);
end
for i = 1:n_r
    row = row + 1;
    series = mat_sig2r(:, i);
    subplot(total, 2, 2*row-1);
    plot(series, 'LineWidth', 0.8); grid on;
    title(['Trace: ', pair_labels{i}]); xlabel('Draw'); xlim([1 N]);

    subplot(total, 2, 2*row);
    stem(0:max_lag, acf_of(series), 'Marker','none'); grid on;
    yline(0.1, 'r--'); yline(-0.1, 'r--');
    title(['ACF: ', pair_labels{i}]); xlabel('Lag'); ylim([-0.3 1.05]);
end
sgtitle('All variance hyperparameters');
exportgraphics(fig3, fullfile(workpath, 'traces_sig2_all.pdf'));
fprintf('Saved traces_sig2_all.pdf\n');

%% ---- local function ----
function a = local_acf(x, max_lag)
    x = x(:) - mean(x);
    v = mean(x.^2);
    N = numel(x);
    a = zeros(max_lag+1, 1);
    if v < 1e-14
        a(1) = 1; return
    end
    for L = 0:max_lag
        a(L+1) = mean(x(1:N-L) .* x(L+1:N)) / v;
    end
end
