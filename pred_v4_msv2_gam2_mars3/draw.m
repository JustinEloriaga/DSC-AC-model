clear all
close all
clc

load msv2_mars3_long_1_pred_vMARS3_t_1990-01-31.mat

X=r.P;
% Compute mean and standard deviation across draws (4th dimension)
X_mean  = mean(X, 4);
X_std   = std(X, 0, 4);

% Construct +1 std and -1 std matrices
X_plus  = X_mean + X_std;
X_minus = X_mean - X_std;

% Horizon
h = 1:25;

% -------- Entry (1,2) --------
m = squeeze(X_mean(1,2,:));
p = squeeze(X_plus(1,2,:));
n = squeeze(X_minus(1,2,:));

figure;
plot(h, m, 'k', 'LineWidth', 2); hold on;
plot(h, p, 'r--', 'LineWidth', 1.5);
plot(h, n, 'r--', 'LineWidth', 1.5);
title('Entry (1,2)');
xlabel('Horizon');
ylabel('Value');
legend('Mean', '+1 std', '-1 std', 'Location', 'best');
grid on;

% -------- Entry (1,3) --------
m = squeeze(X_mean(1,3,:));
p = squeeze(X_plus(1,3,:));
n = squeeze(X_minus(1,3,:));

figure;
plot(h, m, 'k', 'LineWidth', 2); hold on;
plot(h, p, 'r--', 'LineWidth', 1.5);
plot(h, n, 'r--', 'LineWidth', 1.5);
title('Entry (1,3)');
xlabel('Horizon');
ylabel('Value');
legend('Mean', '+1 std', '-1 std', 'Location', 'best');
grid on;

% -------- Entry (2,3) --------
m = squeeze(X_mean(2,3,:));
p = squeeze(X_plus(2,3,:));
n = squeeze(X_minus(2,3,:));

figure;
plot(h, m, 'k', 'LineWidth', 2); hold on;
plot(h, p, 'r--', 'LineWidth', 1.5);
plot(h, n, 'r--', 'LineWidth', 1.5);
title('Entry (2,3)');
xlabel('Horizon');
ylabel('Value');
legend('Mean', '+1 std', '-1 std', 'Location', 'best');
grid on;
