% Elliptical sampling for SV
clc; clear all; close all;

% elliptical sampling for log volatility, univariate

% geweke test


%% DGP
T= 5;

% Random walk process
sig2h = 0.5;
mh0 = 0;
Vh0_scale = 1^2;
Vh0 = sig2h*Vh0_scale; %for convenience

ht = zeros(T,1);
h0 = mh0 + sqrt(Vh0)*randn;
for t=1:T
    h0 = h0 + sqrt(sig2h)*randn;
    ht(t,1) = h0;
end

yt = zeros(T,1);
for t=1:T
    yt(t,1) = exp(ht(t)/2)*randn;
end

%% Mean and covariance vector

% setting up slice sampler (common part)
slice.nobs = T;
var_z = Vh0_scale + (1:1:T)';
cov_z = tril(repmat(var_z', T, 1), 0) + triu(repmat(var_z', T, 1), 1);
slice.chol_cov_z = chol(cov_z,'lower');
slice.mean = mh0*ones(T,1);

%% Volatility sampling exercise
sig2h_old = sig2h;
et_old = yt;


% initialize
ht_old = ht;
slice.fcn_lik = @(ht_old) sum(log(normpdf(et_old, 0, exp(ht_old/2))));
lik_old = slice.fcn_lik(ht_old);

% % actual mcmc
% nsim = 1000;
% mat_ht = zeros(T, nsim);
% mat_lik = zeros(T,1);
% mat_n = zeros(T,1);
% for sind = 1:nsim
%     
%     % draw ht (log volaility)
%     % prep-slice sampler (update parts that change over time)
%     slice.scale_z = sqrt(sig2h_old);
%     slice.fcn_lik = @(ht_old) sum(log(normpdf(et_old, 0, exp(ht_old/2))));
%     [ht_old, lik_old, n_try] = slice_sampling_v02(slice, ht_old, lik_old);
%     
%     % store
%     mat_ht(:,sind) = ht_old;
%     mat_lik(sind,1) = lik_old;
%     mat_n(sind,1) = n_try;
%     
%     % report
%     if rem(sind,100) == 0
%         disp(['sind = ', num2str(sind), ' / ', num2str(nsim)]);
%     end
% end

%% Ploting
% m_ht = quantile(mat_ht, 0.5, 2);
% q1_ht = quantile(mat_ht, 0.1, 2);
% q2_ht = quantile(mat_ht, 0.9, 2);
% 
% plot(ht);
% hold on
% plot(m_ht);
% plot(q1_ht, 'g');
% plot(q2_ht, 'g');
% hold off

%% Geweke test


% method 1
disp('method 1 ....');
% actual mcmc
nsim = 50000;
m1_ht = zeros(T, nsim);
for sind = 1:nsim
    
    % generate h(t) from the prior
    ht = zeros(T,1);
    h0 = mh0 + sqrt(Vh0)*randn;
    for t=1:T
        h0 = h0 + sqrt(sig2h)*randn;
        ht(t,1) = h0;
    end
    
    yt = zeros(T,1);
    for t=1:T
        yt(t,1) = exp(ht(t)/2)*randn;
    end
    % store
    m1_ht(:,sind) = ht;
    
    % report
    if rem(sind,100) == 0
        disp(['sind = ', num2str(sind), ' / ', num2str(nsim)]);
    end
end


% method 2
disp('method 2 ....');
% actual mcmc
m2_ht = zeros(T, nsim);
for sind = 1:nsim
    
    
    % Data conditional on ht
    yt = zeros(T,1);
    for t=1:T
        yt(t,1) = exp(ht_old(t)/2)*randn;
    end
    
    
    % h(t) conditional on data
    % draw ht (log volaility)
    % prep-slice sampler (update parts that change over time)
    et_old = yt;
    slice.scale_z = sqrt(sig2h_old);
    slice.fcn_lik = @(ht_old) sum(log(normpdf(et_old, 0, exp(ht_old/2))));
    lik_old = slice.fcn_lik(ht_old);
    [ht_old, lik_old, n_try] = slice_sampling_v02(slice, ht_old, lik_old);
    
    % store
    m2_ht(:,sind) = ht_old;
    
    % report
    if rem(sind,1000) == 0
        disp(['sind = ', num2str(sind), ' / ', num2str(nsim)]);
    end
end


%%
figure(1)
hist([m1_ht(1,:)', m2_ht(1,:)'],100)

figure(2)
hist([m1_ht(2,:)', m2_ht(2,:)'],100)

figure(3)
hist([m1_ht(3,:)', m2_ht(3,:)'],100)

figure(4)
hist([m1_ht(4,:)', m2_ht(4,:)'],100)

figure(5)
hist([m1_ht(5,:)', m2_ht(5,:)'],100)





