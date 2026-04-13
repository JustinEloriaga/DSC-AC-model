% Elliptical sampling for SV
clc; clear all; close all;

% 07: elliptical sampling for log volatility, bivariate


%% DGP

T= 1000;

% Random walk process
m = 2;
sig2h = 0.1*ones(m,1);
mh0 = ones(m,1);
Vh0_scale = 2^2*ones(m,1);
Vh0 = sig2h.*Vh0_scale; %for convenience
h0 = mh0 + diag(sqrt(Vh0))*randn(m,1);

ht_old = zeros(T,m);
for t=1:T
    h0 = h0 + diag(sqrt(sig2h))*randn(m,1);
    ht_old(t,:) = h0;
end
ht = ht_old;

% Correlation process
A     = [30, -3; -3, 30];
d     = 0.8;
k     = 50;

Q00      = eye(m); %we assume that this is known
Q00(1,2) = 0.9;
Q00(2,1) = 0.9;

[Pt_old,Qt_old,Qinvt_old] = gen_data_sdc(T, A, d, k, Q00);

% Data
yt = zeros(T,m);
for t=1:T
    yt(t,:) = diag(exp(ht_old(t,:)/2))*(chol(Pt_old(:,:,t),'lower')*randn(m,1));
end

%% Mean and covariance vector

% setting up slice sampler (common part)

slice.nobs = T;

% h1
vvind = 1;
var_z = Vh0_scale(vvind) + (1:1:T)';
% cov_z = tril(repmat(var_z', T, 1), 0) + triu(repmat(var_z', T, 1), 1); %old and wrong
cov_z = tril(repmat(var_z', T, 1), 0) + tril(repmat(var_z', T, 1), -1)';
chol_cov_z_h1 = chol(cov_z,'lower');
mean_z_h1 = mh0(vvind)*ones(T,1);

% h2
vvind = 2;
var_z = Vh0_scale(vvind) + (1:1:T)';
% cov_z = tril(repmat(var_z', T, 1), 0) + triu(repmat(var_z', T, 1), 1); %old and wrong
cov_z = tril(repmat(var_z', T, 1), 0) + tril(repmat(var_z', T, 1), -1)';
chol_cov_z_h2 = chol(cov_z,'lower');
mean_z_h2 = mh0(vvind)*ones(T,1);

slice_chol_cov_z = cat(3, chol_cov_z_h1, chol_cov_z_h2);
slice_mean_z = [mean_z_h1, mean_z_h2];

% -------------------
% For sig2h
% Prior for sig2h 
prior_sig2h.T0 = ones(m,1)*5; %first arg (as in Primiceri, but independent prior)
prior_sig2h.V0 = ones(m,1)*(0.01*5); %second arg ((2*0.01^2), as in Primiceri, but independent prior)
% sig2h initial
sig2h_old = prior_sig2h.V0 ./ (prior_sig2h.T0 - 1); % start from the prior mean

%% Volatility sampling exercise
% sig2h_old = sig2h;
et_old = yt;


% initialize
ht_old = ht;
lik_old = loglike_yt_given_ht_Pt(yt,ht_old,Pt_old,1,ht_old(:,1));

% actual mcmc
nsim = 1000;
mat_ht = zeros(T, m, nsim);
mat_lik = zeros(nsim,1);
mat_n = zeros(nsim,1);
mat_sig2h = zeros(nsim,m);

for sind = 1:nsim
    
    % draw ht (log volaility)
    % prep-slice sampler (update parts that change over time)
    for vvind = 1:m
        slice.scale_z = sqrt(sig2h_old(vvind));
        slice.mean    = slice_mean_z(:,vvind);
        slice.chol_cov_z = slice_chol_cov_z(:,:,vvind);
        slice.fcn_lik = @(ht_vvind) loglike_yt_given_ht_Pt(et_old,ht_old,Pt_old,vvind,ht_vvind);
        
        [ht_vvind_old, lik_old, n_try] = slice_sampling_v02(slice, ht_old(:,vvind), lik_old);
        ht_old(:,vvind) = ht_vvind_old;
    end
    
    
	% Drawing var(errH)
    errH = ht_old(2:end,:)-ht_old(1:end-1,:);
    for vvind = 1:m
        EE = errH(:,vvind);
        T1 = size(errH,1)/2 + prior_sig2h.T0(vvind);
        V1 = EE'*EE/2 + prior_sig2h.V0(vvind);
        sig2h_old(vvind) = gamrnd(T1, V1^(-1))^(-1);
    end
    
    % store
    mat_ht(:,:,sind)  = ht_old;
    mat_lik(sind,1) = lik_old;
    mat_n(sind,1)   = n_try;
    
    mat_sig2h(sind,:) = sig2h_old;
    % report
    if rem(sind,100) == 0
        disp(['sind = ', num2str(sind), ' / ', num2str(nsim)]);
    end
end


%% Ploting
m_ht = quantile(mat_ht, 0.5, 3);
q1_ht = quantile(mat_ht, 0.1, 3);
q2_ht = quantile(mat_ht, 0.9, 3);

subplot(2,1,1)
vind = 1;
plot(ht(:,vind));
hold on
plot(m_ht(:,vind));
plot(q1_ht(:,vind), 'g');
plot(q2_ht(:,vind), 'g');
hold off

subplot(2,1,2)
vind = 2;
plot(ht(:,vind));
hold on
plot(m_ht(:,vind));
plot(q1_ht(:,vind), 'g');
plot(q2_ht(:,vind), 'g');
hold off



