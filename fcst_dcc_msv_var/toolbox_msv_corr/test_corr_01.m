% Elliptical sampling for SV and correlation
clc; clear all; close all;

% _corr_01: Peter Hasen's correlation decomposition


%% DGP

T= 100;

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

% [Pt_old,Qt_old,Qinvt_old] = gen_data_sdc(T, A, d, k, Q00);
% Pt_true = Pt_old;
Pt_true = zeros(m,m,T);
for t=1:T
    
    if t<T/2
        temp_P = eye(m);
        temp_P(1,2) = -0.3;
        temp_P(2,1) = -0.3;
    else
        temp_P = eye(m);
        temp_P(1,2) = 0.9;
        temp_P(2,1) = 0.9;
    end
    
    Pt_true(:,:,t) = temp_P;
end
Pt_old = Pt_true;

% Data
yt = zeros(T,m);
for t=1:T
    yt(t,:) = diag(exp(ht_old(t,:)/2))*(chol(Pt_old(:,:,t),'lower')*randn(m,1));
end

%% Mean and covariance vector

% -----------------------------------------------
% log-volatility
% setting up slice sampler (common part)
slice.nobs = T;

% h1
vvind = 1;
var_z = Vh0_scale(vvind) + (1:1:T)';
cov_z = tril(repmat(var_z', T, 1), 0) + tril(repmat(var_z', T, 1), -1)';
chol_cov_z_h1 = chol(cov_z,'lower');
mean_z_h1 = mh0(vvind)*ones(T,1);

% h2
vvind = 2;
var_z = Vh0_scale(vvind) + (1:1:T)';
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


% -----------------------------------------------
% correlation elements
n_r = 1; % # of correlation elements
Vr0_scale = 10*ones(n_r,1);
mr0 = 0*ones(n_r,1);

% r1
vvind = 1;
var_z = Vr0_scale(vvind) + (1:1:T)';
cov_z = tril(repmat(var_z', T, 1), 0) + tril(repmat(var_z', T, 1), -1)';
chol_cov_z_r1 = chol(cov_z,'lower');
mean_z_r1 = mr0(vvind)*ones(T,1);

slice_chol_cov_z_r = cat(3, chol_cov_z_r1);
slice_mean_z_r = [mean_z_r1];

% initialize r
rt_old = 0.3*ones(T,1);
Pt_old = zeros(m,m,T);
for t=1:T
    Pt_old(:,:,t) = veclAtoC(rt_old(t,:));
end

% initialize sig2r
sig2r_old = 0.1*ones(n_r,1);

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


mat_rt = zeros(T, n_r, nsim);
mat_n_r = zeros(nsim,1);
mat_sig2r = zeros(nsim,n_r);
mat_Pt = zeros(m,m,T,nsim);

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
    

    % draw rt (log volaility)
    % prep-slice sampler (update parts that change over time)
    for vvind = 1:n_r
        slice.scale_z = sqrt(sig2r_old(vvind));
        slice.mean    = slice_mean_z_r(:,vvind);
        slice.chol_cov_z = slice_chol_cov_z_r(:,:,vvind);
        slice.fcn_lik = @(rt_vvind) loglike_yt_given_rt_ht(et_old,ht_old,rt_old,vvind,rt_vvind);
        
        [rt_vvind_old, lik_old, n_try_r] = slice_sampling_v02(slice, rt_old(:,vvind), lik_old);
        rt_old(:,vvind) = rt_vvind_old;
    end
       
    % from rt  to Pt
    for t=1:T
        Pt_old(:,:,t) = veclAtoC(rt_old(t,:));
    end

    % store
    mat_ht(:,:,sind)  = ht_old;
    mat_lik(sind,1) = lik_old;
    mat_n(sind,1)   = n_try;
    mat_sig2h(sind,:) = sig2h_old;
    
    mat_rt(:,:,sind) = rt_old;
    mat_sig2r(sind,:) = sig2r_old;
    mat_n_r(sind,1)   = n_try_r;
    mat_lik(sind,1) = lik_old;
    mat_Pt(:,:,:,sind) = Pt_old;
    
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

%% Plotting - correlation

figure(101)

m_Pt = mean(mat_Pt,4);
q1_Pt = quantile(mat_Pt,0.1,4);
q2_Pt = quantile(mat_Pt,0.9,4);
plot(squeeze(m_Pt(1,2,:)))
hold on
plot(squeeze(q1_Pt(1,2,:)), 'g')
plot(squeeze(q2_Pt(1,2,:)), 'g')

plot(squeeze(Pt_true(1,2,:)), 'r');

hold off





