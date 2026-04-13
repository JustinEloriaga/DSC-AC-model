% Elliptical sampling for SV
clc; clear all; close all;

% 08_all: putting everything together (sampling everything)


%% DGP

T= 200;

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
A     = [30, -5; -5, 30];
d     = 0.85;
k     = 100;

Q00      = eye(m); %we assume that this is known
Q00(1,2) = 0.9;
Q00(2,1) = 0.9;

[Pt_old,Qt_old,Qinvt_old] = gen_data_sdc(T, A, d, k, Q00);
Pt = Pt_old;
% plot(squeeze(Pt(1,2,:)));

% Data
yt = zeros(T,m);
for t=1:T
    yt(t,:) = diag(exp(ht_old(t,:)/2))*(chol(Pt_old(:,:,t),'lower')*randn(m,1));
end

%% Prior and tuning

% initial Q00
Q00      = eye(m); %we assume that this is known
Q00(1,2) = 0.9;
Q00(2,1) = 0.9;

% prior for Ainv (following Asai and McAleer)
prior_Ainv.gam  = m + 30;
prior_Ainv.Cinv = ([1, 0.3; 0.3, 1]) * (prior_Ainv.gam);

% prior d -> uniform over [L,U]
prior_d.L = 0.5;
prior_d.U = 0.95;
prior_d.mhtune_a = 0.3; %target acceptance rate
prior_d.mhtune_c = 0.55; %cooling rate (between 0.5 and 1), c = inf means no adaptation
prior_d.mhtune_b = 10000; %bound
prior_d.mhtune_logvar = log(0.1^2); %log proposal variance (initial proposal variance)

% prior k -> truncated exponential
prior_k.lam0 = 0.01;
prior_k.mhtune_a = 0.3; %target acceptance rate
prior_k.mhtune_c = 0.55; %cooling rate (between 0.5 and 1), c = inf means no adaptation
prior_k.mhtune_b = 10000; %bound
prior_k.mhtune_logvar = log(0.1^2); %log proposal variance (initial proposal variance)

% prior W (PRIOR as in Primiceri)
% W ~ IW(kw^2*4*I(m), 4)
% T0H: degrees of freedom for the IW prior distribution of var(errH)
% kH: constant scaling the prior of var(errH)
T0H=4;
kH=.01;

prior_sig2h.T0 = ones(m,1)*2;
prior_sig2h.V0 = ones(m,1)*2;

%% Setting up slice sampler (common part)
% h(t) are independent RW

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

%% Initilization
et_old   = yt;
Ainv_old = wishrnd(eye(m)/prior_Ainv.Cinv, prior_Ainv.gam);
A_old    = eye(m)/Ainv_old;
d_old    = 0.8;
k_old    = 20;
[Pt_old, Qt_old, Qinvt_old] = gen_data_sdc(T, A_old, d_old, k_old, Q00); %initialize
sig2h_old = sig2h;


%% Actual sampler
nsim = 1000;

mat_Pt = zeros(nsim, T);
mat_Ainv = zeros(m,m,nsim);
mat_d = zeros(nsim,1);
mat_d_rej = zeros(nsim,1);
mat_k = zeros(nsim,1);
mat_k_rej = zeros(nsim,1);
mat_Q_rej = zeros(nsim,T);
mat_ht = zeros(T, m, nsim);
mat_lik = zeros(T,1);
% mat_W = zeros(m,m,nsim);
mat_sig2h = zeros(nsim,m);

for sind = 1:nsim
    
    % ---
    % Dynamic correlation
    % standardized data
    zt_old = et_old ./ exp(ht_old/2);
    
    % draw Qt,Pt given zt and Ainv (MH)
    [Qt_old, Qinvt_old, Q_rej] = gen_Qt(Qt_old, Qinvt_old, d_old, k_old, Ainv_old, Q00, zt_old);
    for t=1:T
        Pt_old(:,:,t) = QtoP(Qt_old(:,:,t));
    end
    
    % draw Ainv given data
    Ainv_old = gen_Ainv(d_old, k_old, Q00, Qt_old, Qinvt_old, prior_Ainv);
    
    % drawing d given data
    [d_old, rej_d, prior_d] = gen_d(Qt_old, Qinvt_old, d_old, k_old, Ainv_old, Q00, prior_d, sind);
    
    % drawing k given data
    [k_old, rej_k, prior_k] = gen_k(Qt_old, Qinvt_old, d_old, k_old, Ainv_old, Q00, prior_k, sind);
    
    % ---
    % Dynamic variance
    % draw ht (log volaility)
    % prep-slice sampler (update parts that change over time)
    lik_old = loglike_yt_given_ht_Pt(et_old,ht_old,Pt_old,1,ht_old(:,1)); %have to evaluate likelihood again with new Pt_old
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
    %V1H = errH'*errH + (kH^2)*eye(m)*T0H;
    %W_old = iwishrnd(V1H, T-1+T0H);
    for vvind=1:m
        EE = errH(:,vvind);
        T1 = size(errH,1)/2 + prior_sig2h.T0(vvind);
        V1 = EE'*EE/2 + prior_sig2h.V0(vvind);
        sig2h_old(vvind) = gamrnd(T1, V1^(-1))^(-1);
    end
    

    % report
    if rem(sind,100) == 0
        disp(['sind = ', num2str(sind), ' / ', num2str(nsim)]);
    end
    
    % store correlation and A
    mat_Pt(sind,:) = (squeeze(Pt_old(1,2,:)));
    mat_Ainv(:,:,sind) = Ainv_old;
    mat_Q_rej(sind,:) = Q_rej;
    
    mat_d_rej(sind,:) = rej_d;
    mat_d(sind,1) = d_old;
    
    mat_k_rej(sind,:) = rej_k;
    mat_k(sind,1) = k_old;
    
    mat_ht(:,:,sind)  = ht_old;
    mat_lik(sind,1) = lik_old;
    
    mat_sig2h(sind,:) = sig2h_old;
end


%% Ploting

figure(1)
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

figure(2)
m_Pt = quantile(mat_Pt(100:end,:), 0.5, 1)';
q1_Pt = quantile(mat_Pt(100:end,:), 0.1, 1)';
q2_Pt = quantile(mat_Pt(100:end,:), 0.9, 1)';

plot(squeeze(Pt(1,2,:)));
hold on
plot(m_Pt);
plot(q1_Pt, 'g');
plot(q2_Pt, 'g');
hold off



