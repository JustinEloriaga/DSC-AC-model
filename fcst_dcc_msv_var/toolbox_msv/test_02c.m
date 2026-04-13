% Wishart based volatility model from Asai and McAleer (Joe, 2009)

% test_02b: Qinv generation testing for DGP that is different from the
% model (structural break)


clc; clear all; close all;

%% Unknowns
% A, d, k

A     = inv([30, -3; -3, 30]);
% A = inv([1, 0.3; 0.3, 1]);
d     = 0.80;
k     = 15;


%% Data generation
T = 300;
m = size(A,1);
Q00      = eye(m); %we assume that this is known
Q00(1,2) = 0.9;
Q00(2,1) = 0.9;

% [Pt, Qt, Qinvt] = gen_data_sdc(T, A, d, k, Q00);
Pt = ones(m,m,T);
for t=1:T
    
    if t<100
        temp_r = 0.9;
    else
        temp_r = 0.5;
    end
    Pt(1,2,t) = temp_r;
    Pt(2,1,t) = temp_r;
end
zt = gen_data_zt(Pt);

% mean(squeeze(Pt(1,2,:)))
% corr(zt)
% plot(squeeze(Pt(1,2,:)))

%% Prior
% prior for Ainv (following Asai and McAleer)
prior_Ainv.gam  = m + 20;
prior_Ainv.Cinv = inv([30, -3; -3, 30]) * (prior_Ainv.gam);

% Ainv_new = gen_Ainv(d, k, Q00, Qt, Qinvt, prior_Ainv);

%%
Ainv_old = inv(A);
d_old = d;
k_old = k;

% Qt_old = Qt;
% Qinvt_old = Qinvt;

[~, Qt_old, Qinvt_old] = gen_data_sdc(T, A, d, k, Q00); %initialize

%% MCMC testing...

nsim = 1000;
mat_Q_rej = zeros(nsim,T);
mat_Pt = zeros(nsim,T);

for sind = 1:nsim
    
    
    % MH-generation of Qt
    [Qt_old, Qinvt_old, Q_rej] = gen_Qt(Qt_old, Qinvt_old, d_old, k_old, Ainv_old, Q00, zt);
    
    % draw Ainv given data
    Ainv_old = gen_Ainv(d_old, k_old, Q00, Qt_old, Qinvt_old, prior_Ainv);
    %A_old = eye(m)/Ainv_old;
    
    % store correlation draw
    for t=1:T
        temp_P = QtoP(Qt_old(:,:,t));
        mat_Pt(sind,t) = temp_P(1,2);
    end
    mat_Q_rej(sind,:) = Q_rej;
    
end



disp(mean(mat_Q_rej,1));


%% Figure
plot(squeeze(Pt(1,2,:)))
q1 = quantile(mat_Pt, 0.1, 1)';
q2 = quantile(mat_Pt, 0.9, 1)';
q3 = quantile(mat_Pt, 0.5, 1)';
hold on
plot(q1, 'r');
plot(q2, 'r');
plot(q3, 'g');
hold off
