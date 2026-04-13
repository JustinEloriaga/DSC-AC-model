% Wishart based volatility model from Asai and McAleer (Joe, 2009)

% test 04: for testing generation of d

clc; clear all; close all;

%% Unknowns
% A, d, k

% A     = [30, -3; -3, 30];
A = inv([1, 0.3; 0.3, 1]);
d     = 0.9;
k     = 50;


%% Data generation
T = 20;
m = size(A,1);
Q00      = eye(m); %we assume that this is known
Q00(1,2) = 0.9;
Q00(2,1) = 0.9;

[Pt, Qt, Qinvt] = gen_data_sdc(T, A, d, k, Q00);
zt = gen_data_zt(Pt);

% mean(squeeze(Pt(1,2,:)))
% corr(zt)
% plot(squeeze(Pt(1,2,:)))

%% Prior and tuning

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

%% Initialization
Ainv_old = inv(A);
A_old = A;

d_old = d;
k_old = k;


Qt_old = Qt;
Qinvt_old = Qinvt;
Pt_old = Pt;
% [~, Qt_old, Qinvt_old] = gen_data_sdc(T, A, d, k, Q00); %initialize

%% test Sampling

% nsim = 1000;
% mat_d = zeros(nsim,1);
% mat_rej_d = zeros(nsim,1);
% 
% for sind = 1:nsim
%     
%     % generating d
%     [d_old, rej_d, prior_d] = gen_d(Qt_old, Qinvt_old, d_old, k_old, Ainv_old, Q00, prior_d, sind);
%     
%     % store
%     mat_d(sind,1) = d_old;
%     mat_rej_d(sind,1) = rej_d;
% end
% 
% xxxxxxxxxxxxxxxxxxxxxxxx


%% Geweke test, method 1
nsim = 500000;

disp('implementing method 1 ...');
m1_Pt = zeros(nsim, T);
m1_A = zeros(m,m,nsim);
m1_d = zeros(nsim,1);
for sind = 1:nsim
    
    % draw d from the prior
    d_old = prior_d.L + (prior_d.U-prior_d.L)*rand;
    
%     % draw Ainv from prior
%     Ainv_old = wishrnd(eye(m)/prior_Ainv.Cinv, prior_Ainv.gam);
%     A_old = eye(m)/Ainv_old;
    
    % draw data given new A draw
    [Pt, Qt, Qinvt] = gen_data_sdc(T, A_old, d_old, k_old, Q00);
    
    % store correlation and A
    m1_Pt(sind,:) = (squeeze(Pt(1,2,:)));
    m1_A(:,:,sind) = A_old;
    m1_d(sind,:) = d_old;
end


%% Geweke test, method 2
disp('implementing method 2 ...');
m2_Pt = zeros(nsim, T);
m2_A = zeros(m,m,nsim);
m2_d = zeros(nsim,1);
m2_d_rej = zeros(nsim,1);

% initialize with prior
% Ainv_old = wishrnd(eye(m)/prior_Ainv.Cinv, prior_Ainv.gam);
% A_old = eye(m)/Ainv_old;
% [Pt_old, Qt_old, Qinvt_old] = gen_data_sdc(T, A_old, d_old, k_old, Q00); %initialize

mat_Q_rej = zeros(nsim,T);
for sind = 1:nsim
    
%     % draw data given Pt_old
%     zt = gen_data_zt(Pt_old);
    
%     % draw Qt,Pt given zt and Ainv (MH)
%     [Qt_old, Qinvt_old, Q_rej] = gen_Qt(Qt_old, Qinvt_old, d_old, k_old, Ainv_old, Q00, zt);
%     for t=1:T
%         Pt_old(:,:,t) = QtoP(Qt_old(:,:,t));
%     end
    
%     % draw Ainv given data
%     Ainv_old = gen_Ainv(d_old, k_old, Q00, Qt_old, Qinvt_old, prior_Ainv);
%     A_old = eye(m)/Ainv_old;
    
    % draw data given d_old
    [Pt_old, Qt_old, Qinvt_old] = gen_data_sdc(T, A_old, d_old, k_old, Q00);
    
    % drawing d
    [d_old, rej_d, prior_d] = gen_d(Qt_old, Qinvt_old, d_old, k_old, Ainv_old, Q00, prior_d, sind);
    
    % store correlation and A
    m2_Pt(sind,:) = (squeeze(Pt_old(1,2,:)));
    m2_A(:,:,sind) = A_old;
%     mat_Q_rej(sind,:) = Q_rej;
    
    m2_d_rej(sind,:) = rej_d;
    m2_d(sind,1) = d_old;
end

mean(mat_Q_rej,1)

%% report
disp('method 1, (A), mean and std')
mean(m1_A,3)
std(m1_A,[], 3)

disp('method 2, (A), mean and std')
mean(m2_A,3)
std(m2_A,[], 3)


%% Histograms
figure(1)
hist([squeeze(m1_A(1,1,:)), squeeze(m2_A(1,1,:))],100)
figure(2)
hist([squeeze(m1_A(1,2,:)), squeeze(m2_A(1,2,:))],100)
figure(3)
hist([squeeze(m1_A(2,2,:)), squeeze(m2_A(2,2,:))],100)


figure(101)
hist([m1_Pt(:,1), m2_Pt(:,1)],100)
figure(102)
hist([m1_Pt(:,2), m2_Pt(:,2)],100)
figure(103)
hist([m1_Pt(:,3), m2_Pt(:,3)],100)
% figure(104)
% hist([m1_Pt(:,4), m2_Pt(:,4)],100)
% figure(105)
% hist([m1_Pt(:,5), m2_Pt(:,5)],100)

figure(201)
hist([m1_d, m2_d],100);
mean(m2_d_rej)
%% Figure
% plot(squeeze(Pt(1,2,:)))
% q1 = quantile(mat_Pt, 0.1, 1)';
% q2 = quantile(mat_Pt, 0.9, 1)';
% q3 = quantile(mat_Pt, 0.5, 1)';
% hold on
% plot(q1, 'r');
% plot(q2, 'r');
% plot(q3, 'g');
% hold off



