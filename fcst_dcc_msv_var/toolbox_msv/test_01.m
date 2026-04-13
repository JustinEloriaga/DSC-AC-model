% Wishart based volatility model from Asai and McAleer (Joe, 2009)

% Geweke test for Ainv

clc; clear all; close all;

%% Unknowns
% A, d, k

A     = [30, -3; -3, 30];
% A = inv([1, 0.3; 0.3, 1]);
d     = 0.8;
k     = 50;


%% Data generation
T = 5;
m = size(A,1);
Q00      = eye(m); %we assume that this is known
Q00(1,2) = 0.9;
Q00(2,1) = 0.9;

[Pt, Qt,Qinvt] = gen_data_sdc(T, A, d, k, Q00);

% plot(squeeze(Pt(1,2,:)))

%% Prior
% prior for Ainv (following Asai and McAleer)
prior_Ainv.gam  = m + 100;
prior_Ainv.Cinv = ([30, -3; -3, 30]) * (prior_Ainv.gam);
Ainv_new = gen_Ainv(d, k, Q00, Qt, Qinvt, prior_Ainv);


%% Geweke test, method 1
nsim = 50000;

disp('implementing method 1 ...');
m1_Pt = zeros(nsim, T);
m1_A = zeros(m,m,nsim);
for sind = 1:nsim
    % draw Ainv from prior
    Ainv_old = wishrnd(eye(m)/prior_Ainv.Cinv, prior_Ainv.gam);
    A_old = eye(m)/Ainv_old;
    
    % draw data given new A draw
    [Pt, Qt, Qinvt] = gen_data_sdc(T, A_old, d, k, Q00);
    
    % store correlation and A
    m1_Pt(sind,:) = (squeeze(Pt(1,2,:)));
    m1_A(:,:,sind) = A_old;
end

%% Geweke test, method 2
disp('implementing method 2 ...');
m2_Pt = zeros(nsim, T);
m2_A = zeros(m,m,nsim);

% initialize with prior
Ainv_old = wishrnd(eye(m)/prior_Ainv.Cinv, prior_Ainv.gam);
A_old = eye(m)/Ainv_old;

for sind = 1:nsim
    % draw data given A
    [Pt, Qt, Qinvt] = gen_data_sdc(T, A_old, d, k, Q00);
    
    % draw Ainv given data
    Ainv_old = gen_Ainv(d, k, Q00, Qt, Qinvt, prior_Ainv);
    A_old = eye(m)/Ainv_old;
    
    % store correlation and A
    m2_Pt(sind,:) = (squeeze(Pt(1,2,:)));
    m2_A(:,:,sind) = A_old;
end

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








