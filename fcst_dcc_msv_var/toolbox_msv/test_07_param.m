
%% 
clc; clear all; close all;

T= 200;


%% generate data 
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

%% PRIOR as in Primiceri
% W ~ IW(kw^2*4*I(m), 4)
% T0H: degrees of freedom for the IW prior distribution of var(errH)
% kH: constant scaling the prior of var(errH)
T0H=4;
kH=.01;
% W_old = iwishrnd(kH^2*T0H*eye(m), T0H);


%% Posterior sampling
nsim = 10000;
mat_W = zeros(m,m,nsim);
for sind = 1:nsim
    % Drawing var(errH)
    errH = ht(2:end,:)-ht(1:end-1,:);
    V1H = errH'*errH + (kH^2)*eye(m)*T0H;
    W_old = iwishrnd(V1H, T-1+T0H);
    
    % Store
    mat_W(:,:,sind) = W_old;
end

disp( mean(mat_W,3) );








