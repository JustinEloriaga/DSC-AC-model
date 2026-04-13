close all
clear all

load data

% Primiceri and Delnegro's original setup
y=data;
lags=2;
T0=40;
T0B=40;
T0A=[2 3];
T0H=4;
kB=.01;
kA=.1;
kH=.01;
M=500;
N=100; %#of burns

% ---
% jonas/minchul's adapation to do forecasting exercise
Yest = y(1:end-8,:);
Yact = y(end-7:end,:);
% ---

r=tvsvar(Yest,lags,T0,T0B,T0A,T0H,kB,kA,kH,M);

% int=10; % thinning parameter (choose 1 if want to use all draws for constructing 
        % the final graphs; choose J if want to use one every J draws for 
        % constructing the final graphs) 
% NewGraphs;

% ----
% Primiceri and Del Negro's code that Minchul is not sure about ...
% - I don't know why we put "-x" rather than "+x" in "tria(x)"
% ----

%% Forecast 

% function for the forecasting similar to our "fcst_var_tvpsv.m"
% output of this function is exactly the same as "fcst_var_tvpsv.m"

% r is a structure with
%     B: [153×21×600 double]
%     V: [27×27×601 double]
%     A: [153×3×601 double]
%     H: [153×3×601 double]

% other inputs
info.lags = lags;
info.nburns = N;

hmax = 8;
Xcalest = [];

rst = fcst_var_primiceri(r, hmax, Yest, Yact, Xcalest, info);















