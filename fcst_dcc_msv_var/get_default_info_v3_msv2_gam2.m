function info = get_default_info_v3_msv2_gam2()

% --------------------------
% Initialize info structure
info = [];

% --------------------------
% MCMC SETTING
info.ndraws  = 1000; % number of draws, in a typical application, most likely it needs to be >50000
info.nburn   = 100;  %
info.nthin   = 5;    % we store every {nthin} draws

% info.ndraws  = 50000;
% info.nburn   = 5000;
% info.ndraws  = 100000;
% info.nburn   = 10000;
% info.nthin   = 5;
% info.nreport = info.ndraws + info.nburn;

info.nreport = 100;

info.save_pred_dens = 0;

% ----
% VAR FORM
info.p   = 2; % number of lags
info.nex = 1; % if a constant is included

% ----
% FORECASTING
info.hmax = 8; % maximum forecast horizon
info.stfcst  = 1; % Impose stationarity when computing the forecast: 0 for unrestricted forecast


% temp
% info.eval_T0 = '3/1/2016'; %starting of evaluation period
% info.eval_T1 = '12/1/2016'; %end of evaluation period

% original
info.eval_T0 = '3/1/1987'; %starting of evaluation period
info.eval_T1 = '12/1/2016'; %end of evaluation period

% past
% info.eval_T1 = '3/1/2007'; %end of evaluation period
% info.eval_T1 = '3/1/16'; %end of evaluation period
% info.eval_T1 = '12/1/10'; %end of evaluation period
% 3,6,9,12

% Primiceri
info.primiceri = 3; % 1=original Primiceri model; 2=TVP with Dynamic correlation, 3=TVP with Hansen model, 0=DLM


% --------------------------
% PRIOR 
if info.primiceri == 0 %DLM
% --------------------------
% PRIOR
info.T0           = 40;      % number of training sample
info.kkappa       = 4.0;     % scales C00
info.ggamma       = 1.0/3.0; % scales S10
info.ddelta       = 0.01;    % scales S0
info.ssigma_bbeta = 0.01;    % standard deviation of the proposal density for beta

% Prior for the Beta Distribution. Use a high value for a. Lower Values can generate values for H(t) that complicate chol(H(t)^-1).
info.a = 323.3333;
info.b = 30;

elseif info.primiceri == 1 % Original Primiceri
    info.T0           = 40;      % number of training sample
    
elseif info.primiceri == 2 % Primiceri with dynamic-stochastic correlation
    info.T0           = 40;      % number of training sample
    
elseif info.primiceri == 3 % Primiceri with dynamic-stochastic correlation (Hansen)
    info.T0           = 40;      % number of training sample
end





