% ======================================================================
% Adapted from main_v3c_msv2_gam2.m
% Uses MARS portfolio data (data_JRR.csv) instead of FRED-QD macro data
% 3-variable VAR: [equities, bonds, commodities]
% (drops mars_inflation_portfolio to avoid NaN issues)
% ======================================================================

% ======================================================================
% NOTES (msv2 model with MARS data, 3-variable)
% 1. tvsvar_modified_msv2_gam2 is the main file that runs MCMC
% 2. Some of hyperparameters are controlled in the following script:
%    "get_default_info_v3_msv2_gam2_mars3.m"
% 3. It is based on Elliptical slice sampling as described in our JE paper.
% 4. file is stored under "\pred_v4_msv2_gam2_mars3"
%  - r: contains estimated parameters, to keep the file size small we only
%  store the last period's quantity for B and P
%  - rp: metrics for prediction eval
% ======================================================================


% ===
% _msv: dcc-based msv model (2021/3/24)
% _msv2: dcc-based msv model with Hansen's parameterization (2021/4/5)

% ===

% filename without "driver": single run


%% housekeeping
clc; clear variables; close all;
workpath = pwd;
datapath = pwd;
cd(workpath);
savepath = [workpath, filesep, 'pred_v4_msv2_gam2_mars3'];
rng('default'); % reinitialize the random number generator to its startup configuration
rng(83);         % set seed

addpath([workpath,filesep, 'fcst_Primiceri_v00'])
addpath([workpath,filesep, 'toolbox_msv'])
addpath([workpath,filesep, 'toolbox_msv_corr'])

% Create output directory if it does not exist
if ~exist(savepath, 'dir')
    mkdir(savepath);
end


%% Get default info common to all estimation and prediction
info = get_default_info_v3_msv2_gam2_mars3(); % initialize and put default values

% patch extra path info
info.workpath = workpath;
info.savepath = savepath;
info.datapath = datapath;

%% Setup - VARIABLE LIST

% For MARS data we use 3 portfolios as the VAR variables:
% [equities, bonds, commodities]
var_list = {'MARS3'};


%% Actual computation

tStart   = tic;

ind_t = 1; % ind_t=1 (see how ind_t is used in main_estimation_and_forecast_v3c_msv_gam2_mars3.m)
ind_m = 1; % variable ordering: does not matter for msv model

disp(['Starting ... t = ', num2str(ind_t), ' ... m = ', num2str(ind_m)]);
main_estimation_and_forecast_v3c_msv_gam2_mars3(ind_m,var_list{1},info,[],ind_t);
disp(['Done ... t = ', num2str(ind_t), ' ... m = ', num2str(ind_m)]);

tElapsed = toc(tStart);
disp(['Total time: ', num2str(tElapsed), ' sec ...'])