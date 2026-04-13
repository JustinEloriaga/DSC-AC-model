%% Updates
% diary output_window

% ===
% _msv: dcc-based msv model (2021/3/24)
% _msv2: dcc-based msv model with Hansen's parameterization (2021/4/5)

% ===

% -------------------------------
% Minchul's modification in 2020/2/3
% main_estimation_and_forecast.m

% main_estimation_and_forecast_v2.m : Minchul's adaptation to allow parfor
% over ordering and time (3 variable VAR)

% main_estimation_and_forecast_v3.m : same as v2 but four variables


% -------------------------------


% date: 12/4/2019
% fcst_var_tvpsv.m
% -	Allow to compute pred score for h>1
% -	Modify to count non-stationary draws of B correctly
% -	Allow to compute pred score for inflation rate (second variable) as well

% fcst_var_flat.m
% -	now we have pred score estimate

% get_default_info.m
% -	Allow us to specify beginning and ending evaluation sample

% main_estimation_and_forecast
% - one thing to check: some variables do not have few observations at the
% beginning of the sample


%% housekeeping
clc; clear variables; close all;
workpath = pwd;
datapath = pwd;
cd(workpath);
savepath = [workpath, filesep, 'pred_v4_msv2_gam2'];
rng('default'); % reinitialize the random number generator to its startup configuration
rng(83);         % set seed

addpath([workpath,filesep, 'fcst_Primiceri_v00'])
addpath([workpath,filesep, 'toolbox_msv'])
addpath([workpath,filesep, 'toolbox_msv_corr'])


%% Get default info common to all estimation and prediction
info = get_default_info_v3_msv2_gam2(); % initialize and put default values

% patch extra path info
info.workpath = workpath;
info.savepath = savepath;
info.datapath = datapath;

%% Setup - VARIABLE LIST

% _v3 ONLY ALLOWS UNRATE ( ** 4-variable VAR ** )
var_list = {'UNRATE'};

%% Collecting all possible combinations


t_set = (1:1:120)';

all_set = [ ...
    kron(ones(1,numel(m_set)),t_set'); ...
    kron(m_set', ones(1,numel(t_set)))];
n_all_set = size(all_set,2);


%% TEST

% tStart   = tic;  
% ind_m = 4;
% ind_t = 9;
% main_estimation_and_forecast_v3c(ind_m,var_list{1},info,[],ind_t);
% tElapsed = toc(tStart); 
% disp(['Total time: ', num2str(tElapsed), ' sec ...'])
    
%% Actual computation

tStart   = tic;  

% Pararell session

% delete(gcp('nocreate'))
% parpool(16)

% Actual Loop
ind_m1 = 1;

parfor ind_t=1:30

        disp(['Starting ... t = ', num2str(ind_t), ' ... m = ', num2str(ind_m)]);
        main_estimation_and_forecast_v3c_msv_gam2(ind_m,var_list{1},info,[],ind_t);
        disp(['Done ... t = ', num2str(ind_t), ' ... m = ', num2str(ind_m)]);
    
end

tElapsed = toc(tStart); 
disp(['Total time: ', num2str(tElapsed), ' sec ...'])





