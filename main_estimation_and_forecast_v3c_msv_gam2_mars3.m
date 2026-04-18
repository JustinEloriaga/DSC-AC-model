function [return_flag,log_plike_mean_vec,log_plike_11_mean_vec,log_plike_rao_mean_vec,log_plike_rao_11_mean_vec] = main_estimation_and_forecast_v3c_msv_gam2_mars3(var_ordering,var_name,info,j_var,ind_t)

% Adapted from main_estimation_and_forecast_v3c_msv_gam2.m
% Uses MARS portfolio data (data_JRR.csv) instead of FRED-QD macro data
% 3-variable VAR: [equities, bonds, commodities]
% (drops mars_inflation_portfolio to avoid NaN issues)

% ind_t is determined outside of the function

% ===
% _msv: dcc-based msv model (2021/3/24)
% ===

% 3-variable VAR (equities, bonds, commodities)

% _v3c: tvsvar_modified.m is used

% --------------------------
%% Housekeeping
% --------------------------
% Function to estimate TVP-VAR-SV
% com2: [equities, bonds, commodities]

% --------------------------
%% Data selection
% --------------------------

% Load MARS portfolio data from data_JRR.csv
% The CSV is located one level above the working directory (fcst_dcc_msv_var)
try
    tmp = readtable([info.datapath, '/../data_JRR.csv']);
catch
    tmp = readtable('data_JRR.csv');
end

% Extract date column as datetime
dates_raw = datetime(tmp.obs_date);

% Extract three MARS portfolio return series (drop inflation to avoid NaNs)
% Columns: mars_equities_portfolio, mars_bonds_portfolio,
%          mars_commodities_portfolio
data_raw = [tmp.mars_equities_portfolio, ...
            tmp.mars_bonds_portfolio, ...
            tmp.mars_commodities_portfolio];

% Keep only rows where all three series are available
valid_idx = all(~isnan(data_raw), 2);
data_daily = data_raw(valid_idx, :);
dates_daily = dates_raw(valid_idx);

disp(['MARS daily data (3-var): ', num2str(size(data_daily,1)), ' observations']);
disp(['  Date range: ', datestr(dates_daily(1),'yyyy-mm-dd'), ' to ', datestr(dates_daily(end),'yyyy-mm-dd')]);

% ---------------------------------------------------------------
% Aggregate daily returns to MONTHLY frequency
% Sum daily log-returns within each month (log-returns are additive)
% This reduces ~10700 daily obs to ~500 monthly obs, making the
% TVP-VAR-MSV model feasible (it builds T x T matrices internally).
% ---------------------------------------------------------------
ym = year(dates_daily)*100 + month(dates_daily); % e.g. 198405, 198406, ...
[unique_ym, ~, ic] = unique(ym);

n_months = length(unique_ym);
n_vars   = size(data_daily, 2);
data     = zeros(n_months, n_vars);
dcal     = cell(n_months, 1);

for mm = 1:n_months
    idx = (ic == mm);
    % sum daily log-returns to get monthly return
    data(mm,:) = sum(data_daily(idx,:), 1);
    % use last trading day of the month as the date label
    month_dates = dates_daily(idx);
    dcal{mm} = datestr(month_dates(end), 'yyyy-mm-dd');
end

disp(['MARS monthly data (3-var): ', num2str(size(data,1)), ' observations']);
disp(['  Date range: ', dcal{1}, ' to ', dcal{end}]);

% Scale returns to percentage units (multiply by 100) to be
% comparable in magnitude to the macro variables in the original code
data = data * 100;

% NOTE: no differencing needed since data is already in returns form

% ----
%% Estimation part
% ----
% Transformed data as input (p, nex, X, Y have to be consistent)
[X,Y]  = make_varXY(data,info.p,info.nex);
Xcal   = dcal(info.p+1:end,:);




%% Extra info

hmax = info.hmax;          % forecast horizon

% Picking evaluation period
i0 = find(strcmp(Xcal,info.eval_T0));
i1 = find(strcmp(Xcal,info.eval_T1));

% If exact date not found, find the nearest date after eval_T0 / before eval_T1
if isempty(i0)
    dates_num = datenum(Xcal, 'yyyy-mm-dd');
    t0_num = datenum(info.eval_T0, 'yyyy-mm-dd');
    idx = find(dates_num >= t0_num, 1, 'first');
    if ~isempty(idx)
        i0 = idx;
        disp(['eval_T0 not found exactly, using nearest: ', Xcal{i0}]);
    else
        error('eval_T0 date not found in calendar');
    end
end

if isempty(i1)
    dates_num = datenum(Xcal, 'yyyy-mm-dd');
    t1_num = datenum(info.eval_T1, 'yyyy-mm-dd');
    idx = find(dates_num <= t1_num, 1, 'last');
    if ~isempty(idx)
        i1 = idx;
        disp(['eval_T1 not found exactly, using nearest: ', Xcal{i1}]);
    else
        error('eval_T1 date not found in calendar');
    end
end

% End point of evaluation period
Tmin = i0;
Tmax = min(i1,size(Y,1)-hmax);

% ---
% Choose date
set_t = Tmin:1:Tmax;
n_t   = numel(set_t);

% keyboard

log_plike_mean_vec        = nan(n_t,1);
log_plike_11_mean_vec     = nan(n_t,1);
log_plike_rao_mean_vec    = nan(n_t,1);
log_plike_rao_11_mean_vec = nan(n_t,1);

%% Looping
% for ind_t = 1:1:n_t
%     ind_t
    % Pick a date
    t       = set_t(ind_t);
        %% Sample selection
    Yest = Y(1:t,:);
    Xest = X(1:t,:);

    Yact = Y(t+1:t+hmax,:);
    Xcalest = Xcal(info.T0+1:t,:);

    if info.primiceri==1

        % Ordering for the three variable VAR
        ord_set = flipud(perms([1,2,3]));
        n_ord_set = size(ord_set,1);
        tmpYest = Yest(:,ord_set(var_ordering,:));
        tmpYact = Yact(:,ord_set(var_ordering,:));

        Yest = tmpYest;
        Yact = tmpYact;
        % Primiceri and Delnegro's original setup
        lags=2;
        T0=40;
        T0B=40;

        T0A = [2 3]; % 3-variable VAR (dim+1 is default prior)

        T0H=4;
        kB=.01;
        kA=.1;
        kH=.01;
        N=info.nburn; %#of burns
        M=info.ndraws + info.nburn;

        r =tvsvar_modified(Yest,lags,T0,T0B,T0A,T0H,kB,kA,kH,M);

        % other inputs
        info.lags = lags;
        info.nburns = N;

        hmax = 8;
        rp = fcst_var_primiceri(r, hmax, Yest, Yact, [], info);

        %% Save
        strname = Xcalest(end,:);
        strname = strrep(strname{1},'/','-');
        savefilename = ['primiceri_mars3_long_',num2str(var_ordering),'_pred_v',var_name, '_t_',strname, '.mat'];
        cd(info.savepath);
        save(savefilename, 'rp');
        cd(info.workpath);

  %=============================================================================
    elseif info.primiceri==2 %dynamic-correlation based model

        % Ordering for the three variable VAR
        ord_set = flipud(perms([1,2,3]));
        n_ord_set = size(ord_set,1);
        tmpYest = Yest(:,ord_set(var_ordering,:));
        tmpYact = Yact(:,ord_set(var_ordering,:));

        Yest = tmpYest;
        Yact = tmpYact;

        lags=2;
        T0=40;
        T0B=40;

        T0A = [2 3]; % 3-variable VAR (dim+1 is default prior)

        T0H=4;
        kB=.01;
        kA=.1;
        kH=.01;
        N=info.nburn; %#of burns
        M=info.ndraws + info.nburn;

        nreport = info.nreport;

        if info.prior_original == 1
            r =tvsvar_modified_msv_gam2(Yest,lags,T0,T0B,T0A,T0H,kB,kA,kH,M, nreport, info.nthin);
            % other inputs
            info.lags = lags;
            info.nburns = N;

            hmax = 8;
            rp = fcst_var_primiceri_msv(r, hmax, Yest, Yact, [], info);

            %% Save
            strname = Xcalest(end,:);
            strname = strrep(strname{1},'/','-');
            savefilename = ['msv_gam2_mars3_long_',num2str(var_ordering),'_pred_v',var_name, '_t_',strname, '.mat'];
            cd(info.savepath);
            save(savefilename, 'rp');
            cd(info.workpath);
        else
            r =tvsvar_modified_msv(Yest,lags,T0,T0B,T0A,T0H,kB,kA,kH,M, nreport);
            % other inputs
            info.lags = lags;
            info.nburns = N;

            hmax = 8;
            rp = fcst_var_primiceri_msv(r, hmax, Yest, Yact, [], info);

            %% Save
            strname = Xcalest(end,:);
            strname = strrep(strname{1},'/','-');
            savefilename = ['msv_mars3_long_',num2str(var_ordering),'_pred_v',var_name, '_t_',strname, '.mat'];
            cd(info.savepath);
            save(savefilename, 'rp');
            cd(info.workpath);
        end



  %=============================================================================
    elseif info.primiceri==3 %dynamic-correlation based model, (hansen's parameterization)

        % Ordering for the three variable VAR
        ord_set = flipud(perms([1,2,3]));
        n_ord_set = size(ord_set,1);
        tmpYest = Yest(:,ord_set(var_ordering,:));
        tmpYact = Yact(:,ord_set(var_ordering,:));

        Yest = tmpYest;
        Yact = tmpYact;

        lags=2;
        T0=40;
        T0B=40;

        T0A = [2 3]; % 3-variable VAR (dim+1 is default prior)

        T0H=4;
        kB=.01;
        kA=.1;
        kH=.01;
        N=info.nburn; %#of burns
        M=info.ndraws + info.nburn;

        nreport = info.nreport;

        r =tvsvar_modified_msv2_gam2_gen(Yest,lags,T0,T0B,T0A,T0H,kB,kA,kH,M, nreport , info.nthin);

        % other inputs
        info.lags = lags;
        info.nburns = N;

        hmax = 8;
        rp = fcst_var_primiceri_msv2_gen(r, hmax, Yest, Yact, [], info);

        %% Save
        strname = Xcalest(end,:);
        strname = strrep(strname{1},'/','-');
        savefilename = ['msv2_mars3_long_',num2str(var_ordering),'_pred_v',var_name, '_t_',strname, '.mat'];
        cd(info.savepath);
        save(savefilename, 'rp' ,'r');
        cd(info.workpath);

  %=============================================================================
    else


    %% TVP-SV
    % Estimation - tvpsv
    [gb_B,~,gb_H,gb_bbeta,gb_W,rate_bbeta,log_plike_mean,log_plike_11_mean] = gb_var_tvpsv(Yest,Xest,Yact,info,j_var);

    % Forecasting - tvpsv
    rst = fcst_var_tvpsv(gb_B,gb_H,gb_bbeta,gb_W,hmax, Yest, Xest, Yact,Xcalest,info);

    %% Save
    strname = Xcalest(end,:);
    strname = strrep(strname{1},'/','-');
    savefilename = ['mars3_pred_v',var_name, '_t_',strname, '.mat'];
    cd(info.savepath);
    save(savefilename, 'rst','rate_bbeta');
    cd(info.workpath);

    log_plike_mean_vec(ind_t,1)    = log_plike_mean;
    log_plike_11_mean_vec(ind_t,1) = log_plike_11_mean;

    log_plike_rao_mean_vec(ind_t,1)    = rst.plike_rao_mean(1);
    log_plike_rao_11_mean_vec(ind_t,1) = rst.plike_rao_11_mean(1);

    end

% end

% disp('Done.')
return_flag=1;

end
