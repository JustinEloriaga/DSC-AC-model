function [return_flag,log_plike_mean_vec,log_plike_11_mean_vec,log_plike_rao_mean_vec,log_plike_rao_11_mean_vec] = main_estimation_and_forecast_v3c_msv_gam2(var_ordering,var_name,info,j_var,ind_t)

% to run minchul's job 2020/1/31
% ind_t is determined outside of the function

% ===
% _msv: dcc-based msv model (2021/3/24)
% ===

% _v3: 4-variable VAR

% _v3c: tvsvar_modified.m is used

% --------------------------
%% Housekeeping
% --------------------------
% Function to estimate TVP-VAR-SV
% Code based on main_minchul3_loop_prototype.m
% com2: [base three variables + extra variable]

% Updates (2019/12/3)
% 1) minchul modified transformation part of the code
% -> there is an issue about unit for the extra variable
% 2) moved common setting (info) part of the code to get_default_info.m
% 3) we introduce info.eval_T0, info.eval_T1: starting and end points of evaluation

% --------------------------
%% Data selection
% --------------------------

% QD-FRED dataset: Next eight lines follow QD-FRED matlab code
% Load data from csv file

switch var_name
    case 'OILPRICEx'
        tmp     = importdata([info.datapath,'/data/fred-qd/1974-1-2019-11.csv'],',');
    case 'TWEXMMTH'
        tmp     = importdata([info.datapath,'/data/fred-qd/1973-3-2019-11.csv'],',');
    otherwise
        try
            tmp     = importdata([info.datapath,'/data/fred-qd/2019-11.csv'],',');
        catch
            tmp     = importdata(['data/fred-qd/2019-11.csv'],',');
        end
end

% Variable names
Xnames  = tmp.textdata(1,2:end);
% Transformation numbers
tcode   = tmp.data(2,:);
% Raw data
Xdata   = tmp.data(3:end,:);
% Get calendar info
Xcal    = tmp.textdata(4:end,1);

% Selecting variables var_index = find(strcmp(series,'Variable Name')==1);.
% Examples: 1 = find(strcmp(series,'GDPC1')==1);
%          97 = find(strcmp(series,'GDPCTPI')==1);
% Variables mnenomics are obtained from FRED-QD_appendix.pdf: https://s3.amazonaws.com/files.fred.stlouisfed.org/fred-md/FRED-QD_appendix.pdf


% with growth rate (Basic three variables)
% data = [diff(log(Xdata(:,strcmp(Xnames,'GDPC1')==1)))*400, diff(log(Xdata(:,strcmp(Xnames,'PCEPILFE')==1)))*400, Xdata(2:end,strcmp(Xnames,'TB3MS')==1) ];

data = [diff(log(Xdata(:,strcmp(Xnames,'GDPC1')==1)))*400, ...
    diff(log(Xdata(:,strcmp(Xnames,'PCEPILFE')==1)))*400, ...
    Xdata(2:end,strcmp(Xnames,'TB3MS')==1), ...
    Xdata(2:end,strcmp(Xnames,'UNRATE')==1) ];

% calendar (first observation is gone because of growth)
dcal = Xcal(2:end,1);

% NOTE: we start from 3 because some transformation requires to use first two obs
% -> we cut one more initial obs
data(1,:) = [];
dcal(1,:) = [];

if strcmp(var_name,'UNRATE')==0
    disp(['This function only allows UNRATE ... Foreced to quit']);
    adklfjldjflkd;
end

% Augmenting extra variable
% if strcmp(var_name,'Baseline')==0
%     
%     % Transformation
%     t_var_name = tcode(strcmp(Xnames,var_name)==1); %transformation code
%     temp_x = Xdata(:,strcmp(Xnames,var_name)==1);  %raw data
%     
%     % some observations are not avariable
%     if sum(isnan(temp_x))>0
%         disp('Line 59 Main Estimation and Forecast');
%         j_var
%         return
%     end
%     
%     temp_x = TransData(temp_x, t_var_name); %transform it
%     
%     % Standardize it - probably it is ok because we are not interested inforecasting this variable
%     temp_x = (temp_x - mean(temp_x(3:end,:)))/std(temp_x(3:end,:)) * std(data(:,1)); %to have the same level of standard deviation of output
%     %temp_x = temp_x*100; %JA 12/03/2019 replicate error with smoother_mc
%     
%     % Augment a new series
%     data   = [data, temp_x(3:end,:)]; %it is better to start from 3 because some transformation requires to use first two obs
%     
% end

% OLD way (we may come back to this)
% if strcmp(var_name,'Baseline')==0
%
%     t_var_name = tcode(strcmp(Xnames,var_name)==1);
%     % ---------
%     % Augment additional variable
%     switch t_var_name
%         case 1
%             data = [data, Xdata(2:end,strcmp(Xnames,var_name)==1) ];                % no transformation
%         case 5
%             data = [data, 100*diff(log(Xdata(1:end,strcmp(Xnames,var_name)==1))) ]; % log-differences in percent
%         otherwise
%             disp('NO extra data series ... stop the program ');
%             data = [];
%     end
%
% end

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

        % Ordering for the four variable VAR
        ord_set = flipud(perms([1,2,3,4]));
        n_ord_set = size(ord_set,1);
        tmpYest = Yest(:,ord_set(var_ordering,:));
        tmpYact = Yact(:,ord_set(var_ordering,:));

        Yest = tmpYest;
        Yact = tmpYact;
        % Primiceri and Delnegro's original setup
        %y=Yest;
        lags=2;
        T0=40;
        T0B=40;
        
        
        %T0A=[2 3];
        T0A = [2 3 4]; % to adapt the 4-VAR VAR ( dim+1 is default prior)
        
        
        T0H=4;
        kB=.01;
        kA=.1;
        kH=.01;
        N=info.nburn; %#of burns
        M=info.ndraws + info.nburn;
        
%         % ---
%         % jonas/minchul's adapation to do forecasting exercise
%         Yest = y(1:end-8,:);
%         Yact = y(end-7:end,:);
%         % ---
        
        %r =tvsvar(Yest,lags,T0,T0B,T0A,T0H,kB,kA,kH,M);
        r =tvsvar_modified(Yest,lags,T0,T0B,T0A,T0H,kB,kA,kH,M);
        
        % other inputs
        info.lags = lags;
        info.nburns = N;
        
        hmax = 8;
        %Xcalest = [];
        %rp = fcst_var_primiceri(r, hmax, Yest, Yact, Xcalest, info);
        rp = fcst_var_primiceri(r, hmax, Yest, Yact, [], info);

        %% Save
        strname = Xcalest(end,:);
        strname = strrep(strname{1},'/','-');
        savefilename = ['primiceri_long_',num2str(var_ordering),'_pred_v',var_name, '_t_',strname, '.mat'];
        cd(info.savepath);
        %save(savefilename, 'rst', 'rst_flat_p2','rate_bbeta');
        save(savefilename, 'rp');
        cd(info.workpath);
  
%         log_plike_mean_vec(ind_t,1)    = rp.log_plike_mean;
%         log_plike_11_mean_vec(ind_t,1) = rp.log_plike_11_mean;
    
  %      log_plike_rao_mean_vec(ind_t,1)    = rp.log_plike_rao_mean(1);
  %      log_plike_rao_11_mean_vec(ind_t,1) = rp.log_plike_rao_11_mean(1);
        
        
  %=============================================================================
    elseif info.primiceri==2 %dynamic-correlation based model

        % Ordering for the four variable VAR
        ord_set = flipud(perms([1,2,3,4]));
        n_ord_set = size(ord_set,1);
        tmpYest = Yest(:,ord_set(var_ordering,:));
        tmpYact = Yact(:,ord_set(var_ordering,:));

        Yest = tmpYest;
        Yact = tmpYact;
        
        % Primiceri and Delnegro's original setup
        %y=Yest;
        lags=2;
        T0=40;
        T0B=40;
        
        %T0A=[2 3];
        T0A = [2 3 4]; % to adapt the 4-VAR VAR ( dim+1 is default prior)
        
        T0H=4;
        kB=.01;
        kA=.1;
        kH=.01;
        N=info.nburn; %#of burns
        M=info.ndraws + info.nburn;
        
        %r =tvsvar(Yest,lags,T0,T0B,T0A,T0H,kB,kA,kH,M);
        nreport = info.nreport;
        
        if info.prior_original == 1
            r =tvsvar_modified_msv_gam2(Yest,lags,T0,T0B,T0A,T0H,kB,kA,kH,M, nreport, info.nthin);
            % other inputs
            info.lags = lags;
            info.nburns = N;
            
            hmax = 8;
            %Xcalest = [];
            %rp = fcst_var_primiceri(r, hmax, Yest, Yact, Xcalest, info);
            rp = fcst_var_primiceri_msv(r, hmax, Yest, Yact, [], info);
            
            %% Save
            strname = Xcalest(end,:);
            strname = strrep(strname{1},'/','-');
            savefilename = ['msv_gam2_long_',num2str(var_ordering),'_pred_v',var_name, '_t_',strname, '.mat'];
            cd(info.savepath);
            save(savefilename, 'rp');
            cd(info.workpath);
        else
            r =tvsvar_modified_msv(Yest,lags,T0,T0B,T0A,T0H,kB,kA,kH,M, nreport);
            % other inputs
            info.lags = lags;
            info.nburns = N;
            
            hmax = 8;
            %Xcalest = [];
            %rp = fcst_var_primiceri(r, hmax, Yest, Yact, Xcalest, info);
            rp = fcst_var_primiceri_msv(r, hmax, Yest, Yact, [], info);
            
            %% Save
            strname = Xcalest(end,:);
            strname = strrep(strname{1},'/','-');
            savefilename = ['msv_long_',num2str(var_ordering),'_pred_v',var_name, '_t_',strname, '.mat'];
            cd(info.savepath);
            save(savefilename, 'rp');
            cd(info.workpath);
        end
        
        
  
        
  %=============================================================================
    elseif info.primiceri==3 %dynamic-correlation based model, (hansen's parameterization)

        % Ordering for the four variable VAR
        ord_set = flipud(perms([1,2,3,4]));
        n_ord_set = size(ord_set,1);
        tmpYest = Yest(:,ord_set(var_ordering,:));
        tmpYact = Yact(:,ord_set(var_ordering,:));

        Yest = tmpYest;
        Yact = tmpYact;
        
        % Primiceri and Delnegro's original setup
        %y=Yest;
        lags=2;
        T0=40;
        T0B=40;
        
        %T0A=[2 3];
        T0A = [2 3 4]; % to adapt the 4-VAR VAR ( dim+1 is default prior)
        
        T0H=4;
        kB=.01;
        kA=.1;
        kH=.01;
        N=info.nburn; %#of burns
        M=info.ndraws + info.nburn;
        
        %r =tvsvar(Yest,lags,T0,T0B,T0A,T0H,kB,kA,kH,M);
        nreport = info.nreport;

        %r =tvsvar_modified_msv(Yest,lags,T0,T0B,T0A,T0H,kB,kA,kH,M, nreport);
        %r =tvsvar_modified_msv2(Yest,lags,T0,T0B,T0A,T0H,kB,kA,kH,M, nreport);
        r =tvsvar_modified_msv2_gam2(Yest,lags,T0,T0B,T0A,T0H,kB,kA,kH,M, nreport , info.nthin);
        
        % other inputs
        info.lags = lags;
        info.nburns = N;
        
        hmax = 8;
        %Xcalest = [];
        %rp = fcst_var_primiceri(r, hmax, Yest, Yact, Xcalest, info);
        % rp = fcst_var_primiceri_msv(r, hmax, Yest, Yact, [], info);
        rp = fcst_var_primiceri_msv2(r, hmax, Yest, Yact, [], info);

        %% Save
        strname = Xcalest(end,:);
        strname = strrep(strname{1},'/','-');
        savefilename = ['msv2_long_',num2str(var_ordering),'_pred_v',var_name, '_t_',strname, '.mat'];
        cd(info.savepath);
        save(savefilename, 'rp' ,'r');
        cd(info.workpath);
        
  %=============================================================================
    else
    

    
    %% TVP-SV
    % Estimation - tvpsv
    [gb_B,~,gb_H,gb_bbeta,gb_W,rate_bbeta,log_plike_mean,log_plike_11_mean] = gb_var_tvpsv(Yest,Xest,Yact,info,j_var);
    %[gb_B,~,gb_H,gb_bbeta,gb_W,rate_bbeta,plike_mean,plike_11_mean] = gb_var_tvpsv_mc_v03(Yest,Xest,Yact,info);
    
    % Forecasting - tvpsv
    rst = fcst_var_tvpsv(gb_B,gb_H,gb_bbeta,gb_W,hmax, Yest, Xest, Yact,Xcalest,info);
    
%     %% Const-VAR
%     % Estimation - var-flat
%     [flat_PHI, flat_SIG] = gb_var_flat(Yest,Xest,t,info);
%     % Forecasting - var-flat
%     rst_flat_p2          = fcst_var_flat(flat_PHI,flat_SIG,hmax,info.nex,Yest,Xest,Yact, Xcalest,info);
    
    %% Save
    strname = Xcalest(end,:);
    strname = strrep(strname{1},'/','-');
    savefilename = ['pred_v',var_name, '_t_',strname, '.mat'];
    cd(info.savepath);
    %save(savefilename, 'rst', 'rst_flat_p2','rate_bbeta');
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


