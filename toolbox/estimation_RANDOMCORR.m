function estimation_RANDOMCORR(info)
% TVP-VAR with Hansen's MSV2 model for MARS portfolio data
% 3-variable VAR: [equities, bonds, commodities]

%% Load data
try
    tmp = readtable([info.datapath, '/../data.csv']);
catch
    tmp = readtable('data.csv');
end

dates_raw = datetime(tmp.obs_date);
if isfield(info, 'include_inflation') && info.include_inflation
    data_raw = [tmp.mars_equities_portfolio, ...
                tmp.mars_bonds_portfolio, ...
                tmp.mars_commodities_portfolio, ...
                tmp.mars_inflation_portfolio];
else
    data_raw = [tmp.mars_equities_portfolio, ...
                tmp.mars_bonds_portfolio, ...
                tmp.mars_commodities_portfolio];
end

valid_idx   = all(~isnan(data_raw), 2);
data_daily  = data_raw(valid_idx, :);
dates_daily = dates_raw(valid_idx);

%% Aggregate daily returns to model frequency
n_vars = size(data_daily, 2);

switch lower(info.freq)
    case 'monthly'
        grp = year(dates_daily)*100 + month(dates_daily);
    case 'weekly'
        grp = year(dates_daily)*100 + week(dates_daily, 'weekofyear');
    case 'daily'
        grp = (1:numel(dates_daily))';
    otherwise
        error('info.freq must be ''monthly'', ''weekly'', or ''daily'' (got ''%s'')', info.freq);
end

[~, ~, ic] = unique(grp);
n_periods  = max(ic);
data       = zeros(n_periods, n_vars);
dcal       = cell(n_periods, 1);

for pp = 1:n_periods
    idx          = (ic == pp);
    data(pp, :)  = sum(data_daily(idx, :), 1);
    period_dates = dates_daily(idx);
    dcal{pp}     = datestr(period_dates(end), 'yyyy-mm-dd');
end
data = data * 100;

info.T0 = round(0.1 * n_periods);

disp('=== Data summary ===');
disp(['  Freq    : ', info.freq]);
disp(['  Obs     : ', num2str(n_periods), '  (', dcal{1}, ' to ', dcal{end}, ')']);
disp(['  T0      : ', num2str(info.T0), ' periods (10% of sample)']);
if isfield(info, 'include_inflation') && info.include_inflation
    disp(['  Vars    : equities, bonds, commodities, inflation swaps']);
else
    disp(['  Vars    : equities, bonds, commodities']);
end
disp('====================');

%% Build VAR matrices
[~, Y] = make_varXY(data, info.p, info.nex);
Xcal   = dcal(info.p+1:end, :);

%% Resolve end of estimation sample
eval_T1 = info.eval_T1;
if datenum(eval_T1, 'yyyy-mm-dd') > datenum(Xcal{end}, 'yyyy-mm-dd')
    eval_T1 = Xcal{end};
end

i1 = find(strcmp(Xcal, eval_T1));
if isempty(i1)
    dates_num = datenum(Xcal, 'yyyy-mm-dd');
    i1 = find(dates_num <= datenum(eval_T1, 'yyyy-mm-dd'), 1, 'last');
end

Yest    = Y(1:i1, :);
Xcalest = Xcal(info.T0+1:i1, :);

%% Run MCMC
M = info.ndraws + info.nburn;
r = tvsvar_modified_msv2_gam2_gen(Yest, info.p, info.T0, info.kB, M, info.nreport, info.nthin);

%% Save
strname      = strrep(Xcalest{end}, '/', '-');
savefilename = ['RANDOMCORR_', info.freq, '_', strname, '.mat'];
save(fullfile(info.workpath, savefilename), 'r');
disp(['Saved: ', savefilename]);
end
