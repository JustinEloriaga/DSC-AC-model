function report=dsc_return_ar1(panel)
%DSC_RETURN_AR1 Test weekly return predictability with an intercept and HAC SE.
% Observations are percent log returns on the supplied weekly date grid.
% Missing observations exclude both adjacent regression pairs. Scores stay
% on the calendar grid, including zero rows at invalid pairs, for HAC lags.
required={'dates','returns','observation_mask','tickers'};
if ~isstruct(panel)||~isscalar(panel)||~all(isfield(panel,required))
    error('dsc:ReturnAR1Panel','A dated return panel with observation mask and tickers is required.');
end
dates=panel.dates(:); values=panel.returns; mask=panel.observation_mask;
if ~isdatetime(dates)||isempty(dates)||any(isnat(dates))|| ...
        any(diff(dates)<=seconds(0))|| ...
        ~isnumeric(values)||~isreal(values)||ndims(values)>2|| ...
        size(values,1)~=numel(dates)||size(values,2)<1
    error('dsc:ReturnAR1Panel','Dates must increase strictly and match a real return matrix.');
end
if ~isequal(size(mask),size(values))|| ...
        ~(islogical(mask)||(isnumeric(mask)&&isreal(mask)&&all(ismember(mask(:),[0 1]))))|| ...
        any(isinf(values(:)))||any(~isfinite(values(logical(mask))))
    error('dsc:ReturnAR1Panel','The observation mask must match returns and all observed values must be finite.');
end
if ~(isstring(panel.tickers)||iscellstr(panel.tickers))||~isvector(panel.tickers)
    error('dsc:ReturnAR1Panel','Tickers must identify every return column.');
end
tickers=string(panel.tickers(:)); m=size(values,2);
if numel(tickers)~=m||any(ismissing(tickers))||any(strlength(tickers)==0)|| ...
        numel(unique(tickers))~=m
    error('dsc:ReturnAR1Panel','Tickers must be nonempty, unique and match the return columns.');
end
% Accept omitted calendar weeks without treating the next observation as an
% adjacent week. The production panel already has a complete weekly grid.
week=round(days(dates-dates(1))/7);
if any(dates~=dates(1)+calweeks(week))||any(diff(week)<1)
    error('dsc:ReturnAR1Dates','All return dates must lie on one weekly calendar grid.');
end
calendar_length=week(end)+1;
calendar_values=nan(calendar_length,m); observed=false(calendar_length,m);
calendar_values(week+1,:)=double(values);
observed(week+1,:)=logical(mask)&isfinite(values);
lag=8;
empty=struct('ticker','','n_observations',0,'intercept',NaN,'rho',NaN, ...
    'standard_error',NaN,'statistic',NaN,'p_value',NaN, ...
    'status','unavailable','reason','');
rows=repmat(empty,m,1);
for j=1:m
    row=empty; row.ticker=char(tickers(j));
    valid=observed(2:end,j)&observed(1:end-1,j);
    N=nnz(valid); row.n_observations=N;
    if N<3
        row.reason='Fewer than three observed consecutive-week pairs.';
        rows(j)=row; continue
    end
    previous=calendar_values(1:end-1,j); current=calendar_values(2:end,j);
    X=[ones(N,1),previous(valid)]; y=current(valid);
    [Q,R]=qr(X,0);
    if rank(R)<2||rcond(R)<1e-12
        row.reason='The intercept and lagged return regressors are rank deficient or nearly constant.';
        rows(j)=row; continue
    end
    coefficients=R\(Q'*y); residual=y-X*coefficients;
    row.intercept=coefficients(1); row.rho=coefficients(2);
    scores=zeros(calendar_length-1,2); scores(valid,:)=X.*residual;
    meat=scores'*scores;
    for k=1:min(lag,size(scores,1)-1)
        cross=scores(k+1:end,:)'*scores(1:end-k,:);
        meat=meat+(1-k/(lag+1))*(cross+cross');
    end
    bread=R\(R'\eye(2));
    covariance=(N/(N-2))*(bread*meat*bread);
    variance=covariance(2,2);
    if ~all(isfinite(coefficients))||~isfinite(variance)||variance<=0
        row.reason='The HAC slope variance is zero or nonfinite.';
        rows(j)=row; continue
    end
    row.standard_error=sqrt(variance);
    row.statistic=row.rho/row.standard_error;
    row.p_value=erfc(abs(row.statistic)/sqrt(2));
    row.status='available'; rows(j)=row;
end
report_rows=repmat(struct('ticker','','p_value',NaN),m,1);
for j=1:m
    report_rows(j).ticker=rows(j).ticker;
    report_rows(j).p_value=rows(j).p_value;
end
report=struct('schema','return_ar1_v1', ...
    'method','OLS with intercept; Newey-West HAC covariance with Bartlett weights.', ...
    'equation','x_t = a + rho*x_(t-1) + e_t', ...
    'null_hypothesis','rho = 0', ...
    'rows',report_rows);
end
