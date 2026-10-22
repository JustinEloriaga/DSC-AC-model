function result = analyze_return_ar1(cfg,outdir)
%ANALYZE_RETURN_AR1 Weekly return predictability and volatility diagnostics.
%   Uses the shared panel builder, without estimation or prior recalibration.
%   Files are written only into an existing output directory. The main sample
%   excludes the initial prior_weeks; the full history is a sensitivity check.
%   AR(1): r(t)=Intercept+Phi*r(t-1)+u(t). Inference uses Bartlett Newey-West
%   HAC(8), n/(n-k) correction, and asymptotic normal reference distributions.
%   A zero Phi does not establish iid returns or test mean-state V=0.
if nargin<1 || isempty(cfg), cfg=weekly_config(); end
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'toolbox'));
if nargin<2 || isempty(outdir)
    outdir=fullfile(root,'outputs','weekly_research','data');
end
assert(isfolder(outdir),'ar1:OutputDirectory', ...
    'The output directory must already exist: %s',outdir);
assert(strcmp(cfg.freq,'weekly'),'ar1:Frequency','This diagnostic uses weekly returns.');
validateattributes(cfg.prior_weeks,{'numeric'},{'scalar','integer','nonnegative'});
baselineCfg=cfg; baselineCfg.closure_policy='asof';
panel=prepare_weekly_panel(baselineCfg);
last=numel(panel.dates);
if isfield(cfg,'estimation_end_date') && strlength(string(cfg.estimation_end_date))>0
    cutoff=datetime(cfg.estimation_end_date,'InputFormat','yyyy-MM-dd');
    assert(isscalar(cutoff)&&~isnat(cutoff)&&cutoff<=panel.dates(end), ...
        'ar1:Cutoff','The cutoff must be a valid date within the available history.');
    last=find(panel.dates<=cutoff,1,'last');
end
assert(~isempty(last)&&last-cfg.prior_weeks>=40,'ar1:ShortSample', ...
    'At least 40 modeled weeks are required after initial calibration.');
samples={cfg.prior_weeks+1:last,1:last};
sampleNames=["estimation","full_history"];
rows=struct([]); maskedRows=struct([]);
for s=1:2
    ix=samples{s};
    for j=1:numel(panel.tickers)
        row=diagnostics(panel.returns(ix,j));
        row.Sample=sampleNames(s); row.Policy="asof";
        row.Ticker=panel.tickers(j); row.Label=panel.labels(j);
        row.StartDate=string(panel.dates(ix(1)),'yyyy-MM-dd');
        row.EndDate=string(panel.dates(ix(end)),'yyyy-MM-dd');
        rows=[rows;row]; %#ok<AGROW>
        if panel.tickers(j)=="NKYTR Index"
            x=panel.returns(ix,j);
            x(panel.closure_return_mask(ix,j))=NaN;
            masked=diagnostics(x);
            masked.Sample=sampleNames(s); masked.Policy="mask";
            masked.Ticker=row.Ticker; masked.Label=row.Label;
            masked.StartDate=row.StartDate; masked.EndDate=row.EndDate;
            maskedRows=[maskedRows;masked]; %#ok<AGROW>
        end
    end
end
results=struct2table(rows);
results.PhiPHolm=nan(height(results),1);
results.AR12PHolm=nan(height(results),1);
results.ARCH12PHolm=nan(height(results),1);
results.MeanPositivePHolm=nan(height(results),1);
for s=1:2
    ix=results.Sample==sampleNames(s);
    results.PhiPHolm(ix)=holm(results.PhiP(ix));
    results.AR12PHolm(ix)=holm(results.AR12P(ix));
    results.ARCH12PHolm(ix)=holm(results.ARCH12P(ix));
    results.MeanPositivePHolm(ix)=holm(results.MeanPositiveP(ix));
end
if isempty(maskedRows), maskedResults=table(); else, maskedResults=struct2table(maskedRows); end
metadata=struct();
metadata.created_at=char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss'));
metadata.source_file=panel.source_file;
metadata.source_hash=panel.source_hash;
metadata.frequency=panel.frequency;
metadata.tickers=panel.tickers;
metadata.prior_weeks=cfg.prior_weeks;
metadata.config=baselineCfg;
metadata.estimation_start=char(string(panel.dates(cfg.prior_weeks+1),'yyyy-MM-dd'));
metadata.estimation_end=char(string(panel.dates(last),'yyyy-MM-dd'));
metadata.estimation_returns=last-cfg.prior_weeks;
metadata.hac_lags=8;
metadata.hac_sensitivity_lags=[0,4,13];
metadata.hac_method='Bartlett Newey-West; n/(n-k) finite-sample correction; normal p values and 95% confidence intervals.';
metadata.multiple_testing='Holm across the selected tickers, separately for each sample and each test family; masked NKYTR sensitivity is unadjusted.';
metadata.ar12_method='Wald test of all 12 lagged-return coefficients with HAC(8), asymptotic chi-square(12).';
metadata.arch12_method='Conventional Engle ARCH LM: n*R-squared from squared AR(1) residuals on intercept and 12 calendar lags, asymptotic chi-square(12); requires suitable moment and regularity conditions. Omitted mean dynamics can also cause squared-residual dependence.';
metadata.calendar='Masked values stay at their original dates; regressions use complete lag windows, and missing score rows are zero in HAC sums.';
metadata.lag_boundary='Each subsample constructs its own lags; calibration returns are not used as estimation-sample regressors.';
metadata.interpretation='Failure to reject Phi=0 does not establish independence, constant mean, constant variance, iid returns, or V=0. ARCH rejection indicates squared-residual dependence. All references are asymptotic.';
metadata.outputs=struct('main',fullfile(outdir,'return_ar1.csv'), ...
    'masked',fullfile(outdir,'return_ar1_nky_masked.csv'), ...
    'metadata',fullfile(outdir,'return_ar1_metadata.json'));
writetable(results,metadata.outputs.main);
if ~isempty(maskedRows), writetable(maskedResults,metadata.outputs.masked); end
fid=fopen(metadata.outputs.metadata,'w');
assert(fid>=0,'ar1:Write','Cannot write metadata.');
guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(metadata,PrettyPrint=true));
result=struct('results',results,'masked_results',maskedResults,'metadata',metadata);
disp(results(:,{'Sample','Ticker','NRegression','MeanWeeklyPct','Phi','PhiSE','PhiP','PhiPHolm','AR12P','ARCH12P'}));
end

function row=diagnostics(x)
x=x(:); T=numel(x);
meanFit=fit_hac(x,ones(T,1),8);
lag1=lag_matrix(x,1);
fit=fit_hac(x,[ones(T,1),lag1],8);
se=sqrt(diag(fit.Covariance));
z=sqrt(2)*erfcinv(.05);
row=struct('NCalendarWeeks',T,'NReturns',sum(isfinite(x)), ...
    'NRegression',fit.N,'MeanWeeklyPct',meanFit.Beta(1), ...
    'MeanSE',sqrt(meanFit.Covariance(1)), ...
    'MeanP',normal_p(meanFit.Beta(1)/sqrt(meanFit.Covariance(1))), ...
    'MeanPositiveP',.5*erfc(meanFit.Beta(1)/sqrt(2*meanFit.Covariance(1))), ...
    'Intercept',fit.Beta(1),'InterceptSE',se(1), ...
    'InterceptP',normal_p(fit.Beta(1)/se(1)), ...
    'Phi',fit.Beta(2),'PhiSE',se(2), ...
    'PhiCILow',fit.Beta(2)-z*se(2),'PhiCIHigh',fit.Beta(2)+z*se(2), ...
    'PhiP',normal_p(fit.Beta(2)/se(2)),'R2',fit.R2);
for h=[0,4,13]
    sensitivity=fit_hac(x,[ones(T,1),lag1],h);
    row.(sprintf('PhiP_HAC%d',h))=normal_p(sensitivity.Beta(2)/sqrt(sensitivity.Covariance(2,2)));
end
joint=fit_hac(x,[ones(T,1),lag_matrix(x,12)],8);
v=joint.Covariance(2:end,2:end); b=joint.Beta(2:end);
assert(rcond(v)>1e-12,'ar1:SingularWald','AR(12) covariance is numerically singular.');
row.AR12N=joint.N;
row.AR12Wald=b'*(v\b);
row.AR12P=gammainc(max(row.AR12Wald,0)/2,6,'upper');
squared=fit.Residuals.^2;
arch=fit_hac(squared,[ones(T,1),lag_matrix(squared,12)],0);
row.ARCH12N=arch.N;
row.ARCH12LM=arch.N*arch.R2;
row.ARCH12P=gammainc(max(row.ARCH12LM,0)/2,6,'upper');
end

function fit=fit_hac(y,X,maxLag)
% Retain the original calendar in the score sequence, including missing rows.
valid=isfinite(y)&all(isfinite(X),2);
N=sum(valid); K=size(X,2);
assert(N>K+2,'ar1:ShortRegression','Too few complete rows.');
A=X(valid,:); b=y(valid);
bread=A'*A;
assert(rcond(bread)>1e-12,'ar1:SingularDesign','Regression design is numerically singular.');
beta=A\b;
residuals=nan(size(y)); residuals(valid)=b-A*beta;
scores=zeros(size(X)); scores(valid,:)=A.*residuals(valid);
meat=scores'*scores;
for h=1:min(maxLag,size(X,1)-1)
    cross=scores(h+1:end,:)'*scores(1:end-h,:);
    meat=meat+(1-h/(maxLag+1))*(cross+cross');
end
covariance=(bread\meat)/bread*(N/(N-K));
covariance=(covariance+covariance')/2;
assert(all(isfinite(covariance(:)))&&all(diag(covariance)>0), ...
    'ar1:Covariance','Invalid HAC covariance.');
ssTotal=sum((b-mean(b)).^2);
assert(ssTotal>0,'ar1:ConstantSeries','The dependent variable is constant.');
fit=struct('Beta',beta,'Covariance',covariance,'N',N, ...
    'Residuals',residuals,'R2',max(0,1-sum(residuals(valid).^2)/ssTotal));
end

function X=lag_matrix(x,p)
X=nan(numel(x),p);
for lag=1:p, X(lag+1:end,lag)=x(1:end-lag); end
end

function p=normal_p(z)
p=erfc(abs(z)/sqrt(2));
end

function adjusted=holm(p)
adjusted=nan(size(p)); valid=find(isfinite(p));
[sorted,order]=sort(p(valid)); m=numel(sorted);
sortedAdjusted=min(1,cummax(sorted.*(m-(1:m)'+1)));
adjusted(valid(order))=sortedAdjusted;
end
