function summary = analyze_weekly_panel(panel,cfg)
%ANALYZE_WEEKLY_PANEL Export reproducible descriptive and stationarity results.
%   Outputs are written beneath cfg.output_root/data. Raw input files are
%   never edited. Rolling estimates use trailing 52/104-week windows only;
%   an incomplete pair-window is explicitly NaN, never silently shortened.
%   Stationarity tests concern the mean/unit-root specification. They do not
%   establish constant variance or constancy of the TVP model's parameters.

arguments
    panel (1,1) struct
    cfg (1,1) struct
end
if ~isfield(cfg,'output_root') || strlength(string(cfg.output_root)) == 0
    error('weekly:OutputRoot','Set cfg.output_root for reproducible analysis outputs.');
end
outdir = fullfile(cfg.output_root,'data');
if ~isfolder(outdir), mkdir(outdir); end
[T,m] = size(panel.returns);
q = numel(panel.pair_i);
assert(size(panel.levels,1)==T+1 && numel(panel.dates)==T, ...
    'weekly:Dimensions','Return and level date dimensions disagree.');
assert(all(diff(panel.dates)==days(7)), ...
    'weekly:Calendar','The weekly calendar must not have omitted periods.');
if ~isfield(panel,'observation_mask'), panel.observation_mask=true(T,m); end
save(fullfile(outdir,'weekly_panel.mat'),'panel','-v7.3');

dates = string(panel.dates,'yyyy-MM-dd');
levelDates = string(panel.level_dates,'yyyy-MM-dd');
returnTable = [table(dates,'VariableNames',{'Date'}), ...
    array2table(panel.returns,'VariableNames',cellstr(panel.tickers))];
levelTable = [table(levelDates,'VariableNames',{'Date'}), ...
    array2table(panel.levels,'VariableNames',cellstr(panel.tickers))];
writetable(returnTable,fullfile(outdir,'weekly_returns.csv'));
writetable(levelTable,fullfile(outdir,'weekly_levels.csv'));
nL = T+1;
quality = table(repmat(levelDates,m,1),repelem(panel.tickers',nL), ...
    reshape(string(panel.source_dates,'yyyy-MM-dd'),[],1), ...
    panel.age_days(:),panel.carried(:),panel.closure_level_mask(:), ...
    'VariableNames',{'LevelDate','Ticker','SourceDate','AgeDays','Carried','ClosureLevel'});
writetable(quality,fullfile(outdir,'weekly_quality.csv'));
returnQuality = table(repmat(dates,m,1),repelem(panel.tickers',T), ...
    panel.closure_return_mask(:),panel.observation_mask(:), ...
    'VariableNames',{'Date','Ticker','ClosureReturn','ObservedForModel'});
writetable(returnQuality,fullfile(outdir,'weekly_return_quality.csv'));

X = panel.returns;
X(~panel.observation_mask) = NaN;
seriesRows = cell(m,16);
for j=1:m
    x=X(isfinite(X(:,j)),j);
    raw=panel.raw_levels(:,j);
    seriesRows(j,:)={panel.tickers(j),panel.labels(j),numel(x), ...
        sum(isnan(raw)),sum(isnan(raw))/numel(raw), ...
        sum(panel.carried(:,j)),max(panel.age_days(:,j)), ...
        mean(x),std(x),52*mean(x),sqrt(52)*std(x),skewness(x),kurtosis(x), ...
        min(x),max(x),sum(x==0)};
end
series = cell2table(seriesRows,'VariableNames',{'Ticker','Label','N', ...
    'RawMissingCells','RawMissingFraction','CarriedEndpoints','MaxAgeDays', ...
    'MeanWeeklyPct','StdWeeklyPct','AnnualizedMeanLogPct','AnnualizedVolPct', ...
    'Skewness','Kurtosis','MinWeeklyPct','MaxWeeklyPct','ZeroReturns'});
writetable(series,fullfile(outdir,'series_summary.csv'));

% Export annual raw coverage without treating the partial final year as full.
yrs=unique(year(panel.raw_dates));
annualRows=cell(numel(yrs)*m,5); row=0;
for k=1:numel(yrs)
    ix=year(panel.raw_dates)==yrs(k);
    for j=1:m
        row=row+1; nMissing=sum(isnan(panel.raw_levels(ix,j)));
        annualRows(row,:)={yrs(k),panel.tickers(j),sum(ix),nMissing,nMissing/sum(ix)};
    end
end
annual=cell2table(annualRows,'VariableNames', ...
    {'Year','Ticker','RawRows','MissingCells','MissingFraction'});
writetable(annual,fullfile(outdir,'raw_annual_coverage.csv'));

pairIndex=(1:q)'; tickerI=panel.tickers(panel.pair_i)';
tickerJ=panel.tickers(panel.pair_j)';
fullCorr=nan(q,1); counts=zeros(q,1);
for k=1:q
    a=X(:,panel.pair_i(k)); b=X(:,panel.pair_j(k));
    good=isfinite(a)&isfinite(b); counts(k)=sum(good);
    if counts(k)>=3, c=corrcoef(a(good),b(good)); fullCorr(k)=c(1,2); end
end
correlations=table(pairIndex,tickerI,tickerJ,fullCorr,counts, ...
    'VariableNames',{'PairIndex','TickerI','TickerJ','Correlation','N'});
writetable(correlations,fullfile(outdir,'full_sample_correlations.csv'));

windows=[52,104];
rollingRows=q*sum(max(T-windows+1,0));
rDate=strings(rollingRows,1); rWindow=zeros(rollingRows,1);
rPair=zeros(rollingRows,1); rCorr=nan(rollingRows,1); rN=zeros(rollingRows,1);
row=0;
for w=windows
    for t=w:T
        x=X(t-w+1:t,:);
        complete=all(isfinite(x),1);
        C=nan(m,m);
        if sum(complete)>=2, C(complete,complete)=corrcoef(x(:,complete)); end
        ix=row+(1:q);
        rDate(ix)=dates(t); rWindow(ix)=w; rPair(ix)=pairIndex;
        rCorr(ix)=C(sub2ind([m,m],panel.pair_i,panel.pair_j));
        rN(ix)=sum(isfinite(x(:,panel.pair_i)) & isfinite(x(:,panel.pair_j)),1)';
        row=row+q;
    end
end
rolling=table(rDate,rWindow,rPair,tickerI(rPair),tickerJ(rPair),rCorr,rN, ...
    rN==rWindow,'VariableNames', ...
    {'Date','Window','PairIndex','TickerI','TickerJ','Correlation','N','CompleteWindow'});
writetable(rolling,fullfile(outdir,'rolling_correlations.csv'));

changeRows=cell(numel(windows)*q,15); row=0;
for w=windows
    for k=1:q
        ix=find(rWindow==w & rPair==k & isfinite(rCorr));
        row=row+1;
        if isempty(ix)
            changeRows(row,:)={w,k,tickerI(k),tickerJ(k),NaN,NaN,NaN, ...
                NaN,"",NaN,"",NaN,NaN,"",false};
            continue
        end
        rr=rCorr(ix); dd=rDate(ix); [lo,imin]=min(rr); [hi,imax]=max(rr);
        % A change must compare consecutive calendar weeks, never bridge gaps.
        dr=diff(rr); adjacent=diff(datetime(dd,'InputFormat','yyyy-MM-dd'))==days(7);
        dr(~adjacent)=NaN;
        if any(isfinite(dr))
            [big,iBig]=max(abs(dr),[],'omitnan'); signedStep=dr(iBig); stepDate=dd(iBig+1);
        else
            big=NaN; signedStep=NaN; stepDate="";
        end
        changeRows(row,:)={w,k,tickerI(k),tickerJ(k),rr(1),rr(end),rr(end)-rr(1), ...
            lo,dd(imin),hi,dd(imax),big,signedStep,stepDate,lo<0 && hi>0};
    end
end
changes=cell2table(changeRows,'VariableNames',{'Window','PairIndex','TickerI', ...
    'TickerJ','FirstCorrelation','LastCorrelation','NetChange','MinCorrelation', ...
    'MinDate','MaxCorrelation','MaxDate','MaxAbsStepChange','SignedLargestStep', ...
    'LargestStepDate','CrossesZero'});
writetable(changes,fullfile(outdir,'correlation_changes.csv'));

% Diagnostics use complete, regular marked-return data. Masks are reported,
% not removed: deleting missing observations would compress the time axis.
% These descriptive tests therefore refer to the as-of panel even when the
% optional model observation mask is active.
[stationarity,diagnostics,adfCandidates]=stationarityTables(panel);
writetable(stationarity,fullfile(outdir,'stationarity.csv'));
writetable(diagnostics,fullfile(outdir,'diagnostics.csv'));
writetable(adfCandidates,fullfile(outdir,'adf_candidates.csv'));

summary=struct();
summary.schema_version=1;
summary.n_returns=T; summary.n_levels=T+1; summary.n_series=m; summary.n_pairs=q;
summary.first_return=char(dates(1)); summary.last_return=char(dates(end));
summary.first_level=char(levelDates(1)); summary.last_level=char(levelDates(end));
summary.source_file=panel.source_file; summary.source_hash=panel.source_hash;
summary.tickers=cellstr(panel.tickers);
summary.currency_basis=panel.currency_basis; summary.frequency=panel.frequency;
summary.raw_rows=size(panel.raw_levels,1);
summary.raw_missing_cells=sum(isnan(panel.raw_levels(:)));
summary.raw_missing_fraction=mean(isnan(panel.raw_levels(:)));
summary.carried_endpoint_cells=sum(panel.carried(:));
summary.carried_endpoint_fraction=mean(panel.carried(:));
summary.fridays_with_carry=sum(any(panel.carried,2));
summary.max_age_days=max(panel.age_days(:));
summary.closure_level_cells=sum(panel.closure_level_mask(:));
summary.closure_return_cells=sum(panel.closure_return_mask(:));
summary.model_masked_return_cells=sum(~panel.observation_mask(:));
summary.closure_policy=panel.closure_policy;
summary.closure_reference=panel.closure_reference;
summary.excluded_final_raw_rows=panel.excluded_final_rows;
summary.stationarity_data='complete as-of marks, including flagged closure';
summary.adf_note='BIC lags 0:13 on common dependent sample; MATLAB p-values bounded [0.001,0.999]';
summary.kpss_note='Newey-West bandwidths 4,13,26; MATLAB p-values bounded [0.01,0.10]';
summary.mean_stationarity_note='Mean/unit-root diagnostics do not establish constant variance or invariant dynamic coefficients.';
summary.rolling_note='Trailing complete pair-windows only; no centered windows or future observations.';
summary.ages=struct('days',unique(panel.age_days(:))', ...
    'counts',arrayfun(@(a)sum(panel.age_days(:)==a),unique(panel.age_days(:))'));
ir=stationarity.Transform=="weekly_log_return_pct";
summary.adf_return_rejections=sum(stationarity.Reject5pct(ir & stationarity.Test=="ADF"));
summary.adf_invalid_selected_candidates=sum(adfCandidates.Selected & ~adfCandidates.ValidStatistic);
summary.kpss_return_rejections=struct();
for lag=[4,13,26]
    summary.kpss_return_rejections.(sprintf('lags%d',lag))= ...
        sum(stationarity.Reject5pct(ir & stationarity.Test=="KPSS" & stationarity.Lags==lag));
end
summary.ljung_box_return_rejections=sum(diagnostics.Reject5pct(diagnostics.Test=="LjungBoxReturns"));
summary.ljung_box_squared_return_rejections=sum(diagnostics.Reject5pct(diagnostics.Test=="LjungBoxSquaredReturns"));
summary.arch_rejections=sum(diagnostics.Reject5pct(diagnostics.Test=="ARCH"));
[~,order]=sort(fullCorr,'descend');
summary.top_positive_pairs=table2struct(correlations(order(1:min(5,q)),:));
[~,order]=sort(fullCorr,'ascend');
summary.top_negative_pairs=table2struct(correlations(order(1:min(5,q)),:));
ix=find(changes.Window==104 & isfinite(changes.NetChange));
[~,order]=sort(abs(changes.NetChange(ix)),'descend');
summary.largest_104_week_net_changes=table2struct(changes(ix(order(1:min(5,numel(order)))),:));
summary.series=table2struct(series);
fid=fopen(fullfile(outdir,'data_summary.json'),'w');
if fid<0, error('weekly:WriteSummary','Cannot write data_summary.json.'); end
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(summary,'PrettyPrint',true));
clear cleanup
save(fullfile(outdir,'panel_analysis.mat'),'summary','series','correlations', ...
    'changes','stationarity','diagnostics','adfCandidates','-v7.3');
end

function [stationarity,diagnostics,adfCandidates]=stationarityTables(panel)
m=numel(panel.tickers);
rows=cell(15*m,11); row=0;
diagnosticRows=cell(3*m,8); drow=0;
candidateRows=cell(3*m*14,12); crow=0;
for j=1:m
    % Both plausible deterministic specifications are exposed for levels.
    ys={log(panel.levels(:,j)),log(panel.levels(:,j)),panel.returns(:,j)};
    transforms=["log_level","log_level","weekly_log_return_pct"];
    models=["ARD","TS","ARD"];
    for v=1:3
        y=ys{v}; maxLag=13;
        if numel(y)<=2*maxLag+5
            error('weekly:DiagnosticSample','Diagnostics require more than %d observations.',2*maxLag+5);
        end
        bic=nan(maxLag+1,1); hh=false(maxLag+1,1);
        pp=nan(maxLag+1,1); ss=nan(maxLag+1,1); nn=zeros(maxLag+1,1);
        valid=true(maxLag+1,1); warningId=strings(maxLag+1,1);
        for lag=0:maxLag
            common=y(maxLag-lag+1:end);
            lastwarn('');
            [hh(lag+1),pp(lag+1),ss(lag+1),~,reg]= ...
                adftest(common,'Model',models(v),'Lags',lag,'Alpha',0.05);
            [~,warningId(lag+1)]=lastwarn;
            valid(lag+1)=warningId(lag+1)~="econ:adftest:InvalidStatistic";
            bic(lag+1)=reg.BIC; nn(lag+1)=reg.size;
        end
        assert(all(nn==nn(1)),'weekly:ADFCommonSample','ADF lag samples differ.');
        [bestBIC,chosen]=min(bic); lag=chosen-1;
        if v<=2, testDates=panel.level_dates; else, testDates=panel.dates; end
        for candidate=0:maxLag
            crow=crow+1; k=candidate+1;
            candidateRows(crow,:)={panel.tickers(j),transforms(v),models(v),candidate, ...
                bic(k),nn(k),k==chosen,valid(k),warningId(k),ss(k), ...
                string(testDates(maxLag+2),'yyyy-MM-dd'),string(testDates(end),'yyyy-MM-dd')};
        end
        if ~valid(chosen)
            error('weekly:InvalidADF','BIC-selected ADF statistic is invalid for %s/%s/%s.', ...
                panel.tickers(j),transforms(v),models(v));
        end
        row=row+1;
        rows(row,:)={panel.tickers(j),transforms(v),"ADF",models(v),lag,ss(chosen), ...
            pp(chosen),boundLabel(pp(chosen),0.001,0.999),hh(chosen),nn(chosen),bestBIC};
        trend=models(v)=="TS";
        for bandwidth=[4,13,26]
            [h,p,stat]=kpsstest(y,'Trend',trend,'Lags',bandwidth,'Alpha',0.05);
            if trend, name="trend"; else, name="level"; end
            row=row+1;
            rows(row,:)={panel.tickers(j),transforms(v),"KPSS",name,bandwidth, ...
                stat,p,boundLabel(p,0.01,0.10),h,numel(y),NaN};
        end
    end
    y=panel.returns(:,j); residual=y-mean(y);
    [h1,p1,s1]=lbqtest(y,'Lags',13,'Alpha',0.05);
    [h2,p2,s2]=lbqtest(residual.^2,'Lags',13,'Alpha',0.05);
    [h3,p3,s3]=archtest(residual,'Lags',13,'Alpha',0.05);
    names=["LjungBoxReturns","LjungBoxSquaredReturns","ARCH"];
    hs=[h1,h2,h3]; stats=[s1,s2,s3];
    % Evaluate the chi-square survival function directly: 1-CDF loses tiny
    % positive probabilities to cancellation in MATLAB's test wrappers.
    ps=chi2cdf(stats,13,'upper');
    for k=1:3
        drow=drow+1;
        % Numeric zero is explicitly labeled as numerical underflow, never an
        % exact probability of zero. Publication should print '< 1e-12'.
        if ps(k)<1e-12, label="<1e-12"; else, label="numeric"; end
        diagnosticRows(drow,:)={panel.tickers(j),names(k),13,stats(k),ps(k),label,hs(k),numel(y)};
    end
end
stationarity=cell2table(rows(1:row,:),'VariableNames',{'Ticker','Transform','Test', ...
    'Model','Lags','Statistic','PValue','PValueBound','Reject5pct','N','BIC'});
diagnostics=cell2table(diagnosticRows,'VariableNames',{'Ticker','Test','Lags', ...
    'Statistic','PValue','PValueBound','Reject5pct','N'});
adfCandidates=cell2table(candidateRows,'VariableNames',{'Ticker','Transform','Model', ...
    'Lags','BIC','N','Selected','ValidStatistic','WarningID','Statistic', ...
    'FirstDependentDate','LastDependentDate'});
end

function label=boundLabel(p,lower,upper)
if p<=lower
    label="<="+string(lower);
elseif p>=upper
    label=">="+string(upper);
else
    label="interpolated";
end
end
