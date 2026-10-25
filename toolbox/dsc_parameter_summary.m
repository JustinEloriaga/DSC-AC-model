function summary = dsc_parameter_summary(draws,priors,tickers,pair_i,pair_j)
%DSC_PARAMETER_SUMMARY Marginal prior moments and retained posterior moments.
%   Rows follow mu, sig2h, sig2r, preserving ticker and strict lower-triangle
%   column-major pair order. PosteriorStd is the sample SD (n-1).
%   Prior moments are analytic; this function never draws random numbers.
tickers=string(tickers(:)); m=numel(tickers); nr=m*(m-1)/2;
if m<2||any(ismissing(tickers))||any(strlength(tickers)==0)|| ...
        numel(unique(tickers))~=m
    error('dsc:ParameterSummaryLabels','Provide distinct nonempty ticker labels.');
end
[expected_i,expected_j]=find(tril(true(m),-1));
if ~isequal(pair_i(:),expected_i)||~isequal(pair_j(:),expected_j)
    error('dsc:ParameterSummaryPairs','Pair indices must use column-major strict lower triangle order.');
end
if ~isstruct(draws)||~isscalar(draws)||~all(isfield(draws,{'mu','sig2h','sig2r'}))
    error('dsc:ParameterSummaryDraws','Provide retained mu, sig2h and sig2r draws.');
end
n=size(draws.mu,1);
for field={'mu','sig2h','sig2r'}
    name=field{1}; values=draws.(name); width=m;
    if strcmp(name,'sig2r'), width=nr; end
    if ~isnumeric(values)||~isreal(values)||n<2|| ...
            ~isequal(size(values),[n width])||any(~isfinite(values(:)))|| ...
            (~strcmp(name,'mu')&&any(values(:)<=0))
        error('dsc:ParameterSummaryDraws','At least two finite, dimensionally consistent retained draws are required; variances must be positive.');
    end
end
if ~isstruct(priors)||~isscalar(priors)
    error('dsc:ParameterSummaryPriors','Provide a scalar calibrated-prior structure.');
end
if isfield(priors,'mu0'), mu0=priors.mu0;
elseif isfield(priors,'Bbar'), mu0=priors.Bbar;
else, error('dsc:ParameterSummaryPriors','Missing normal prior mean mu0.'); end
if isfield(priors,'mu_cov0'), covariance=priors.mu_cov0;
elseif isfield(priors,'VBbar'), covariance=4*priors.VBbar;
else, error('dsc:ParameterSummaryPriors','Missing normal prior covariance mu_cov0.'); end
if ~isnumeric(mu0)||~isreal(mu0)||numel(mu0)~=m||any(~isfinite(mu0(:)))|| ...
        ~isnumeric(covariance)||~isreal(covariance)|| ...
        ~isequal(size(covariance),[m m])||any(~isfinite(covariance(:)))|| ...
        norm(covariance-covariance','fro')>1e-10*max(1,norm(covariance,'fro'))
    error('dsc:ParameterSummaryPriors','Normal prior dimensions and finite symmetric covariance must match the series.');
end
[~,bad]=chol((double(covariance)+double(covariance'))/2);
if bad, error('dsc:ParameterSummaryPriors','Normal prior covariance must be positive definite.'); end
if ~all(isfield(priors,{'ig_shape','ig_scale_h','ig_scale_r'}))
    error('dsc:ParameterSummaryPriors','Missing inverse-gamma shape or scale.');
end
[hmean,hstd]=inverse_gamma_moments(priors.ig_shape,priors.ig_scale_h,m);
[rmean,rstd]=inverse_gamma_moments(priors.ig_shape,priors.ig_scale_r,nr);
Parameter=[repmat("mu",m,1);repmat("sig2h",m,1);repmat("sig2r",nr,1)];
Series=[tickers;tickers;tickers(pair_i(:))+" / "+tickers(pair_j(:))];
PriorDistribution=[repmat("Normal",m,1);repmat("InverseGamma(shape,scale)",m+nr,1)];
PriorMean=[double(mu0(:));hmean;rmean];
PriorStd=[sqrt(diag(double(covariance)));hstd;rstd];
values=[double(draws.mu) double(draws.sig2h) double(draws.sig2r)];
PosteriorMean=mean(values,1)';
PosteriorStd=std(values,0,1)';
summary=table(Parameter,Series,PriorDistribution,PriorMean,PriorStd,PosteriorMean,PosteriorStd);
end

function [prior_mean,prior_std]=inverse_gamma_moments(shape,scale,width)
if ~isnumeric(shape)||~isreal(shape)||~isscalar(shape)||~isfinite(shape)||shape<=1|| ...
        ~isnumeric(scale)||~isreal(scale)||~isvector(scale)|| ...
        ~(isscalar(scale)||numel(scale)==width)||any(~isfinite(scale(:)))||any(scale(:)<=0)
    error('dsc:ParameterSummaryPriors','Inverse-gamma shape must be a finite scalar exceeding one; scales must be positive finite scalars or match their state coordinates.');
end
shape=double(shape); scale=double(scale(:));
% Historical saved priors used one shared scale. Broadcast only locally;
% coordinate-specific scales retain their saved series/pair ordering.
if isscalar(scale), scale=repmat(scale,width,1); end
prior_mean=scale/(shape-1);
if shape>2
    prior_std=scale/((shape-1)*sqrt(shape-2));
else
    prior_std=Inf(width,1); % The mean exists but the variance is infinite.
end
end
