function bands = dsc_credible_bands(chain_draws, iterations, burnin)
%DSC_CREDIBLE_BANDS Pointwise 68% equal-tailed intervals for actual P entries.
% chain_draws{c}: dates x pairs x draws, containing P(i,j,t), NOT r or v.
% iterations{c}: original iteration identifiers for those draws.
% burnin(c): original-iteration warm-up cutoff for each chain (scalar allowed).
% Packed sampler chunks are already post-warm-up: keep their original
% iteration identifiers to avoid subtracting burn-in a second time.
% Call on date/pair blocks if necessary to bound memory. This routine does
% not run MCMC, assess convergence, or authorize publication of its output.
if ~iscell(chain_draws) || isempty(chain_draws) || ...
        ~iscell(iterations) || numel(iterations)~=numel(chain_draws)
    error('dsc:BandsInputs','Supply matching nonempty cells of chains and iteration IDs.');
end
nchains=numel(chain_draws);
if isscalar(burnin), burnin=repmat(burnin,1,nchains); end
if numel(burnin)~=nchains || any(~isfinite(burnin(:))) || ...
        any(burnin(:)<0 | burnin(:)~=floor(burnin(:)))
    error('dsc:BandsBurnin','Provide one nonnegative integer cutoff per chain.');
end
T=size(chain_draws{1},1); K=size(chain_draws{1},2);
if T==0 || K==0, error('dsc:BandsShape','Dates and pairs cannot be empty.'); end
keep=cell(1,nchains); retained=cell(1,nchains); counts=zeros(1,nchains);
for c=1:nchains
    x=chain_draws{c}; ids=iterations{c}(:)';
    if ~isnumeric(x) || ~isreal(x) || ndims(x)>3 || ...
            size(x,1)~=T || size(x,2)~=K || size(x,3)~=numel(ids)
        error('dsc:BandsShape','All chains must have matching dates/pairs and draw IDs.');
    end
    if any(~isfinite(ids)) || any(ids<1 | ids~=floor(ids)) || any(diff(ids)<=0)
        error('dsc:BandsIterations','Iteration IDs must be strictly increasing positive integers.');
    end
    keep{c}=find(ids>burnin(c)); retained{c}=ids(keep{c}); counts(c)=numel(keep{c});
    if counts(c)==0
        error('dsc:NoPosteriorDraws','Every supplied chain needs retained post-warm-up draws.');
    end
    if any(~isfinite(x(:))) || any(abs(x(:))>1)
        error('dsc:BandsCorrelations','Inputs must be finite actual correlations in [-1,1].');
    end
end
if sum(counts)<2, error('dsc:TooFewPosteriorDraws','At least two retained draws are required.'); end
bands.lower=nan(T,K); bands.median=nan(T,K); bands.upper=nan(T,K);
% Pool only one pair at a time; do not allocate a second full draw history.
for k=1:K
    samples=zeros(T,sum(counts)); offset=0;
    for c=1:nchains
        samples(:,offset+(1:counts(c)))=reshape(chain_draws{c}(:,k,keep{c}),T,counts(c));
        offset=offset+counts(c);
    end
    q=prctile(samples,[16 50 84],2);
    bands.lower(:,k)=q(:,1); bands.median(:,k)=q(:,2); bands.upper(:,k)=q(:,3);
end
bands.probability=0.68;
bands.percentiles=[16 50 84];
bands.interval_type='Pointwise equal-tailed; not simultaneous and not HPD';
bands.estimand='Actual conditional return-innovation correlations P(i,j,t)';
bands.retained_iterations=retained;
bands.draws_per_chain=counts;
bands.burnin=burnin;
bands.convergence_assessed=false;
bands.publication_ready=false;
bands.warning=['Quantiles alone do not establish convergence. Before publication require four chains, ' ...
    'rank-normalized split R-hat <1.01, bulk/tail ESS >=400 for reported quantities, ' ...
    'and model/numerical review. Intervals condition on empirical-Bayes calibrated priors.'];
end
