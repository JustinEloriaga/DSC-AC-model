function diagnostic = dsc_mcmc_diagnostics(draws)
%DSC_MCMC_DIAGNOSTICS Rank/split diagnostics for one scalar quantity.
% Input is retained iterations-by-ORIGINAL-chains, never pooled chains.
% R-hat is max(rank-normalized split, rank-normalized folded split).
% Bulk/tail ESS and raw-mean MCSE follow stan-dev/posterior R/convergence.R
% (inspected 2026-09-26), including Geyer's positive AND monotone sequence,
% biased FFT autocovariances, and the antithetic ESS cap N*log10(N):
% https://github.com/stan-dev/posterior/blob/master/R/convergence.R
% https://github.com/stan-dev/posterior/blob/master/R/split_chains.R
% Vehtari et al. (2021), Bayesian Analysis 16:667-718, doi:10.1214/20-BA1221.
% With an odd length the middle iteration is omitted ONLY when splitting.
% MCSE uses RAW-mean ESS, not bulk ESS. Single-chain values are exploratory.
% A constant original chain is conservatively unavailable; binary tail
% indicator subchains may be constant, as permitted by posterior's ESS.
diagnostic=struct('rhat',NaN,'ess_bulk',NaN,'ess_tail',NaN, ...
    'mcse_mean',NaN,'mcse_sd_ratio',NaN,'available',false, ...
    'reason','','original_chains',size(draws,2),'draws_per_chain',size(draws,1));
if ~isnumeric(draws)||~isreal(draws)||~ismatrix(draws)
    error('dsc:DiagnosticShape','Draws must be a real numeric iterations-by-chains matrix.');
end
if size(draws,1)<6||size(draws,2)<1
    diagnostic.reason='fewer_than_six_retained_draws'; return
end
draws=double(draws);
if any(~isfinite(draws(:)))
    diagnostic.reason='nonfinite_draws'; return
end
if any(max(draws,[],1)==min(draws,[],1))
    diagnostic.reason='constant_original_chain'; return
end
x=split_chains(draws);
z=rank_normalize(x);
folded=split_chains(abs(draws-median(draws(:))));
rhat_bulk=basic_rhat(z);
rhat_folded=basic_rhat(rank_normalize(folded));
% Do not let a finite component hide an unavailable folded diagnostic.
if isfinite(rhat_bulk)&&isfinite(rhat_folded)
    diagnostic.rhat=max(rhat_bulk,rhat_folded);
end
diagnostic.ess_bulk=effective_sample_size(z);
q=quantile_type7(draws(:),[.05 .95]);
tails=[effective_sample_size(split_chains(double(draws<=q(1)))), ...
    effective_sample_size(split_chains(double(draws<=q(2))))];
if all(isfinite(tails)), diagnostic.ess_tail=min(tails); end
ess_mean=effective_sample_size(x);
posterior_sd=std(draws(:),0);
if isfinite(ess_mean)&&ess_mean>0&&isfinite(posterior_sd)&&posterior_sd>0
    diagnostic.mcse_mean=posterior_sd/sqrt(ess_mean);
    diagnostic.mcse_sd_ratio=diagnostic.mcse_mean/posterior_sd;
end
metrics=[diagnostic.rhat diagnostic.ess_bulk diagnostic.ess_tail ...
    diagnostic.mcse_mean diagnostic.mcse_sd_ratio];
diagnostic.available=all(isfinite(metrics));
if ~diagnostic.available, diagnostic.reason='degenerate_or_unavailable_statistic'; end
end

function x=split_chains(x)
n=size(x,1); half=floor(n/2);
x=[x(1:half,:) x(n-half+1:n,:)];
end

function z=rank_normalize(x)
% Average pooled ranks for ties; Blom offset (rank-3/8)/(S+1/4).
[values,order]=sort(x(:)); count=numel(values);
first=[1;find(diff(values)~=0)+1]; last=[first(2:end)-1;count];
groups=cumsum([1;diff(values)~=0]);
average=(first+last)/2; ranks=zeros(count,1);
ranks(order)=average(groups);
z=reshape(-sqrt(2)*erfcinv(2*(ranks-3/8)/(count+1/4)),size(x));
end

function q=quantile_type7(x,probabilities)
% R's default quantile interpolation (MATLAB's default differs).
x=sort(x); position=1+(numel(x)-1)*probabilities;
lower=floor(position); upper=ceil(position); weight=position-lower;
q=(1-weight).*x(lower)'+weight.*x(upper)';
end

function r=basic_rhat(x)
n=size(x,1); within=mean(var(x,0,1)); between=n*var(mean(x,1),0);
r=NaN;
if isfinite(within)&&within>0&&isfinite(between)
    r=sqrt((between/within+n-1)/n);
end
end

function ess=effective_sample_size(x)
[n,m]=size(x); ess=NaN;
if n<3||any(~isfinite(x(:)))||max(x(:))==min(x(:)), return; end
centered=x-mean(x,1);
nfft=2^nextpow2(2*n);
transformed=fft(centered,nfft,1);
autocov=real(ifft(abs(transformed).^2,[],1))/n;
autocov_mean=mean(autocov(1:n,:),2);
within=autocov_mean(1)*n/(n-1);
var_plus=within*(n-1)/n+var(mean(x,1),0);
if ~isfinite(var_plus)||var_plus<=0, return; end
% Indices t below denote zero-based lags, matching posterior's algorithm.
rho=zeros(n,1); t=0; rho_even=1;
rho(1)=rho_even; rho_odd=1-(within-autocov_mean(2))/var_plus;
rho(2)=rho_odd;
while t<n-5&&isfinite(rho_even+rho_odd)&&(rho_even+rho_odd)>0
    t=t+2;
    rho_even=1-(within-autocov_mean(t+1))/var_plus;
    rho_odd=1-(within-autocov_mean(t+2))/var_plus;
    if rho_even+rho_odd>=0
        rho(t+1)=rho_even; rho(t+2)=rho_odd;
    end
end
max_t=t;
if rho_even>0, rho(max_t+1)=rho_even; end
% Initial monotone sequence of adjacent even/odd lag sums.
t=0;
while t<=max_t-4
    t=t+2;
    if rho(t+1)+rho(t+2)>rho(t-1)+rho(t)
        rho(t+1)=(rho(t-1)+rho(t))/2; rho(t+2)=rho(t+1);
    end
end
total=n*m;
% posterior uses R's inclusive 1:max_t indexing; when max_t is zero,
% 1:0 still selects rho[1]. Preserve that short-chain boundary explicitly.
if max_t==0, truncated_sum=rho(1); else, truncated_sum=sum(rho(1:max_t)); end
tau=-1+2*truncated_sum+rho(max_t+1);
tau=max(tau,1/log10(total));
if isfinite(tau)&&tau>0, ess=total/tau; end
end
