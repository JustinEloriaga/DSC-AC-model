function priors = dsc_calibrate_priors(panel, cfg)
% Documented empirical Bayes: initial states use the first prior_weeks;
% Each evolution variance uses its own rolling-path changes over the complete
% historical panel; no pooling across volatility or correlation coordinates.
Y = panel.returns;
[T,m] = size(Y);
mask = logical(panel.observation_mask) & isfinite(Y);
N0 = min(cfg.prior_weeks,T);
X0 = Y(1:N0,:); X0 = X0(all(mask(1:N0,:),2),:);
if size(X0,1) <= m+1
    error('dsc:PriorRank','Initial window needs more than m+1 complete observations.');
end
S0 = shrink_cov(X0,cfg.shrinkage);
sd = sqrt(diag(S0)); C0 = S0./(sd*sd');
nr = m*(m-1)/2; sel = tril(true(m),-1);
[Q,L] = eig((C0+C0')/2,'vector');
A = Q*diag(log(L))*Q';
priors.mean_model = 'constant';
priors.mu0 = mean(X0,1)';
mu_std_multiplier=2;
if isfield(cfg,'mu_prior_std_multiplier'), mu_std_multiplier=cfg.mu_prior_std_multiplier; end
validateattributes(mu_std_multiplier,{'numeric'},{'real','scalar','finite','positive'});
priors.mu_prior_std_multiplier=mu_std_multiplier;
priors.mu_cov0 = mu_std_multiplier^2*S0/size(X0,1);
priors.mh0 = log(diag(S0));
priors.mr0 = A(sel);
priors.h0_scale = 10;
priors.r0_scale = 10;
windows = cfg.calibration_windows(:)';
cal = repmat(struct('window',0,'usable_windows',0,'sig2h',[], ...
    'sig2r',[]),numel(windows),1);
for widx=1:numel(windows)
    w=windows(widx); cal(widx).window=w;
    H=nan(T,m); R=nan(T,nr);
    for t=w:T
        Xi=Y(t-w+1:t,:);
        Xi=Xi(all(mask(t-w+1:t,:),2),:);
        if size(Xi,1)<=m+1, continue; end
        S=shrink_cov(Xi,cfg.shrinkage);
        si=sqrt(diag(S)); C=S./(si*si');
        [Q,L]=eig((C+C')/2,'vector');
        A=Q*diag(log(L))*Q';
        H(t,:)=log(diag(S))'; R(t,:)=A(sel)';
    end
    dH=diff(H); dR=diff(R);
    keep=all(isfinite(dH),2)&all(isfinite(dR),2);
    cal(widx).usable_windows=sum(all(isfinite(H),2));
    if sum(keep)>1
        cal(widx).sig2h=var(dH(keep,:),0,1);
        cal(widx).sig2r=var(dR(keep,:),0,1);
    end
end
chosen=find(windows==cfg.calibration_window,1);
if isempty(chosen)||numel(cal(chosen).sig2h)~=m||numel(cal(chosen).sig2r)~=nr|| ...
        any(~isfinite(cal(chosen).sig2h))||any(~isfinite(cal(chosen).sig2r))
    error('dsc:CalibrationWindow','Selected calibration window has insufficient observations.');
end
innovation_mean_multiplier=2;
if isfield(cfg,'innovation_prior_mean_multiplier')
    innovation_mean_multiplier=cfg.innovation_prior_mean_multiplier;
end
validateattributes(innovation_mean_multiplier,{'numeric'}, ...
    {'real','scalar','finite','positive'});
priors.sig2h_mean=innovation_mean_multiplier*max(cal(chosen).sig2h,1e-8);
priors.sig2r_mean=innovation_mean_multiplier*max(cal(chosen).sig2r,1e-8);
priors.innovation_prior_mean_multiplier=innovation_mean_multiplier;
priors.innovation_prior_pooling='none';
% User-facing innovation-variance names. No dynamic mean covariance exists.
priors.delta_mean=priors.sig2h_mean;
priors.gamma_mean=priors.sig2r_mean;
priors.delta_definition='delta_j = Var(h_jt-h_j,t-1) = sig2h_j';
priors.gamma_definition='gamma_k = Var(r_kt-r_k,t-1) = sig2r_k';
priors.ig_shape=cfg.ig_shape;
if priors.ig_shape<=1, error('dsc:PriorShape','Inverse gamma shape must exceed one.'); end
% For each coordinate, twice its raw calibrated mean and a=4 give four times
% the SD of an a=10 prior centered on that same raw coordinate variance.
priors.ig_scale_h=priors.sig2h_mean*(priors.ig_shape-1);
priors.ig_scale_r=priors.sig2r_mean*(priors.ig_shape-1);
priors.calibration=cal;
priors.initial_weeks=N0;
priors.initial_complete_observations=size(X0,1);
priors.initial_dates=panel.dates([1 N0]);
priors.calibration_window=cfg.calibration_window;
priors.shrinkage=cfg.shrinkage;
priors.source_hash=panel.source_hash;
priors.tickers=panel.tickers;
priors.scope='Initial states: first configured weeks; separate evolution priors: own full-sample rolling-path change variances, without pooling.';
end

function S=shrink_cov(X,a)
if ~(isscalar(a)&&a>0&&a<=1), error('dsc:Shrinkage','Shrinkage must lie in (0,1].'); end
S=cov(X,1); S=(1-a)*S+a*diag(diag(S)); S=(S+S')/2;
if any(diag(S)<=0)||any(~isfinite(S(:)))
    error('dsc:PriorVariance','All variables need positive finite variance.');
end
[~,bad]=chol(S);
if bad, error('dsc:PriorCovariance','Shrunk covariance is not positive definite.'); end
end
