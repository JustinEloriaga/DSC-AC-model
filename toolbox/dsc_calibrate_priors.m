function priors = dsc_calibrate_priors(panel, cfg)
% Documented empirical Bayes: initial states use the first prior_weeks;
% evolution scales use rolling windows over the complete historical panel.
Y = panel.returns;
[T,m] = size(Y);
if cfg.p ~= 0, error('dsc:UnsupportedVAR','This sampler supports p=0.'); end
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
priors.Bbar = mean(X0,1)';
priors.VBbar = S0/size(X0,1);
priors.nuB = size(X0,1);
priors.V0B = priors.VBbar*priors.nuB*cfg.kB^2;
priors.mh0 = log(diag(S0));
priors.mr0 = A(sel);
priors.h0_scale = 10;
priors.r0_scale = 10;
windows = cfg.calibration_windows(:)';
cal = repmat(struct('window',0,'usable_windows',0,'sig2h',[], ...
    'sig2r',[],'sig2h_median',NaN,'sig2r_median',NaN),numel(windows),1);
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
        cal(widx).sig2h_median=median(cal(widx).sig2h);
        cal(widx).sig2r_median=median(cal(widx).sig2r);
    end
end
chosen=find(windows==cfg.calibration_window,1);
if isempty(chosen)||~isfinite(cal(chosen).sig2h_median)||~isfinite(cal(chosen).sig2r_median)
    error('dsc:CalibrationWindow','Selected calibration window has insufficient observations.');
end
priors.sig2h_mean=max(cal(chosen).sig2h_median,1e-8);
priors.sig2r_mean=max(cal(chosen).sig2r_median,1e-8);
priors.ig_shape=cfg.ig_shape;
if priors.ig_shape<=1, error('dsc:PriorShape','Inverse gamma shape must exceed one.'); end
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
priors.scope='Initial states: first configured weeks; evolution scales: full-sample rolling empirical Bayes.';
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
