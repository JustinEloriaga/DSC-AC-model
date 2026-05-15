function r = tvsvar_modified_msv2_gam2_gen(y, lags, T0, kB, M, nreport, nthins)
% TVP-VAR with Hansen's MSV2 parameterization of the correlation matrix.
% Generalized to any number of variables.

if nargin ~= 7; error('Wrong # of arguments'); end

if lags < 0 || lags ~= floor(lags)
    error('lags must be a nonnegative integer');
end

% Build regressor matrix including lags and constant
nobs = size(y, 1);
x = ones(nobs, 1);
for i = 1:lags
    x = [x, lag(y,i)];
end

x  = x(lags+1:end, :);
y  = y(lags+1:end, :);
[T, n] = size(y);

x0 = x(1:T0, :);
x  = x(T0+1:end, :);
y0 = y(1:T0, :);
y  = y(T0+1:end, :);
Y  = reshape(y, n*(T-T0), 1);
Z  = kron(eye(n), x);

T  = T - T0;
k  = n + lags*n^2;

% Storage for state and covariance
B_old = zeros(T, k);
V_old = zeros(k, k);

% Prior for B from OLS on training sample
priorB = olsblock(y0, x0);
Bbar   = reshape(priorB.bhatols, k, 1)';
resB   = priorB.resols;
VBbar  = kron(resB'*resB/length(resB), eye(k/n)/priorB.XX);

% Prior variance for B evolution
iwishmean = @(V, T) (V / (T + size(V,1) + 1));
V0B = VBbar * T0 * kB^2;
V_old(1:k, 1:k) = iwishmean(V0B, T0);

% --- MSV setup ---
m = n;

% Load empirical calibration of sig2h, sig2r prior means if present;
% otherwise fall back to literature defaults.
calib_file = which('calibrated_priors.mat');
if ~isempty(calib_file)
    calib = load(calib_file, 'sig2h_prior_mean', 'sig2r_prior_mean');
    sig2h_pm = calib.sig2h_prior_mean;
    sig2r_pm = calib.sig2r_prior_mean;
    fprintf('  Calibrated priors loaded: sig2h_pm=%.4g, sig2r_pm=%.4g\n', sig2h_pm, sig2r_pm);
else
    sig2h_pm = 0.01;
    sig2r_pm = 0.005;
    fprintf('  No calibration file; default priors: sig2h_pm=%.4g, sig2r_pm=%.4g\n', sig2h_pm, sig2r_pm);
end
T0_pseudo = 10;     % prior strength: ~10 pseudo-observations

% Prior for sig2h: inverse-gamma with mean V0/(T0-1) = sig2h_pm
prior_sig2h.T0 = ones(m,1) * T0_pseudo;
prior_sig2h.V0 = ones(m,1) * sig2h_pm * (T0_pseudo - 1);

% Prior for h(0)
mh0       = zeros(m, 1);
Vh0_scale = 10 * ones(m, 1);

% Slice sampler covariance for H(t)
slice.nobs       = T;
slice_chol_cov_z = zeros(T, T, m);
slice_mean_z     = zeros(T, m);
for vvind = 1:m
    var_z = Vh0_scale(vvind) + (1:T)';
    cov_z = tril(repmat(var_z', T, 1), 0) + tril(repmat(var_z', T, 1), -1)';
    slice_chol_cov_z(:,:,vvind) = chol(cov_z, 'lower');
    slice_mean_z(:,vvind)       = mh0(vvind) * ones(T, 1);
end

% Slice sampler covariance for R(t)
n_r       = m*(m-1)/2;
Vr0_scale = 10 * ones(n_r, 1);
mr0       = zeros(n_r, 1);

var_z          = Vr0_scale(1) + (1:T)';
cov_z          = tril(repmat(var_z', T, 1), 0) + tril(repmat(var_z', T, 1), -1)';
slice_chol_cov_z_r = chol(cov_z, 'lower');
slice_mean_z_r     = mr0(1) * ones(T, 1);

for vvind = 2:n_r
    var_z = Vr0_scale(vvind) + (1:T)';
    cov_z = tril(repmat(var_z', T, 1), 0) + tril(repmat(var_z', T, 1), -1)';
    slice_chol_cov_z_r = cat(3, slice_chol_cov_z_r, chol(cov_z, 'lower'));
    slice_mean_z_r     = cat(2, slice_mean_z_r, mr0(vvind)*ones(T,1));
end

% Prior for sig2r: inverse-gamma with mean V0/(T0-1) = sig2r_pm
prior_sig2r.T0 = ones(n_r, 1) * T0_pseudo;
prior_sig2r.V0 = ones(n_r, 1) * sig2r_pm * (T0_pseudo - 1);

%% Initialization
yhat      = olsblock(y, x).resols;

% Initialize sig2h, sig2r at their prior means (calibrated from data).
% Avoids starting with the trivial var(yhat)/T values which were ~10-100x
% too small and forced the slice sampler into a long warm-up drift.
sig2h_old = sig2h_pm * ones(1, m);

ht_old = zeros(T, m);
for t = 1:T
    ht_old(t,:) = log(var(yhat)) + sqrt(sig2h_old) .* randn(1, m);
end

[temp_v1, temp_v2] = eig(corr(yhat), 'vector');
temp_v3    = temp_v1 * diag(log(temp_v2)) * temp_v1';
temp_r     = temp_v3(tril(ones(m), -1) == 1);
sig2r_old  = sig2r_pm * ones(1, n_r);

rt_old = zeros(T, n_r);
for t = 1:T
    rt_old(t,:) = temp_r' + sqrt(sig2r_old) .* randn(1, n_r);
end

Pt_old = zeros(m, m, T);
for t = 1:T
    Pt_old(:,:,t) = veclAtoC(rt_old(t,:));
end

%% Storage
if nthins < 1 || nthins ~= floor(nthins)
    error('nthins must be a positive integer');
end
nsave = floor(M / nthins);
if nsave < 1
    error('M and nthins imply zero saved draws; increase M or reduce nthins');
end

mat_ht   = zeros(T, m, nsave);
mat_sig2h = zeros(nsave, m);
mat_rt   = zeros(T, n_r, nsave);
mat_sig2r = zeros(nsave, n_r);
mat_P    = zeros(m, m, T, nsave);
mat_B    = zeros(1, k, nsave);
mat_V    = zeros(k, k, nsave);

%% Gibbs sampler
time_all    = tic;
time_report = tic;
saveind     = 1;

for sind = 1:M

    % --- Step 1: Draw B (Kalman filter + simulation smoother) ---
    SHAT = zeros(T, k);
    SIG  = zeros(T, k, k);
    sig  = 4 * VBbar;
    shat = Bbar';
    for t = 1:T
        OMEGA_t    = diag(exp(ht_old(t,:)/2)) * Pt_old(:,:,t) * diag(exp(ht_old(t,:)/2));
        [shat, sig] = kfilter(y(t,:)', Z([0:n-1]*T+t,:), eye(k), shat, sig, OMEGA_t, squeeze(V_old(1:k,1:k)));
        SHAT(t,:)   = shat';
        SIG(t,:,:)  = sig;
    end

    B_old(T,:) = mvnrnd_modified(shat, sig/2+sig'/2, 1);
    for t = T-1:-1:1
        [btTp, StTp] = kback(SHAT(t,:), squeeze(SIG(t,:,:)), squeeze(B_old(t+1,:)), eye(k), squeeze(V_old(1:k,1:k)));
        B_old(t,:)   = mvnrnd_modified(btTp, StTp/2+StTp'/2, 1);
    end

    % --- Step 2: Draw V ---
    errB = B_old(2:end,:) - B_old(1:end-1,:);
    V1B  = errB'*errB + V0B;
    V_old(1:k,1:k) = iwishrnd((V1B + V1B')/2, T-1+T0);

    % --- Step 3: Draw R(t), H(t) and variances ---
    Yhat   = Y - sum((Z .* repmat(squeeze(B_old(:,:)), n, 1))')';
    et_old = reshape(Yhat, T, n);

    % Dynamic correlation
    lik_old = loglike_yt_given_rt_ht(et_old, ht_old, rt_old, n_r, rt_old(:,n_r));
    for vvind = 1:n_r
        slice.scale_z    = sqrt(sig2r_old(vvind));
        slice.mean       = slice_mean_z_r(:,vvind);
        slice.chol_cov_z = slice_chol_cov_z_r(:,:,vvind);
        slice.fcn_lik    = @(rt_vvind) loglike_yt_given_rt_ht(et_old, ht_old, rt_old, vvind, rt_vvind);
        [rt_old(:,vvind), lik_old] = slice_sampling_v02(slice, rt_old(:,vvind), lik_old);
    end
    for t = 1:T
        Pt_old(:,:,t) = veclAtoC(rt_old(t,:));
    end

    % Draw sig2r
    errR = rt_old(2:end,:) - rt_old(1:end-1,:);
    for vvind = 1:n_r
        EE              = errR(:,vvind);
        T1              = size(errR,1)/2 + prior_sig2r.T0(vvind);
        V1              = EE'*EE/2 + prior_sig2r.V0(vvind);
        sig2r_new       = gamrnd(T1, V1^(-1))^(-1);
        logalp = lognormpdf(rt_old(1,vvind), mr0(vvind), sqrt((Vr0_scale(vvind)+1)*sig2r_new)) ...
               - lognormpdf(rt_old(1,vvind), mr0(vvind), sqrt((Vr0_scale(vvind)+1)*sig2r_old(vvind)));
        if log(rand) < logalp
            sig2r_old(vvind) = sig2r_new;
        end
    end

    % Dynamic variance
    lik_old = loglike_yt_given_ht_Pt(et_old, ht_old, Pt_old, m, ht_old(:,m));
    for vvind = 1:m
        slice.scale_z    = sqrt(sig2h_old(vvind));
        slice.mean       = slice_mean_z(:,vvind);
        slice.chol_cov_z = slice_chol_cov_z(:,:,vvind);
        slice.fcn_lik    = @(ht_vvind) loglike_yt_given_ht_Pt(et_old, ht_old, Pt_old, vvind, ht_vvind);
        [ht_old(:,vvind), lik_old] = slice_sampling_v02(slice, ht_old(:,vvind), lik_old);
    end

    % Draw sig2h
    errH = ht_old(2:end,:) - ht_old(1:end-1,:);
    for vvind = 1:m
        EE        = errH(:,vvind);
        T1        = size(errH,1)/2 + prior_sig2h.T0(vvind);
        V1        = EE'*EE/2 + prior_sig2h.V0(vvind);
        sig2h_new = gamrnd(T1, V1^(-1))^(-1);
        logalp = lognormpdf(ht_old(1,vvind), mh0(vvind), sqrt((Vh0_scale(vvind)+1)*sig2h_new)) ...
               - lognormpdf(ht_old(1,vvind), mh0(vvind), sqrt((Vh0_scale(vvind)+1)*sig2h_old(vvind)));
        if log(rand) < logalp
            sig2h_old(vvind) = sig2h_new;
        end
    end

    % --- Progress report ---
    if rem(sind, nreport) == 0
        elapsed_total  = toc(time_all);
        elapsed_report = toc(time_report);
        est_left       = elapsed_report * (M - sind) / nreport;
        disp(['sind = ', num2str(sind), ' / ', num2str(M), ...
              '  |  block: ', num2str(elapsed_report,'%.1f'), 's', ...
              '  (total: ', num2str(elapsed_total,'%.1f'), 's)', ...
              '  |  remaining: ~', num2str(est_left/60,'%.1f'), ' min']);
        time_report = tic;
    end

    % --- Save draw ---
    if rem(sind, nthins) == 0
        mat_B(:,:,saveind)   = B_old(end,:);
        mat_V(:,:,saveind)   = V_old;
        mat_rt(:,:,saveind)  = rt_old;
        mat_P(:,:,:,saveind) = Pt_old;
        mat_sig2r(saveind,:) = sig2r_old;
        mat_ht(:,:,saveind)  = ht_old;
        mat_sig2h(saveind,:) = sig2h_old;
        saveind = saveind + 1;
    end
end

%% Output
r.B     = mat_B;
r.V     = mat_V;
r.r     = mat_rt;
r.P     = mat_P;
r.sig2r = mat_sig2r;
r.h     = mat_ht;
r.sig2h = mat_sig2h;
