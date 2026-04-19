function r=tvsvar_modified_msv2_gam2_gen(y,lags,T0,T0B,T0A,T0H,kB,kA,kH,M,nreport,nthins);

% Generalized version of tvsvar_modified_msv2_gam2.m
% Works for any number of variables (not hardcoded to 4)
%
% _msv, option1 : original D*R*D' decomposition, Wishart process for R
% _msv2, option2 : Hansen's parameterization of R, random-walk for each element

if (nargin ~= 12); error('Wrong # of arguments'); end

% =========================================================================
% PRELIMINARY STEPS
% putting the VAR in the form Y=ZB+e and adding the constant to the regressors
if lags < 0 || lags ~= floor(lags)
    error('lags must be a nonnegative integer');
end

nobs = size(y,1);
x = ones(nobs,1);
for i=1:lags
    x = [x, lag(y,i)];
end

x = x(lags+1:end,:);
y = y(lags+1:end,:);
[T,n]=size(y);
x0 =[x(1:T0,:)];
x  =[x(T0+1:end,:)];
y0 =y(1:T0,:);
y  =y(T0+1:end,:);
Y0 =reshape(y0,n*T0,1);
Y  =reshape(y,n*(T-T0),1);
Z0 =kron(eye(n),x0);
Z  =kron(eye(n),x);

% defining the relevant dimensions
T  =T-T0;         % length of estimation sample
k  =n+lags*n^2;   % # of VAR coefficients
na =n*(n-1)/2;   % # of free coefficients in the matrix A
ns =k+na+n;      % # of shocks with time-varying variance

% =========================================================================
% storage matrices
B_old = zeros(T,k,1);
V_old = zeros(k,k,1);


% =========================================================================
% SETTING UP THE PRIORS
% time varying coefficients in B
priorB=olsblock(y0,x0);
Bbar=reshape(priorB.bhatols,k,1)';
resB=priorB.resols;
VBbar=kron(resB'*resB/length(resB),eye(k/n)/priorB.XX);

% auxiliary vector to be used below
auxMnew=reshape([1:(n-1)^2],n-1,n-1)'; auxMnew=tril(auxMnew);
auxVnew=reshape(auxMnew',(n-1)^2,1); auxVnew(auxVnew==0)=[]; auxVnew=auxVnew';

% =========================================================================
% INIZIALIZATION OF THE ALGORITHM
% -----
% Minchul's modification -> start from the prior mode for V
iwishmean = @(V,T) ( V/(T+size(V,1)+1) ); %mode
% -----
% Drawing the initial V from the prior
V0B=VBbar*T0B*kB^2;
V_old(1:k,1:k) = iwishmean(V0B, T0B); %mean


% =========================================================================
% =========================================================================
%% Setting up for msv-related parameters

m = n; %notation changes...

% ---
% Prior for sig2h

% Cogley and Sargent prior (correct one)
prior_sig2h.T0 = ones(m,1)*(1 / 2);
prior_sig2h.V0 = ones(m,1)*(0.01^2 / 2 );


% ---
% Prior for h(0)
mh0 = zeros(m,1);
Vh0_scale = 10*ones(m,1);

% ---
% Setting up the slice sampler for H(t) — generalized loop over m variables

slice.nobs = T;

% Build slice_chol_cov_z and slice_mean_z for all m variables
slice_chol_cov_z = zeros(T, T, m);
slice_mean_z = zeros(T, m);
for vvind = 1:m
    var_z = Vh0_scale(vvind) + (1:1:T)';
    cov_z = tril(repmat(var_z', T, 1), 0) + tril(repmat(var_z', T, 1), -1)';
    slice_chol_cov_z(:,:,vvind) = chol(cov_z,'lower');
    slice_mean_z(:,vvind) = mh0(vvind)*ones(T,1);
end

% -----
% Dynamic Stochastic Correlation parameters

n_r = m*(m-1)/2; % # of correlation elements
Vr0_scale = 10*ones(n_r,1);
mr0 = 0*ones(n_r,1);

% r1
vvind = 1;
var_z = Vr0_scale(vvind) + (1:1:T)';
cov_z = tril(repmat(var_z', T, 1), 0) + tril(repmat(var_z', T, 1), -1)';
chol_cov_z_r1 = chol(cov_z,'lower');
mean_z_r1 = mr0(vvind)*ones(T,1);

% beyond r1
slice_chol_cov_z_r = chol_cov_z_r1;
slice_mean_z_r = mean_z_r1;
for vvind = 2:n_r
    var_z = Vr0_scale(vvind) + (1:1:T)';
    cov_z = tril(repmat(var_z', T, 1), 0) + tril(repmat(var_z', T, 1), -1)';
    chol_cov_z_r0 = chol(cov_z,'lower');
    mean_z_r0 = mr0(vvind)*ones(T,1);

    % stack
    slice_chol_cov_z_r = cat(3, slice_chol_cov_z_r, chol_cov_z_r0);
    slice_mean_z_r = cat(2, slice_mean_z_r, mean_z_r0);
end

% -------------------
% For sig2r
% Prior for sig2r
prior_sig2r.T0 = ones(n_r,1)*5; %first arg (as in Primiceri, but independent prior)
prior_sig2r.V0 = ones(n_r,1)*(0.01*5); %second arg ((2*0.01^2), as in Primiceri, but independent prior)


%% Initialization

% OLS residual
yhat = olsblock(y,x).resols; %residual

% sig2h
sig2h_old = var(olsblock(y,x).resols) / T ;

% Or, setting the initial H as the log-variance of the OLS residual
ht_old = zeros(T,m);
for t= 1:T
    ht_old(t,:) = log(var(yhat)) + sqrt(sig2h_old).*randn(1,m);
end

% ---
% Initialize correlation at OLS estimate
[temp_v1, temp_v2] = eig(corr(yhat), 'vector');
temp_v3 = temp_v1*diag(log(temp_v2))*temp_v1';
temp_r  = temp_v3(tril(ones(m), -1)==1);

sig2r_old = abs(temp_r')/T;
rt_old = zeros(T,n_r);
for t=1:T
    rt_old(t,:) = temp_r' + sqrt(sig2r_old).*randn(1,n_r);
end

Pt_old = zeros(m,m,T);
for t=1:T
    Pt_old(:,:,t) = veclAtoC(rt_old(t,:));
end



%% Storage matrix

saveind = 1;
if nthins < 1 || nthins ~= floor(nthins)
    error('nthins must be a positive integer');
end

% We save only on iterations where rem(sind, nthins) == 0, so the
% storage size must match that exact count even when M is not divisible by nthins.
nsave = floor(M / nthins);

if nsave < 1
    error('M and nthins imply zero saved draws; increase M or reduce nthins');
end

mat_ht = zeros(T, m, nsave);
mat_lik = zeros(nsave,1);
mat_sig2h = zeros(nsave,m);
mat_rt = zeros(T, n_r, nsave);
mat_sig2r = zeros(nsave,n_r);
mat_PT = zeros(m,m,1,nsave);

mat_B=zeros(1,k,nsave);
mat_V=zeros(k,k,nsave); %V = var(errB)

%% Actual run
% GIBBS SAMPLING ALGORITHM
time_all = tic;
time_report = tic;

i = [];
for sind = 1:M

    % STEP 1: DRAWS of B (Primiceri and Del Negro)

    % Kalman filter
    SHAT = zeros(T, k);
    SIG  = zeros(T, k, k);
    sig  = 4*VBbar;
    shat = Bbar';
    for t=1:T

        % =================================
        % reduced-form covariance matrix
        OMEGA_t = diag(exp(ht_old(t,:)/2))*Pt_old(:,:,t)*diag(exp(ht_old(t,:)/2));
        % =================================
        [shat,sig] = kfilter(y(t,:)',Z([0:n-1]*T+t,:),eye(k),shat,sig,...
            OMEGA_t, squeeze(V_old(1:k,1:k)),t);
        SHAT(t,:)=shat';
        SIG(t,:,:)=sig;
    end

    % simulation smoother
    B_old(T,:)=mvnrnd_modified(shat,sig/2+sig'/2,1);
    for t=T-1:-1:1
        [btTp,StTp]=kback(SHAT(t,:),squeeze(SIG(t,:,:)),squeeze(B_old(t+1,:)),eye(k),squeeze(V_old(1:k,1:k)));
        B_old(t,:)=mvnrnd_modified(btTp,StTp/2+StTp'/2,1);
    end

    % =====================================================================
    % STEP 2: DRAWS OF V (i.e., var(errB), Primiceri and Del Negro)

    % Drawing var(errB)
    errB=[B_old(2:end,:)-B_old(1:end-1,:)];
    V1B=errB'*errB+V0B;
    V_old(1:k,1:k)=iwishrnd( (V1B + V1B')/2.0 , T-1+T0B);


    % =====================================================================
    % STEP 3: DRAWS of R(t), H(t), all other related parameters
    Yhat=(Y-sum((Z.*repmat(squeeze(B_old(:,:)),n,1))')');
    yhat=reshape(Yhat,T,n);

    et_old = yhat;

    % ---
    % Dynamic correlation
    % prep-slice sampler (update parts that change over time)
    lik_old = loglike_yt_given_rt_ht(et_old, ht_old, rt_old, n_r, rt_old(:,n_r)); %have to evaluate it again as et_old has changed
    for vvind = 1:n_r
        slice.scale_z = sqrt(sig2r_old(vvind));
        slice.mean    = slice_mean_z_r(:,vvind);
        slice.chol_cov_z = slice_chol_cov_z_r(:,:,vvind);
        slice.fcn_lik = @(rt_vvind) loglike_yt_given_rt_ht(et_old,ht_old,rt_old,vvind,rt_vvind);

        [rt_vvind_old, lik_old, n_try_r] = slice_sampling_v02(slice, rt_old(:,vvind), lik_old);
        rt_old(:,vvind) = rt_vvind_old;
    end
    % from rt  to Pt
    for t=1:T
        Pt_old(:,:,t) = veclAtoC(rt_old(t,:));
    end

    % Drawing var(errR)
    errR = rt_old(2:end,:)-rt_old(1:end-1,:);
    for vvind = 1:n_r

        % auxiliary part
        rt_old_vvind    = rt_old(1,vvind);
        sig2r_old_vvind = sig2r_old(vvind);

        % proposal
        EE = errR(:,vvind);
        T1 = size(errR,1)/2 + prior_sig2r.T0(vvind);
        V1 = EE'*EE/2 + prior_sig2r.V0(vvind);
        sig2r_new_vvind = gamrnd(T1, V1^(-1))^(-1);

        % mh ratio in log
        logalp = lognormpdf(rt_old_vvind, mr0(vvind), sqrt( (Vr0_scale(vvind)+1)*sig2r_new_vvind) ) ...
            - lognormpdf(rt_old_vvind, mr0(vvind), sqrt( (Vr0_scale(vvind)+1)*sig2r_old_vvind) );
        if log(rand) < logalp %accept
            sig2r_old(vvind) = sig2r_new_vvind;
        else % reject
            sig2r_old(vvind) = sig2r_old_vvind;
        end

    end


    % ---
    % Dynamic variance
    % draw ht (log volatility)
    % preparing slice sampler (update parts that change over time)
    lik_old = loglike_yt_given_ht_Pt(et_old, ht_old, Pt_old, m, ht_old(:,m)); %have to evaluate likelihood again with new Pt_old
    for vvind = 1:m
        slice.scale_z = sqrt(sig2h_old(vvind));
        slice.mean = slice_mean_z(:,vvind);
        slice.chol_cov_z = slice_chol_cov_z(:,:,vvind);
        slice.fcn_lik = @(ht_vvind) loglike_yt_given_ht_Pt(et_old, ht_old, Pt_old, vvind, ht_vvind);

        [ht_vvind_old, lik_old, n_try] = slice_sampling_v02(slice, ht_old(:,vvind), lik_old);
        ht_old(:,vvind) = ht_vvind_old;
    end

    % Drawing var(errH)
    errH = ht_old(2:end,:)-ht_old(1:end-1,:);
    for vvind = 1:m

        % auxiliary part
        ht_old_vvind    = ht_old(1,vvind);
        sig2h_old_vvind = sig2h_old(vvind);

        % proposal
        EE = errH(:,vvind);
        T1 = size(errH,1)/2 + prior_sig2h.T0(vvind);
        V1 = EE'*EE/2 + prior_sig2h.V0(vvind);
        sig2h_new_vvind = gamrnd(T1, V1^(-1))^(-1);

        % mh ratio in log
        logalp = lognormpdf(ht_old_vvind, mh0(vvind), sqrt( (Vh0_scale(vvind)+1)*sig2h_new_vvind) ) ...
            - lognormpdf(ht_old_vvind, mh0(vvind), sqrt( (Vh0_scale(vvind)+1)*sig2h_old_vvind) );
        if log(rand) < logalp %accept
            sig2h_old(vvind) = sig2h_new_vvind;
        else % reject
            sig2h_old(vvind) = sig2h_old_vvind;
        end

    end

    % ===================================
    % Store matrices
    if rem(sind, nreport) == 0

        elapsed_total = toc(time_all);
        elapsed_report = toc(time_report);
        remaining_draws = M - sind;
        est_time_left = elapsed_report * remaining_draws / nreport;

        disp(['sind = ', num2str(sind), ' / ', num2str(M)]);
        disp(['time per report (', num2str(nreport), ' draws, sec.) = ', num2str(elapsed_report)]);
        disp(['elapsed total (sec.) = ', num2str(elapsed_total)]);
        disp(['estimated time left (sec.) = ', num2str(est_time_left), ...
            ' (~', num2str(est_time_left/60), ' min)']);

        time_report = tic;
    end

    if rem(sind,nthins) == 0

        mat_B(:,:,saveind) = B_old(end,:);
        mat_V(:,:,saveind) = V_old;

        mat_rt(:,:,saveind) = rt_old;
        mat_PT(:,:,:,saveind) = Pt_old(:,:,end);

        mat_PT_new(:,:,:,saveind) = Pt_old;


        mat_sig2r(saveind,:) = sig2r_old;

        mat_ht(:,:,saveind)  = ht_old;
        mat_lik(saveind,1) = lik_old;
        mat_sig2h(saveind,:) = sig2h_old;

        saveind = saveind + 1;
    end
end

% Build a output variables
r.B = mat_B;
r.V = mat_V;

r.r = mat_rt;
%r.P = mat_PT(:,:,end,:);
r.P = mat_PT_new;
r.sig2r = mat_sig2r;

r.h = mat_ht;
r.sig2h = mat_sig2h;
