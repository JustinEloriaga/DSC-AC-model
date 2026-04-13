function r=tvsvar_modified_msv_gam0(y,lags,T0,T0B,T0A,T0H,kB,kA,kH,M,nreport);

% _msv, option1 : original D*R*D' decomposition, Wishart process for R
% _msv_gam0: original prior of Asai and McAleer
% _msv2, option2 : Hansen's parameterization of R, random-walk for each element

if (nargin ~= 11); error('Wrong # of arguments'); end

% =========================================================================
% PRELIMINARY STEPS
% putting the VAR in the form Y=ZB+e and adding the constant to the regressors
x=[];
for i=1:lags
    x=[x lag(y,i)];
end
x=[ones(length(x),1) x];
x=x(lags+1:end,:);
y=y(lags+1:end,:);
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
B=zeros(T,k,M);
V=zeros(k,k,M+1); %V = var(errB)

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
% iwishmean = @(V,T) ( V/(T-size(V,1)-1) ); %mean (may not be well-defined depending
% on df)
iwishmean = @(V,T) ( V/(T+size(V,1)+1) ); %mode
% -----
% Drawing the initial V from the prior
V0B=VBbar*T0B*kB^2;
% V(1:k,1:k,1)=iwishrnd(V0B,T0B);
V(1:k,1:k,1) = iwishmean(V0B, T0B); %mean


% =========================================================================
% =========================================================================
%% Setting up for msv-related parameters

m = n; %notation changes...

% ---
% Prior for sig2h
% prior_sig2h.T0 = ones(m,1)*2; %first arg (2, as in Primiceri, but independent prior)
% prior_sig2h.V0 = ones(m,1)*(2*0.01^2); %second arg ((2*0.01^2), as in Primiceri, but independent prior)


sig2h0 = var(olsblock(y0,x0).resols) / T0;

prior_sig2h.T0 = ones(m,1)*2; %first arg (2, as in Primiceri, but independent prior)
% prior_sig2h.V0 = ones(m,1)*(sig2h0*(2-1)); %second arg ((2*0.01^2), as in Primiceri, but independent prior)
prior_sig2h.V0 = (sig2h0*(2-1))'; %second arg ((2*0.01^2), as in Primiceri, but independent prior)

% ---
% Prior for h(0)
mh0 = zeros(m,1);
Vh0_scale = 10*ones(m,1);

% ---
% Setting up the slice sampler for H(t)

slice.nobs = T;

% h1
vvind = 1;
var_z = Vh0_scale(vvind) + (1:1:T)';
cov_z = tril(repmat(var_z', T, 1), 0) + tril(repmat(var_z', T, 1), -1)';
chol_cov_z_h1 = chol(cov_z,'lower');
mean_z_h1 = mh0(vvind)*ones(T,1);

% h2
vvind = 2;
var_z = Vh0_scale(vvind) + (1:1:T)';
cov_z = tril(repmat(var_z', T, 1), 0) + tril(repmat(var_z', T, 1), -1)';
chol_cov_z_h2 = chol(cov_z,'lower');
mean_z_h2 = mh0(vvind)*ones(T,1);

% h3
vvind = 3;
var_z = Vh0_scale(vvind) + (1:1:T)';
cov_z = tril(repmat(var_z', T, 1), 0) + tril(repmat(var_z', T, 1), -1)';
chol_cov_z_h3 = chol(cov_z,'lower');
mean_z_h3 = mh0(vvind)*ones(T,1);

% h4
vvind = 4;
var_z = Vh0_scale(vvind) + (1:1:T)';
cov_z = tril(repmat(var_z', T, 1), 0) + tril(repmat(var_z', T, 1), -1)';
chol_cov_z_h4 = chol(cov_z,'lower');
mean_z_h4 = mh0(vvind)*ones(T,1);

% putting everything together
slice_chol_cov_z = cat(3, chol_cov_z_h1, chol_cov_z_h2, chol_cov_z_h3, chol_cov_z_h4);
slice_mean_z = [mean_z_h1, mean_z_h2, mean_z_h3, mean_z_h4];

% -----
% Dynamic Stochastic Correlation parameters
% initial Q00
temp_ols = olsblock(y0,x0);
H0hat = diag(std(temp_ols.resols));
Q00 = (H0hat)*corr(temp_ols.resols)*(H0hat);

% prior for Ainv (following Asai and McAleer)
prior_Ainv.gam  = m;
prior_Ainv.Cinv = (eye(m)) * (prior_Ainv.gam);

% prior d -> uniform over [L,U]
prior_d.L = -1;
prior_d.U = 1;
prior_d.mhtune_a = 0.3; %target acceptance rate
prior_d.mhtune_c = 0.55; %cooling rate (between 0.5 and 1), c = inf means no adaptation
prior_d.mhtune_b = 10000; %bound
prior_d.mhtune_logvar = log(0.1^2); %log proposal variance (initial proposal variance)

% prior k -> truncated exponential
prior_k.lam0 = 5;
prior_k.mhtune_a = 0.3; %target acceptance rate
prior_k.mhtune_c = 0.55; %cooling rate (between 0.5 and 1), c = inf means no adaptation
prior_k.mhtune_b = 10000; %bound
prior_k.mhtune_logvar = log(0.1^2); %log proposal variance (initial proposal variance)

%% Initialization


% % drawing the initial H from the prior
% ht_old = zeros(T,m);
% h0 = mh0 + diag(sqrt(Vh0_scale .* sig2h_old))*randn(m,1);
% for t = 1:T
%     h0 = h0 + diag(sqrt(sig2h_old))*randn(m,1);
%     ht_old(t,:) = h0;
% end

% OLS residual
yhat = olsblock(y,x).resols; %residual

% sig2h
% sig2h_old = prior_sig2h.V0 ./ (prior_sig2h.T0 - 1); % start from the prior mean
sig2h_old = var(olsblock(y,x).resols) / T ;

% Or, setting the initial H as the log-variance of the OLS residual
ht_old = zeros(T,m);
for t= 1:T
    ht_old(t,:) = log(var(yhat)) + sqrt(sig2h_old).*randn(1,m);
end

% ---
% Initialize correlation at OLS estimate
% Pt_old = zeros(m,m,T);
% Qt_old = zeros(m,m,T);
% Qinvt_old = zeros(m,m,T);
% for t=1:T
%     Pt_temp = corr(yhat);
%     Pt_old(:,:,t) = Pt_temp;
%     Qt_old(:,:,t) = diag(exp(ht_old(t,:)/2)) * Pt_temp * diag(exp(ht_old(t,:)/2));
%     Qinvt_old(:,:,t) = eye(m) / Qt_old(:,:,t);
% end
% % ---

d_old = 0.5;
k_old = 15;
Ainv_old = eye(m);
[Pt_old, Qt_old, Qinvt_old] = gen_data_sdc(T, inv(Ainv_old), d_old, k_old, Q00); %initialize

%% Storage matrix
mat_ht = zeros(T, m, M);
mat_lik = zeros(M,1);
mat_sig2h = zeros(M,m);

% mat_Pt = zeros(m,m,T,M);
% mat_Qt = zeros(m,m,T,M);
mat_QT = zeros(m,m,1,M);
mat_Q_rej = zeros(M,T);


mat_Ainv = zeros(m,m,M);
mat_d = zeros(M,1);
mat_d_rej = zeros(M,1);
mat_k = zeros(M,1);
mat_k_rej = zeros(M,1);


%% Actual run
% GIBBS SAMPLING ALGORITHM
time_all = tic;
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
        % OMEGA_t = tria(A(t,:,i))\diag(exp(2*squeeze(H(t,:,i))))/(tria(A(t,:,i))');
        OMEGA_t = diag(exp(ht_old(t,:)/2))*Pt_old(:,:,t)*diag(exp(ht_old(t,:)/2));
        % =================================
        [shat,sig] = kfilter(y(t,:)',Z([0:n-1]*T+t,:),eye(k),shat,sig,...
            OMEGA_t, squeeze(V(1:k,1:k,sind)));
        SHAT(t,:)=shat';
        SIG(t,:,:)=sig;
    end
    
    % simulation smoother
    B(T,:,sind)=mvnrnd_modified(shat,sig/2+sig'/2,1);
    for t=T-1:-1:1
        [btTp,StTp]=kback(SHAT(t,:),squeeze(SIG(t,:,:)),squeeze(B(t+1,:,sind)),eye(k),squeeze(V(1:k,1:k,sind)));
        B(t,:,sind)=mvnrnd_modified(btTp,StTp/2+StTp'/2,1);
    end
    
    % =====================================================================
    % STEP 2: DRAWS OF V (i.e., var(errB), Primiceri and Del Negro)
    
    % Drawing var(errB)
    errB=[B(2:end,:,sind)-B(1:end-1,:,sind)];
    V1B=errB'*errB+V0B;
    %V(1:k,1:k,i+1)=iwishrnd(V1B,T-1+T0B);
    V(1:k,1:k,sind+1)=iwishrnd( (V1B + V1B')/2.0 , T-1+T0B);
    
    
    % =====================================================================
    % STEP 3: DRAWS of R(t), H(t), all other related parameters
    Yhat=(Y-sum((Z.*repmat(squeeze(B(:,:,sind)),n,1))')');
    yhat=reshape(Yhat,T,n);
    
    %yhat = olsblock(y,x).resols; % OLS
    
    et_old = yhat;
    zt_old = et_old ./ exp(ht_old/2);
    
    % draw Qt,Pt given zt and Ainv (MH)
    [Qt_old, Qinvt_old, Q_rej] = gen_Qt(Qt_old, Qinvt_old, d_old, k_old, Ainv_old, Q00, zt_old);
    for t=1:T
        Pt_old(:,:,t) = QtoP(Qt_old(:,:,t));
    end
    
    % draw Ainv given data
    Ainv_old = gen_Ainv(d_old, k_old, Q00, Qt_old, Qinvt_old, prior_Ainv);
    
    % drawing d given data
    [d_old, rej_d, prior_d] = gen_d(Qt_old, Qinvt_old, d_old, k_old, Ainv_old, Q00, prior_d, sind);
%         d_old = 0.3;
%         rej_d = 1;
    
    % drawing k given data
    [k_old, rej_k, prior_k] = gen_k(Qt_old, Qinvt_old, d_old, k_old, Ainv_old, Q00, prior_k, sind);
%         k_old = 20;
%         rej_k = 1;
    
    
    % ---
    % Dynamic variance
    % draw ht (log volatility)
    % preparing slice sampler (update parts that change over time)
    lik_old = loglike_yt_given_ht_Pt(et_old, ht_old, Pt_old, 4, ht_old(:,4)); %have to evaluate likelihood again with new Pt_old
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
        EE = errH(:,vvind);
        T1 = size(errH,1)/2 + prior_sig2h.T0(vvind);
        V1 = EE'*EE/2 + prior_sig2h.V0(vvind);
        sig2h_old(vvind) = gamrnd(T1, V1^(-1))^(-1);
    end
    
    % ===================================
    % Store matrices
    if rem(sind, nreport) == 0
        disp(['sind = ', num2str(sind), ' / ', num2str(M)]);
        time_per_draw = toc(time_all)/sind;
        disp(['time per draw (sec.) = ', num2str(time_per_draw)]);
    end
    
%     mat_Pt(:,:,:,sind) = Pt_old;
%     mat_Qt(:,:,:,sind) = Qt_old;
    mat_QT(:,:,:,sind) = Qt_old(:,:,end);
    
    mat_Ainv(:,:,sind) = Ainv_old;
    mat_Q_rej(sind,:) = Q_rej;
    
    mat_d_rej(sind,:) = rej_d;
    mat_d(sind,1) = d_old;
    
    mat_k_rej(sind,:) = rej_k;
    mat_k(sind,1) = k_old;
    
    mat_ht(:,:,sind)  = ht_old;
    mat_lik(sind,1) = lik_old;
    mat_sig2h(sind,:) = sig2h_old;
    
end

% Build a output variables
r.B = B(end,:,:);
r.V = V;

r.Q = mat_QT(:,:,end,:);
r.k = mat_k;
r.d = mat_d;
r.Ainv = mat_Ainv;
r.h = mat_ht(end,:,:);
r.sig2h = mat_sig2h;

























