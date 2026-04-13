function [k_old, rej_k, prior_k] = gen_k(Qt_old, Qinvt_old, d_old, k_old, Ainv_old, Q00, prior_k, sind)

% MH generation of d

% size info
[m, ~, T] = size(Qt_old);

% proposal
k_prop = k_old + exp(prior_k.mhtune_logvar/2)*randn;

% proposal density at the proposal and the old
Q0 = Q00;
Qinv1 = Qinvt_old(:,:,1);

Qdh_old = matrix_power(Q0, d_old/2);
QQQ_old = Qdh_old * Qinv1 * Qdh_old';
Cinv_old = QQQ_old;
logsumQQQ_old = LogAbsDet(QQQ_old);

for t=2:T
    Q0 = Qt_old(:,:,t-1);
    Qinv1 = Qinvt_old(:,:,t);
    
    Qdh_old   = matrix_power(Q0, d_old/2);
    QQQ_old = Qdh_old * Qinv1 * Qdh_old';
    Cinv_old  = Cinv_old + QQQ_old;
    logsumQQQ_old = logsumQQQ_old + LogAbsDet(QQQ_old);
    
end
Cinv_prop = k_prop*Cinv_old;
Cinv_old  = k_old*Cinv_old;

logf_k_prop = -Inf;
if k_prop > m
    logf_k_prop = -prior_k.lam0*k_prop + T*k_prop*m/2*log(k_prop/2) ...
        +T*k_prop/2*LogAbsDet(Ainv_old) ...
        -T*sum(gammaln( (k_prop+1-(1:1:m)')/2) ) ...
        +k_prop/2*logsumQQQ_old ...
        -1/2*trace(Ainv_old*Cinv_prop);
end

logf_k_old = -prior_k.lam0*k_old + T*k_old*m/2*log(k_old/2) ...
    +T*k_old/2*LogAbsDet(Ainv_old) ...
    -T*sum(gammaln( (k_old+1-(1:1:m)')/2) ) ...
    +k_old/2*logsumQQQ_old ...
    -1/2*trace(Ainv_old*Cinv_old);


% mhratio
logmhr_k = logf_k_prop-logf_k_old;

% accept/reject
rej_k = 0;
if (log(rand) < logmhr_k) %accept
    k_old = k_prop;
else %reject
    rej_k = 1;
end

% Adaptation d
mhratio_k = min(exp(logmhr_k), 1);
prior_k.mhtune_logvar = adaptiveMCMCvar(prior_k.mhtune_logvar, ...
    mhratio_k, prior_k.mhtune_a, prior_k.mhtune_b, prior_k.mhtune_c, sind);
