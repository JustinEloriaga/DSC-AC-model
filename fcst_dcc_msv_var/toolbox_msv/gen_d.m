function [d_old, rej_d, prior_d] = gen_d(Qt_old, Qinvt_old, d_old, k_old, Ainv_old, Q00, prior_d, sind)

% MH generation of d

% size info
[~, ~, T] = size(Qt_old);

% proposal
d_prop = d_old + exp(prior_d.mhtune_logvar/2)*randn;

% proposal density at the proposal and the old
Q0 = Q00;
Qinv1 = Qinvt_old(:,:,1);
psi = LogAbsDet(Q0);

Qdh_prop = matrix_power(Q0, d_prop/2);
Cinv_prop = Qdh_prop * Qinv1 * Qdh_prop';
Qdh_old = matrix_power(Q0, d_old/2);
Cinv_old = Qdh_old * Qinv1 * Qdh_old';
for t=2:T
    Q0 = Qt_old(:,:,t-1);
    Qinv1 = Qinvt_old(:,:,t);
    psi = psi + LogAbsDet(Q0);
    
    Qdh_prop  = matrix_power(Q0, d_prop/2);
    Cinv_prop = Cinv_prop + Qdh_prop * Qinv1 * Qdh_prop';
    Qdh_old   = matrix_power(Q0, d_old/2);
    Cinv_old  = Cinv_old + Qdh_old * Qinv1 * Qdh_old';
    
end
psi = (k_old/2) * psi;
Cinv_prop = k_old*Cinv_prop;
Cinv_old = k_old*Cinv_old;

logf_d_prop = -Inf;
if (d_prop > prior_d.L) && (d_prop < prior_d.U)
    logf_d_prop = psi*d_prop -1/2*trace(Ainv_old*Cinv_prop);
end
logf_d_old = psi*d_old -1/2*trace(Ainv_old*Cinv_old);

% mhratio
logmhr_d = logf_d_prop-logf_d_old;

% accept/reject
rej_d = 0;
if (log(rand) < logmhr_d) %accept
    d_old = d_prop;
else %reject
    rej_d = 1;
end

% Adaptation d
mhratio_d = min(exp(logmhr_d), 1);
prior_d.mhtune_logvar = adaptiveMCMCvar(prior_d.mhtune_logvar, ...
    mhratio_d, prior_d.mhtune_a, prior_d.mhtune_b, prior_d.mhtune_c, sind);