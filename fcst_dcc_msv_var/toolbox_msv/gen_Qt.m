function [Qt_old, Qinvt_old, Q_rej] = gen_Qt(Qt_old, Qinvt_old, d_old, k_old, Ainv_old, Q00, zt)

% dimension
[T,m] = size(zt);

% MH-step for Qt matrices
Q_rej = zeros(T,1);
Q0 = Q00; %initial Q (is fixed as in Asai and McAleer)
for t=1:T
    
    % construct proposal density
    z1    = zt(t,:)';
    zz1   = z1*z1';
    Qdh0  = matrix_power(Q0, d_old/2);
    Sinv0 = k_old * Qdh0 * Ainv_old * Qdh0';
    Shat0 = eye(m)/(Sinv0  + zz1);
    khat  = k_old + 1;
    
    % generate from the proposal
    Shat0 = (Shat0 + Shat0')/2;
    Qinv1_prop = wishrnd(Shat0, khat);
    Q1_prop    = eye(m)/Qinv1_prop;
    Qinv1_old  = Qinvt_old(:,:,t);
    Q1_old     = Qt_old(:,:,t);
    
    % MH ratio
    if t<T
        % mh ratio for t<T
        Qinv2 = Qinvt_old(:,:,t+1);
        logfQinv_prop = logfQinv(Qinv1_prop, Q1_prop, Qinv2, d_old, k_old, Ainv_old, zz1);
        logfQinv_old  = logfQinv(Qinv1_old, Q1_old, Qinv2, d_old, k_old, Ainv_old, zz1);
        
    else
        % mh ratio for t=T
        logfQinv_prop = logfQinv_T(Qinv1_prop, Q1_prop, d_old, k_old, Ainv_old, zz1);
        logfQinv_old  = logfQinv_T(Qinv1_old, Q1_old, d_old, k_old, Ainv_old, zz1);
    end
    
    mhr = logfQinv_prop-logfQinv_old;
    if log(rand) < mhr %accept
        
        Qt_old(:,:,t) = Q1_prop;
        Qinvt_old(:,:,t) = Qinv1_prop;
    else % reject
        Q_rej(t,1) = 1;
    end
    
    % Update Q0 for tomorrow
    Q0 = Qt_old(:,:,t);
end

