function updated = dsc_rt_update(previous,y_new)
% The only observation supplied to this function is the new one.
saved_rng=rng; cleanup=onCleanup(@()rng(saved_rng));
rng(previous.rng_state);
particles=dsc_rt_append(previous.particles);
T=size(particles.B,3);
new_loglik=dsc_rt_loglik(particles,y_new,T);
[weights,log_weights,log_predictive]=dsc_rt_normalize(previous.log_weights+new_loglik);

% Update each conditional component using just its new state increments.
hyper=previous.hyper;
dB=particles.B(:,:,T)-particles.B(:,:,T-1);
hyper(:,1)=hyper(:,1)+1;
hyper(:,2:4)=hyper(:,2:4)+[dB(:,1).^2,dB(:,1).*dB(:,2),dB(:,2).^2];
hyper(:,[5 8])=hyper(:,[5 8])+1/2;
hyper(:,6:7)=hyper(:,6:7)+(particles.h(:,:,T)-particles.h(:,:,T-1)).^2/2;
hyper(:,9)=hyper(:,9)+(particles.r(:,T)-particles.r(:,T-1)).^2/2;
updated=struct('particles',particles,'weights',weights,'log_weights',log_weights, ...
    'hyper',hyper,'prior',previous.prior,'rng_state',rng, ...
    'log_evidence',previous.log_evidence+log_predictive,'log_predictive',log_predictive);
end
