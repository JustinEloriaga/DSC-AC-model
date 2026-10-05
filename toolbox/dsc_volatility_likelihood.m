function ll=dsc_volatility_likelihood(x,h,residual,z,cross,precision,observed,base_ll)
% Exact Gaussian log-likelihood change when one log-volatility path moves.
delta=residual(observed).*exp(-x(observed)/2)-z(observed);
change=sum(x(observed)-h(observed))+ ...
    sum(2*delta.*cross(observed)+delta.^2.*precision(observed));
ll=base_ll-.5*change;
if ~isfinite(ll), ll=-Inf; end
end
