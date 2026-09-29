function [weights,log_weights,log_total] = dsc_rt_normalize(log_values)
peak=max(log_values);
if ~isfinite(peak)
    error('rt:Weights','All particle likelihoods are numerically zero.');
end
log_total=peak+log(sum(exp(log_values-peak)));
log_weights=log_values-log_total;
weights=exp(log_weights);
end
