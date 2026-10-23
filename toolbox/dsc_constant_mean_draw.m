function [mu,posterior_mean,posterior_covariance] = dsc_constant_mean_draw(Y,h,P,mask,prior_mean,prior_covariance,z,check)
%DSC_CONSTANT_MEAN_DRAW Exact Gaussian update for one time-invariant mean.
%   The prior is N(prior_mean,prior_covariance). Each observed date contributes
%   its own observed covariance submatrix diag(exp(h/2))*P*diag(exp(h/2)).
%   Missing coordinates and dates contribute no likelihood information.
%   Optional z fixes the standard-normal draw for deterministic validation.
if nargin<7, z=[]; end
if nargin<8 || isempty(check), check=@() []; end
[T,m]=size(Y);
assert(isreal(Y)&&isreal(h)&&isreal(P)&&isequal(size(h),[T m])&&isequal(size(mask),[T m])&& ...
    size(P,1)==m&&size(P,2)==m&&size(P,3)==T, ...
    'dsc:ConstantMeanShape','Mean-update observations and covariance dimensions differ.');
prior_mean=prior_mean(:);
assert(isreal(prior_mean)&&isreal(prior_covariance)&&numel(prior_mean)==m&&all(isfinite(prior_mean))&& ...
    isequal(size(prior_covariance),[m m])&&all(isfinite(prior_covariance(:)))&& ...
    norm(prior_covariance-prior_covariance','fro')<=1e-12*max(1,norm(prior_covariance,'fro')), ...
    'dsc:ConstantMeanPrior','The mean prior must have a finite mean and symmetric covariance.');
[prior_factor,bad]=chol((prior_covariance+prior_covariance')/2);
assert(bad==0,'dsc:ConstantMeanPrior','The mean prior covariance must be positive definite.');
precision=prior_factor\(prior_factor'\eye(m));
rhs=precision*prior_mean;
for t=1:T
    if mod(t,16)==1, check(); end
    ix=find(mask(t,:));
    if isempty(ix), continue; end
    assert(all(isfinite(Y(t,ix)))&&all(isfinite(h(t,ix))), ...
        'dsc:ConstantMeanObserved','Observed mean-update inputs must be finite.');
    sd=exp(h(t,ix)/2);
    covariance=(sd'*sd).*P(ix,ix,t);
    [factor,bad]=chol((covariance+covariance')/2,'lower');
    assert(bad==0&&all(isfinite(factor(:))), ...
        'dsc:ConstantMeanCovariance','An observed covariance is not finite positive definite.');
    whitened=factor\eye(numel(ix));
    precision(ix,ix)=precision(ix,ix)+whitened'*whitened;
    rhs(ix)=rhs(ix)+whitened'*(factor\Y(t,ix)');
end
[factor,bad]=chol((precision+precision')/2);
assert(bad==0,'dsc:ConstantMeanPrecision','The posterior mean precision is not positive definite.');
posterior_mean=factor\(factor'\rhs);
if isempty(z), z=randn(m,1); end
assert(isnumeric(z)&&isreal(z)&&numel(z)==m&&all(isfinite(z(:))), ...
    'dsc:ConstantMeanDraw','The normal draw must be a finite real m-vector.');
mu=posterior_mean+factor\z(:);
assert(all(isfinite(mu)),'dsc:ConstantMeanDraw','The sampled mean must be finite.');
if nargout>2
    posterior_covariance=factor\(factor'\eye(m));
    posterior_covariance=(posterior_covariance+posterior_covariance')/2;
end
check();
end
