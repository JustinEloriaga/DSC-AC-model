function [C,stats,x] = dsc_corr_matrix(r,x)
% Inverse off-diagonal matrix-log map, with explicit convergence checks.
n=(1+sqrt(1+8*numel(r)))/2;
if n~=floor(n), error('dsc:CorrelationDimension','Invalid triangular dimension.'); end
if nargin<2||isempty(x), x=zeros(n,1); end
sel=tril(true(n),-1); diagonal=1:n+1:n*n;
A=zeros(n); A(sel)=r; A=A+A';
tol=1e-8*sqrt(n); converged=false;
for iteration=1:200
    A(diagonal)=x;
    [Q,L]=eig(A,'vector'); E=Q*diag(exp(L))*Q';
    delta=log(real(diag(E)));
    if any(~isfinite(delta)), error('dsc:CorrelationNumerical','Nonfinite correlation transform.'); end
    x=x-delta;
    if norm(delta)<tol, converged=true; break; end
end
if ~converged, error('dsc:CorrelationConvergence','Correlation inverse did not converge after 200 iterations.'); end
A(diagonal)=x;
[Q,L]=eig(A,'vector'); C=real(Q*diag(exp(L))*Q');
C=(C+C')/2; C(diagonal)=1;
[~,bad]=chol(C);
if bad, error('dsc:CorrelationPD','Correlation matrix is not numerically positive definite.'); end
stats=struct('iterations',iteration,'converged',true,'residual',norm(delta));
end
