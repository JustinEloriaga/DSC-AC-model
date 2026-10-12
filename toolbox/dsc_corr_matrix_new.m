function [C,stats,x] = dsc_corr_matrix_new(r,x)
%DSC_CORR_MATRIX_NEW Inverse GFT using analytic log-diagonal Newton steps.
% Solve ell(x)=log(diag(expm(A[x;r])))=0 using the exponential's spectral
% Frechet derivative. The dense n-by-n Jacobian reuses the current
% eigensystem, avoiding extra exponential evaluations. This is a small-
% dimension reference inspired by the 2026 Newton method, not its official
% matrix-free PCG implementation. An initial AH phase enters the local Newton
% region. Exact trace normalization uses the same eigenpairs. The
% dsc_corr_matrix contract is preserved, without Newton trial evaluations.
n=(1+sqrt(1+8*numel(r)))/2;
if n~=floor(n), error('dsc:CorrelationDimension','Invalid triangular dimension.'); end
if nargin<2||isempty(x), x=zeros(n,1); end
x=x(:);
sel=tril(true(n),-1); diagonal=1:n+1:n*n;
A=zeros(n); A(sel)=r; A=A+A';
tol=1e-8*sqrt(n); converged=false; residual=Inf;
maxit=200; newton_phase=false;
for iteration=1:maxit
    A(diagonal)=x;
    [Q,L]=eig(A,'vector');
    largest=max(L);
    shift=largest+log(sum(exp(L-largest))/n);
    x=x-shift; L=L-shift;
    w=exp(L);
    d=real((Q.^2)*w);
    ell=log(d);
    if any(~isfinite(ell)), error('dsc:CorrelationNumerical','Nonfinite correlation transform.'); end
    residual=norm(ell);
    if residual<tol, converged=true; break; end
    if ~newton_phase && max(abs(ell))>.25
        x=x-ell;
        continue
    end
    newton_phase=true;
    % H is the derivative of diag(expm(A)) with respect to diag(A).
    % Its log-residual Jacobian is diag(1./d)*H, so solve the symmetric
    % system H*step=-d.*ell without forming that nonsymmetric Jacobian.
    pair_vectors=zeros(n,n*(n+1)/2);
    pair_weights=zeros(1,n*(n+1)/2);
    k=0;
    for a=1:n
        k=k+1;
        pair_vectors(:,k)=Q(:,a).^2;
        pair_weights(k)=w(a);
        for b=a+1:n
            delta=abs(L(a)-L(b));
            if delta==0
                divided_difference=w(a);
            else
                divided_difference=max(w(a),w(b))*(-expm1(-delta)/delta);
            end
            k=k+1;
            pair_vectors(:,k)=Q(:,a).*Q(:,b);
            pair_weights(k)=2*divided_difference;
        end
    end
    H=(pair_vectors.*pair_weights)*pair_vectors';
    step=-H\(d.*ell);
    if any(~isfinite(step)), error('dsc:CorrelationNumerical','Nonfinite correlation transform.'); end
    x=x+step;
end
if ~converged, error('dsc:CorrelationConvergence','New correlation inverse did not converge.'); end
C=real((Q.*w')*Q');
C=(C+C')/2; C(diagonal)=1;
[~,bad]=chol(C);
if bad, error('dsc:CorrelationPD','Correlation matrix is not numerically positive definite.'); end
stats=struct('iterations',iteration,'converged',true,'residual',residual,'backend','matlab_new');
end
