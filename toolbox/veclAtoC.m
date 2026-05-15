function C = veclAtoC(veclA)
% from veclA  to Correlation matrix
% Peter Hansen's paper

% tuning parameter for the fixed-point
tolX = 1e-8 * sqrt(4); % as in the original paper
Kmax = 200;

% auxiliary
m = numel(veclA);
n = (1+sqrt(1+8*m))/2;
sel_diag = eye(n)==1;

% initialization
A = zeros(n);
A(logical(tril(ones(n),-1))) = veclA;
A = A + A';

% fixed-point algorithm
x0 = zeros(n,1);

i=1;
finish = false;
while (~finish)
    A(sel_diag) = x0;
    [Q,L] = eig(A, 'vector');
    eAx = Q*diag(exp(L))*Q';

    diffx0 = log(diag(eAx));
    diffx0 = real(diffx0);

    x0 = x0 - diffx0;

    if ( i > Kmax ) || (norm(diffx0)<tolX)
        finish = true;
    end

    i = i + 1;
end

% final correlation matrix
A(sel_diag) = x0;
[Q,L] = eig(A, 'vector');
eAx = Q*diag(exp(L))*Q';

% clean-up for numerical reason
C = real(eAx);
C(sel_diag) = 1;
C = (C + C')/2;
