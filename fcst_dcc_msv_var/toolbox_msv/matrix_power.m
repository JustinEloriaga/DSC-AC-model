function B = matrix_power(A, a)
% B = A^a
% where a can be a fraction

% A = [1,0.3;0.3,1];
% a = 0.3;

[U, S, V] = svd(A);
B = U * diag(diag(S).^(a)) * V';

