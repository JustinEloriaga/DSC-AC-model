function crps = compute_crps2(y,X)
% second way to compute crps
% y: actual
% X: is a vector that collects draws from the predictive distribution
m = size(X,1);
X = sort(X, 'ascend');
crps = 2/(m^2) * sum( (X-y).*(m*(y<X)-(1:1:m)'+1/2) );