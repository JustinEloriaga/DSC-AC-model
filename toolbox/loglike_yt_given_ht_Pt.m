function loglik = loglike_yt_given_ht_Pt(yt,ht,Pt,hhind,ht_prop)
% for yt ~ N(0, diag(exp(ht/2)) * Pt * diag(exp(ht/2)))
%
% Optimized:
%  - precompute sd = exp(ht/2) once outside the t-loop
%  - build Sigma_t via outer product (sd_t' * sd_t) .* Pt(:,:,t)
%    (avoids two diagonal-matrix multiplications per t)
%  - one batched log_mvnpdf call with 3D Sigma instead of T per-row calls

% update hhind-th column of ht with ht_prop
ht(:,hhind) = ht_prop;

% size info
[m, ~, T] = size(Pt);

% precompute exp(ht/2) once
sd = exp(ht/2);    % T x m

% vectorized Sigma(:,:,t) = (sd_t' * sd_t) .* Pt(:,:,t)
sdcol = reshape(sd', m, 1, T);   % m x 1 x T
sdrow = reshape(sd', 1, m, T);   % 1 x m x T
Sigma = (sdcol .* sdrow) .* Pt;
Sigma = (Sigma + permute(Sigma, [2 1 3])) / 2;

% one batched log_mvnpdf call
loglik = sum(log_mvnpdf(yt, zeros(1, m), Sigma));
