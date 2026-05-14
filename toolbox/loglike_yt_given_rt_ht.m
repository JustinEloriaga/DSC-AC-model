function loglik = loglike_yt_given_rt_ht(yt,ht,rt,rrind,rt_prop)
% for yt ~ N(0, diag(exp(ht/2)) * Pt(rt) * diag(exp(ht/2)))
%
% Optimized:
%  - precompute sd = exp(ht/2) once outside the t-loop
%  - build Sigma_t via outer product (sd_t' * sd_t) .* Pt_t
%    (avoids two diagonal-matrix multiplications per t)
%  - warm-start veclAtoC's fixed-point iteration across t (rt is a
%    random walk so adjacent t solutions are close)
%  - one batched log_mvnpdf call with 3D Sigma instead of T per-row calls

% update rrind-th column of rt with rt_prop
rt(:,rrind) = rt_prop;

% size info
[T, m] = size(ht);

% precompute exp(ht/2) once
sd = exp(ht/2);    % T x m

% build Pt stack with warm-start across t (random-walk rt -> close adjacent solutions)
Pt     = zeros(m, m, T);
x0_wrm = [];
for t = 1:T
    [Pt(:,:,t), x0_wrm] = veclAtoC(rt(t,:), x0_wrm);
end

% vectorized Sigma(:,:,t) = (sd_t' * sd_t) .* Pt(:,:,t)
sdcol = reshape(sd', m, 1, T);   % m x 1 x T
sdrow = reshape(sd', 1, m, T);   % 1 x m x T
Sigma = (sdcol .* sdrow) .* Pt;
Sigma = (Sigma + permute(Sigma, [2 1 3])) / 2;

% one batched log_mvnpdf call (uses 3D Sigma branch internally)
try
    loglik = sum(log_mvnpdf(yt, zeros(1, m), Sigma));
catch
    disp('**** given_rt_ht, log_mvnpdf error ...');
    loglik = -Inf;
end
