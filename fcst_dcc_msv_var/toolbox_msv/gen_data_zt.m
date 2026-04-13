function zt = gen_data_zt(Pt)
% zt ~ N(0, Pt)

[m, ~, T] = size(Pt);

zt = zeros(T,m);
for t=1:T
    zt(t,:) = chol(Pt(:,:,t), 'lower') *randn(m,1);
end