function Ainv_new = gen_Ainv(d, k, Q00, Qt, Qinvt, prior_Ainv)

% size info
[m,~,T] = size(Qt);

% posterior sampling of A^(-1)
Qdh = matrix_power(Q00, d/2);
Cinv = Qdh*Qinvt(:,:,1)*Qdh;

for t=2:T
    Qdh = matrix_power(Qt(:,:,t-1), d/2);
    Cinv = Cinv + Qdh*Qinvt(:,:,t)*Qdh;
end
Cinv = k*Cinv;

Cinvhat = prior_Ainv.Cinv + Cinv;
gamhat = k*T + prior_Ainv.gam;

% generate a draw
Chat = eye(m)/Cinvhat;
Chat = (Chat + Chat')/2;
Ainv_new = wishrnd(Chat, gamhat);