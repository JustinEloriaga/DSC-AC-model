function Pinv = QinvtoPinv(Qinv, Q)
% correlation matrix P is defined by normalizing the positive definite
% matrix Q
Pinv = diag(diag(Q).^(1/2)) *Qinv * diag(diag(Q).^(1/2));