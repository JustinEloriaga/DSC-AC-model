function P = QtoP(Q)
% correlation matrix P is defined by normalizing the positive definite
% matrix Q
P = diag(diag(Q).^(-1/2)) *Q * diag(diag(Q).^(-1/2));