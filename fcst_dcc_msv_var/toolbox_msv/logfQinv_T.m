function logf = logfQinv_T(Qinv1,Q1,d,k,Ainv,zz1)

% _T: for t=T

% log f(Qinv)
% timing
% 0 = t-1
% 1 = t
% 2 = t+2

% See Asai-McAleer (2009) page 191

% For t=T

Pinv1 = QinvtoPinv(Qinv1,Q1);
% Qdh1 = matrix_power(Q1, d/2);

% Sinv1 = k*Qdh1*Ainv*Qdh1;
% term1 = 1/2*LogAbsDet(Pinv1);
% term2 = -(1+k*d)/2 * LogAbsDet(Qinv1);
% term3 = -1/2 * trace(Sinv1*Qinv2);
% term4 = -1/2 * trace( (Pinv1-Qinv1) * zz1);

% logf = term1+term2+term3+term4;

% Wrong correction term in the paper
% logf = -1/2*trace( (Pinv1-Qinv1) * zz1);

% Typo in the paper is corrected by Minchul
% correction term should be:
term0 = 1/2*sum(log(diag(Q1)));
logf = term0 -1/2*trace( (Pinv1-Qinv1) * zz1);
