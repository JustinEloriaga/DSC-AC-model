function [Pt, Qt,Qinvt] = gen_data_sdc(T, A, d, k, Q00)
% stochastic dynamic correlation

% T = 1000;
% m = size(A,1);
% 
% Q00      = eye(m); %we assume that this is known
% Q00(1,2) = 0.9;
% Q00(2,1) = 0.9;

% correlation matrix
m = size(A,1);
Pt = zeros(m,m,T); %correlation matrix
Qt = zeros(m,m,T);
Qinvt = zeros(m,m,T);
Q0 = Q00;
for t=1:1:T
    
    % update
    Qmdh = matrix_power(Q0, -d/2);
    S1 = (1/k) * Qmdh * A * Qmdh';
    Qinv1 = wishrnd(S1, k);
    
    % correlation
    Q1 = eye(m)/Qinv1;
    P1 = QtoP(Q1);
    
    % store
    Pt(:,:,t) = P1;
    Qt(:,:,t) = Q1;
    Qinvt(:,:,t) = Qinv1;
    
    % update for tomorrow
    Q0 = Q1;
end