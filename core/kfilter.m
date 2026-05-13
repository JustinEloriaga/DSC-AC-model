function [shatnew,signew]=kfilter(y,H,F,shat,sig,R,Q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL
% y(t)=H*s(t)+e(t)
% s(t)=F*s(t-1)+v(t)
% V(e(t))=R
% V(v(t))=Q
%
% At t=1 the inputs (shat, sig) represent the prior on the previous
% state (e.g. B_0), so sfor = F*shat = E[B_1|0] and
% omega = F*sig*F' + Q = Var[B_1|0].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=y(:);

% forecasting (single branch, valid for all t)
sfor  = F*shat;
omega = F*sig*F' + Q;

% updating
sigma=H*omega*H'+R;

% original
% sigma = (sigma + sigma')/2.0;
k=omega*H'/sigma; 

% % modified
% k=omega*H'*(invChol_mex(sigma)); 

ferr=y-H*sfor; 
shatnew=sfor+k*ferr;
signew=omega-k*H*omega;