function [shatnew,signew]=kfilter(y,H,F,shat,sig,R,Q,cur_t)

% for r02
% to fix the bug in the original paper (Primiceri),"initial B(t) depends on
% V as well; to update V correctly we need to take care of the initial
% value as well";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL
% y(t)=H*s(t)+e(t)
% s(t)=F*s(t-1)+v(t)
% V(e(t))=R
% V(v(t))=Q
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=y(:); 

% forecasting
if cur_t == 1
    sfor = shat;
    omega = sig;
else
    sfor=F*shat;
    omega=F*sig*F'+Q;
end

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