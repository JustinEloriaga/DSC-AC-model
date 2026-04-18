function [btT,StT]=kback(btt,Stt,bnT,A,omega)
%[btT StT]=ksmooth(btt,Stt,bnT,SnT,A,omega)
% Smoothing recursion.  State evolution equation is
%    bn=A*bt+e,  Var(e)=omega
%    bt|t ~ N(btt,Stt) -- from Kalman Filter
%    bn|T ~ N(bnT,SnT) -- distribution of bn given full sample. From 
%                         KF if n=T, otherwise from this recursion
%    bt|T ~ N(btT,StT)
AS=A*Stt;
G=AS*A'+omega;


% original
% G = (G + G')/2.0;
SAGI=AS'/G;

% modified
% opts.POSDEF = true; 
% opts.SYM = true; 
% SAGI = ( linsolve(G, AS, opts) )';

% % another trial
% SAGI = AS'*(invChol_mex(G));


btT=(SAGI*(bnT'-A*btt'))'+btt;
StT=Stt-SAGI*AS;%SAGI*SnT*SAGI';


