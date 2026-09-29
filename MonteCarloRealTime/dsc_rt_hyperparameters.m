function [hyper,names] = dsc_rt_hyperparameters(particles,prior)
% Conditional on each path: IW(nu,Psi), IG(a_h,b_hj), IG(a_r,b_r).
N=size(particles.B,1); T=size(particles.B,3);
names=["nu_V","Psi_V11","Psi_V12","Psi_V22","a_h","b_h1","b_h2","a_r","b_r"];
base=[prior.nu+T-1,prior.Psi(1,1),prior.Psi(1,2),prior.Psi(2,2), ...
    prior.a+T/2,prior.bh,prior.bh,prior.a+T/2,prior.br];
hyper=repmat(base,N,1);
hyper(:,6:7)=hyper(:,6:7)+(particles.h(:,:,1)-prior.mh).^2/(2*(prior.ch+1));
hyper(:,9)=hyper(:,9)+(particles.r(:,1)-prior.mr).^2/(2*(prior.cr+1));
for t=2:T
    dB=particles.B(:,:,t)-particles.B(:,:,t-1);
    hyper(:,2:4)=hyper(:,2:4)+[dB(:,1).^2,dB(:,1).*dB(:,2),dB(:,2).^2];
    hyper(:,6:7)=hyper(:,6:7)+(particles.h(:,:,t)-particles.h(:,:,t-1)).^2/2;
    hyper(:,9)=hyper(:,9)+(particles.r(:,t)-particles.r(:,t-1)).^2/2;
end
end
