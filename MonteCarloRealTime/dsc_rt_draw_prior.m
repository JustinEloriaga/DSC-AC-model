function particles = dsc_rt_draw_prior(N,T,prior,seed)
% Draw the same DSC prior as the sampler, specialized to two series.
rng(seed,'twister');
particles.V=zeros(2,2,N);
for s=1:N
    particles.V(:,:,s)=iwishrnd(prior.Psi,prior.nu);
end
particles.sig2h=1./gamrnd(prior.a,1/prior.bh,N,2);
particles.sig2r=1./gamrnd(prior.a,1/prior.br,N,1);
particles.B=randn(N,2)*chol(prior.CB)+prior.mB;
particles.h=prior.mh+randn(N,2).*sqrt((prior.ch+1)*particles.sig2h);
particles.r=prior.mr+randn(N,1).*sqrt((prior.cr+1)*particles.sig2r);
for t=2:T
    particles=dsc_rt_append(particles);
end
end
