function particles = dsc_rt_append(particles)
% Append one transition while keeping all old states and static parameters.
N=size(particles.B,1); T=size(particles.B,3);
v11=reshape(particles.V(1,1,:),N,1);
v12=reshape(particles.V(1,2,:),N,1);
v22=reshape(particles.V(2,2,:),N,1);
l11=sqrt(v11); l21=v12./l11; l22=sqrt(v22-l21.^2);
if any(~isfinite(l22)|l22<=0)
    error('rt:Covariance','Mean innovation covariance must be positive definite.');
end
z=randn(N,2);
particles.B(:,:,T+1)=particles.B(:,:,T)+[l11.*z(:,1),l21.*z(:,1)+l22.*z(:,2)];
particles.h(:,:,T+1)=particles.h(:,:,T)+randn(N,2).*sqrt(particles.sig2h);
particles.r(:,T+1)=particles.r(:,T)+randn(N,1).*sqrt(particles.sig2r);
end
