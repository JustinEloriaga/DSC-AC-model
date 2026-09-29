function values = dsc_rt_loglik(particles,Y,times)
% Each row of Y corresponds to one entry in times; no other data are read.
if size(Y,2)~=2||size(Y,1)~=numel(times)||any(~isfinite(Y(:)))
    error('rt:Data','Supply one finite bivariate observation per requested date.');
end
values=zeros(size(particles.B,1),numel(times));
for k=1:numel(times)
    t=times(k); h=particles.h(:,:,t); r=particles.r(:,t);
    z=(Y(k,:)-particles.B(:,:,t)).*exp(-h/2);
    rho=tanh(r);
    log_one_minus_rho2=2*(log(2)-abs(r)-log1p(exp(-2*abs(r))));
    quadratic=z(:,1).^2+(z(:,2)-rho.*z(:,1)).^2.*exp(-log_one_minus_rho2);
    values(:,k)=-log(2*pi)-sum(h,2)/2-log_one_minus_rho2/2-quadratic/2;
end
if any(isnan(values(:)))||any(values(:)==Inf)
    error('rt:Likelihood','Undefined bivariate likelihood.');
end
end
