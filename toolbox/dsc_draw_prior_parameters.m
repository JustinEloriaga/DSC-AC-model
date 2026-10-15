function draws=dsc_draw_prior_parameters(priors,m,count,seed)
%DSC_DRAW_PRIOR_PARAMETERS Independent static parameter draws, without data.
% V ~ IW(V0B,nuB); each innovation variance ~ IG(shape,scale), with
% density proportional to x^(-shape-1)*exp(-scale/x). Caller RNG is restored.
for entry={{m,'m',2},{count,'count',1},{seed,'seed',0}}
    value=entry{1}{1}; name=entry{1}{2}; minimum=entry{1}{3};
    if ~isnumeric(value)||~isreal(value)||~isscalar(value)||~isfinite(value)|| ...
            value<minimum||value~=floor(value)
        error('dsc:PriorDrawConfig','%s must be an integer at least %d.',name,minimum);
    end
end
if seed>2^32-1
    error('dsc:PriorDrawConfig','seed must be between 0 and 2^32-1.');
end
required={'V0B','nuB','ig_shape','ig_scale_h','ig_scale_r'};
if ~isstruct(priors)||~isscalar(priors)||~all(isfield(priors,required))
    error('dsc:PriorDrawConfig','Missing inverse-Wishart or inverse-gamma prior hyperparameters.');
end
V0=priors.V0B;
if ~isnumeric(V0)||~isreal(V0)||~isequal(size(V0),[m m])|| ...
        any(~isfinite(V0(:)))||norm(V0-V0','fro')>1e-12*max(1,norm(V0,'fro'))
    error('dsc:PriorDrawConfig','V0B must be a finite symmetric m-by-m matrix.');
end
V0=(double(V0)+double(V0)')/2;
[~,bad]=chol(V0);
nu=priors.nuB;
if bad||~isnumeric(nu)||~isreal(nu)||~isscalar(nu)||~isfinite(nu)||nu<=m-1
    error('dsc:PriorDrawConfig','V0B must be positive definite and nuB must exceed m-1.');
end
q=m*(m-1)/2;
shape_h=positive_row(priors.ig_shape,m,'ig_shape');
shape_r=positive_row(priors.ig_shape,q,'ig_shape');
scale_h=positive_row(priors.ig_scale_h,m,'ig_scale_h');
scale_r=positive_row(priors.ig_scale_r,q,'ig_scale_r');
caller_rng=rng;
restore_rng=onCleanup(@()rng(caller_rng)); %#ok<NASGU>
rng(double(seed),'twister');
draws=struct('V',zeros(m,m,count),'sig2h',zeros(count,m), ...
    'sig2r',zeros(count,q),'weights',ones(1,count)/count,'draw_ids',1:count);
for d=1:count
    V=iwishrnd(V0,double(nu));
    draws.V(:,:,d)=(V+V')/2;
    draws.sig2h(d,:)=scale_h./randg(shape_h);
    draws.sig2r(d,:)=scale_r./randg(shape_r);
    [~,bad]=chol(draws.V(:,:,d));
    if bad||any(~isfinite(V(:)))|| ...
            any(~isfinite([draws.sig2h(d,:) draws.sig2r(d,:)]))|| ...
            any([draws.sig2h(d,:) draws.sig2r(d,:)]<=0)
        error('dsc:PriorDrawFailure', ...
            'Prior draw %d is not finite and positive definite; it was not redrawn.',d);
    end
end
end

function value=positive_row(value,width,name)
if ~isnumeric(value)||~isreal(value)||(~isscalar(value)&&numel(value)~=width)|| ...
        any(~isfinite(value(:)))||any(value(:)<=0)
    error('dsc:PriorDrawConfig','%s must be a finite positive scalar or %d-vector.',name,width);
end
if isscalar(value), value=repmat(double(value),1,width);
else, value=reshape(double(value),1,width); end
end
