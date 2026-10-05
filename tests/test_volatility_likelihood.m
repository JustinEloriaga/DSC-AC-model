function tests=test_volatility_likelihood
tests=functiontests(localfunctions);
end

function setupOnce(~)
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'toolbox'));
end

function testCoordinateChangeMatchesFullLikelihood(testCase)
rng(166,'twister'); T=27; m=4; E=randn(T,m); h=randn(T,m);
mask=true(T,m); mask(2,2:4)=false; mask(9,:)=false; mask(15,[1 3])=false;
E(~mask)=NaN; P=zeros(m,m,T); directions=zeros(m,m,T); u=zeros(m,1,T);
z=E.*exp(-h/2); z(~mask)=0;
for t=1:T
    X=randn(m); S=X*X'+eye(m); sd=sqrt(diag(S)); P(:,:,t)=S./(sd*sd');
    ix=find(mask(t,:)); n=numel(ix); if n==0, continue; end
    L=chol(P(ix,ix,t),'lower'); directions(1:n,ix,t)=L\eye(n);
    u(1:n,1,t)=L\z(t,ix)';
end
base=full_likelihood(E,h,P,mask);
for j=1:m
    direction=directions(:,j,:); cross=reshape(sum(direction.*u,1),T,1);
    precision=reshape(sum(direction.^2,1),T,1);
    for scale=[0 .2 2]
        x=h(:,j)+scale*randn(T,1); candidate=h; candidate(:,j)=x;
        actual=dsc_volatility_likelihood(x,h(:,j),E(:,j),z(:,j),cross,precision,mask(:,j),base);
        verifyEqual(testCase,actual,full_likelihood(E,candidate,P,mask),'AbsTol',2e-10);
    end
end
end

function testNoObservationsAndOverflow(testCase)
verifyEqual(testCase,dsc_volatility_likelihood(NaN,NaN,NaN,0,0,0,false,3),3);
verifyEqual(testCase,dsc_volatility_likelihood(-2000,0,1,1,1,1,true,-1),-Inf);
end

function ll=full_likelihood(E,h,P,mask)
ll=0;
for t=1:size(E,1)
    ix=find(mask(t,:)); if isempty(ix), continue; end
    L=chol(P(ix,ix,t),'lower'); z=E(t,ix)'.*exp(-h(t,ix)'/2); u=L\z;
    ll=ll-.5*(numel(ix)*log(2*pi)+sum(h(t,ix))+2*sum(log(diag(L)))+u'*u);
end
end
