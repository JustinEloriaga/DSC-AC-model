function tests=test_mean_kernel
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'toolbox'));
assumeTrue(testCase,exist('dsc_mean_mex','file')==3,'Build the native kernels first.');
end

function testMaskedNativeDrawMatchesReference(testCase)
rng(840,'twister'); T=21;
for m=[2 4 14]
    Y=randn(T,m); h=randn(T,m); P=zeros(m,m,T);
    for t=1:T
        X=randn(m); S=X*X'+eye(m); d=sqrt(diag(S)); P(:,:,t)=S./(d*d');
    end
    mask=true(T,m); mask(4,1)=false; mask(7,:)=false; mask(12,2:end)=false;
    Y(~mask)=NaN; h(~mask)=NaN;
    X=randn(m); V=.001*(X*X'+eye(m)); mu=randn(m,1); S=eye(m);
    noise=randn(m,T); rng_before=rng;
    actual=dsc_mean_mex(Y,h,P,mask,V,mu,S,noise);
    expected=reference_draw(Y,h,P,mask,V,mu,S,noise);
    verifyEqual(testCase,actual,expected,'AbsTol',2e-11);
    verifyEqual(testCase,rng,rng_before);
end
end

function testOneUnobservedDateUsesInitialPrior(testCase)
mu=[1;2]; S=[2 .5;.5 1]; noise=[.3;-.7];
actual=dsc_mean_mex(nan(1,2),nan(1,2),eye(2),false(1,2),100*eye(2),mu,S,noise);
verifyEqual(testCase,actual,(mu+chol(S,'lower')*noise)','AbsTol',1e-13);
end

function testInvalidDimensionsAndObservedValues(testCase)
Y=zeros(3,2); h=Y; P=repmat(eye(2),1,1,3); mask=true(3,2);
f=@(noise) dsc_mean_mex(Y,h,P,mask,eye(2),zeros(2,1),eye(2),noise);
verifyError(testCase,@()f(zeros(3,2)),'dsc:MeanInput');
bad=Y; bad(1)=NaN;
verifyError(testCase,@()dsc_mean_mex(bad,h,P,mask,eye(2),zeros(2,1),eye(2),zeros(2,3)), ...
    'dsc:MeanInput');
end

function B=reference_draw(Y,h,P,mask,V,mu,S,noise)
[T,m]=size(Y); means=zeros(T,m); covariances=zeros(m,m,T);
for t=1:T
    if t==1, predicted=S; else, predicted=S+V; end
    ix=find(mask(t,:));
    if isempty(ix)
        S=predicted;
    else
        sd=exp(h(t,ix)/2); R=(sd'*sd).*P(ix,ix,t);
        gain=predicted(:,ix)/(predicted(ix,ix)+R);
        mu=mu+gain*(Y(t,ix)'-mu(ix));
        S=predicted-gain*predicted(ix,:); S=(S+S')/2;
    end
    means(t,:)=mu'; covariances(:,:,t)=S;
end
B=zeros(T,m); B(T,:)=(mu+chol(S,'lower')*noise(:,1))';
for t=T-1:-1:1
    S=covariances(:,:,t); gain=S/(S+V);
    mu=means(t,:)'+gain*(B(t+1,:)'-means(t,:)');
    C=S-gain*S; C=(C+C')/2;
    B(t,:)=(mu+chol(C,'lower')*noise(:,T-t+1))';
end
end
