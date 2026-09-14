function tests=test_correlation_kernel
% Independent numerical regression against the pre-optimization algorithm.
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'toolbox'));
assumeTrue(testCase,exist('dsc_corr_batch_mex','file')==3, ...
    'Compile the optional correlation kernel before running these tests.');
end

function testRandomCorrelationRoundTrips(testCase)
rng(618,'twister');
for m=[2 3 14]
    T=7; [r,expected]=random_correlations(T,m);
    Z=randn(T,m); mask=true(T,m);
    [ll,P,L,logdet,iterations]=dsc_corr_batch_mex(r,Z,mask);
    [reference_ll,reference_P,reference_L,reference_logdet]= ...
        reference_batch(r,Z,mask);
    verifyEqual(testCase,ll,reference_ll,'AbsTol',1e-7);
    verifyEqual(testCase,P,reference_P,'AbsTol',1e-7);
    verifyEqual(testCase,P,expected,'AbsTol',2e-7);
    verifyEqual(testCase,L,reference_L,'AbsTol',1e-7);
    verifyEqual(testCase,logdet(:),reference_logdet,'AbsTol',1e-7);
    verifyEqual(testCase,numel(iterations),T);
    verifyTrue(testCase,all(iterations(:)>=1 & iterations(:)<=200));
    for t=1:T
        verifyEqual(testCase,P(:,:,t),P(:,:,t)','AbsTol',1e-13);
        verifyEqual(testCase,diag(P(:,:,t)),ones(m,1),'AbsTol',1e-13);
        verifyGreaterThan(testCase,min(eig(P(:,:,t))),0);
        verifyEqual(testCase,L(:,:,t)*L(:,:,t)',P(:,:,t),'AbsTol',1e-12);
    end
end
end

function testMaskedAndUnobservedDates(testCase)
rng(919,'twister'); T=6; m=4; [r,~]=random_correlations(T,m);
Z=randn(T,m); mask=true(T,m);
mask(2,[1 3])=false;
mask(3,:)=false;
mask(4,1:3)=false;
mask(5,[2 4])=false;
Z(~mask)=NaN;
[ll,P,L,logdet]=dsc_corr_batch_mex(r,Z,mask);
[reference_ll,reference_P,reference_L,reference_logdet]= ...
    reference_batch(r,Z,mask);
verifyEqual(testCase,ll,reference_ll,'AbsTol',1e-7);
verifyEqual(testCase,P,reference_P,'AbsTol',1e-7);
verifyEqual(testCase,L,reference_L,'AbsTol',1e-7);
verifyEqual(testCase,logdet(:),reference_logdet,'AbsTol',1e-7);
verifyEqual(testCase,L(:,:,3),zeros(m));
verifyEqual(testCase,logdet(3),0);
for t=1:T
    ix=find(mask(t,:)); n=numel(ix);
    verifyEqual(testCase,L(1:n,1:n,t)*L(1:n,1:n,t)',P(ix,ix,t),'AbsTol',1e-12);
    verifyEqual(testCase,L(n+1:end,:,t),zeros(m-n,m));
    verifyEqual(testCase,L(:,n+1:end,t),zeros(m,m-n));
end
% Arbitrary values outside the observed subspace cannot affect any output.
Z(~mask)=Inf;
[other_ll,other_P,other_L,other_logdet]=dsc_corr_batch_mex(r,Z,mask);
verifyEqual(testCase,other_ll,ll);
verifyEqual(testCase,other_P,P);
verifyEqual(testCase,other_L,L);
verifyEqual(testCase,other_logdet,logdet);
end

function testAllDatesUnobserved(testCase)
rng(241,'twister'); T=3; m=3; [r,~]=random_correlations(T,m);
[ll,P,L,logdet]=dsc_corr_batch_mex(r,nan(T,m),false(T,m));
verifyEqual(testCase,ll,0);
verifyEqual(testCase,L,zeros(m,m,T));
verifyEqual(testCase,logdet(:),zeros(T,1));
verifyTrue(testCase,all(isfinite(P(:))));
end

function testRepeatedCallsAreDeterministic(testCase)
rng(752,'twister'); T=5; m=3; [r,~]=random_correlations(T,m);
Z=randn(T,m); mask=true(T,m); first=cell(1,5); second=cell(1,5);
rng_before=rng;
[first{:}]=dsc_corr_batch_mex(r,Z,mask);
[second{:}]=dsc_corr_batch_mex(r,Z,mask);
verifyEqual(testCase,second,first);
verifyEqual(testCase,rng,rng_before);
end

function testThreadedLikelihoodsAndDeterminism(testCase)
% Seventeen dates exercise uneven contiguous blocks for 2 and 4 workers.
rng(1821,'twister'); T=17; m=14; [r,~]=random_correlations(T,m);
Z=randn(T,m); mask=true(T,m);
mask(5,[1 3 7])=false;
mask(9,:)=false;
mask(13,2:14)=false;
Z(~mask)=NaN;
[reference_ll,reference_P,reference_L,reference_logdet]= ...
    reference_batch(r,Z,mask);
rng_before=rng;
default=cell(1,5); [default{:}]=dsc_corr_batch_mex(r,Z,mask);
for threads=[1 2 4]
    first=cell(1,5); second=cell(1,5);
    [first{:}]=dsc_corr_batch_mex(r,Z,mask,threads);
    [second{:}]=dsc_corr_batch_mex(r,Z,mask,threads);
    verifyEqual(testCase,second,first);
    verifyEqual(testCase,first{1},reference_ll,'AbsTol',1e-7);
    verifyEqual(testCase,first{2},reference_P,'AbsTol',1e-7);
    verifyEqual(testCase,first{3},reference_L,'AbsTol',1e-7);
    verifyEqual(testCase,first{4}(:),reference_logdet,'AbsTol',1e-7);
    verifyEqual(testCase,numel(first{5}),T);
    verifyTrue(testCase,all(first{5}(:)>=1 & first{5}(:)<=200));
    if threads==1, verifyEqual(testCase,first,default); end
end
verifyEqual(testCase,rng,rng_before);
end

function testInvalidThreadCountsFail(testCase)
r=zeros(4,3); Z=zeros(4,3); mask=true(4,3);
counts={0,-1,1.5,65,NaN,Inf,[1 2],complex(2,1),'two'};
for i=1:numel(counts)
    failed=false;
    try
        ignored=dsc_corr_batch_mex(r,Z,mask,counts{i}); %#ok<NASGU>
    catch
        failed=true;
    end
    verifyTrue(testCase,failed,sprintf('Invalid thread-count case %d must fail.',i));
end
end

function testThreadedTransformFailureRecovers(testCase)
% A finite but extreme coordinate fails inside a worker, after launch.
% Returning the error must join workers and leave later calls usable.
T=17; m=3; r=zeros(T,3); Z=zeros(T,m); mask=true(T,m);
bad=r; bad(10,1)=1000;
rng_before=rng;
for threads=[2 4]
    verifyError(testCase,@()dsc_corr_batch_mex(bad,Z,mask,threads), ...
        'dsc:CorrelationNumerical');
    [ll,P,L,logdet]=dsc_corr_batch_mex(r,Z,mask,threads);
    verifyEqual(testCase,ll,0);
    verifyEqual(testCase,P,repmat(eye(m),1,1,T));
    verifyEqual(testCase,L,P);
    verifyEqual(testCase,logdet,zeros(T,1));
end
verifyEqual(testCase,rng,rng_before);
end

function testInvalidInputsFail(testCase)
r=zeros(4,3); Z=zeros(4,3); mask=true(4,3);
invalid={ ...
    @() dsc_corr_batch_mex(zeros(4,2),Z,mask), ... % non-triangular coordinate count
    @() dsc_corr_batch_mex(r,zeros(3,3),mask), ... % inconsistent dates
    @() dsc_corr_batch_mex(r,Z,true(4,2)), ... % inconsistent mask
    @() dsc_corr_batch_mex(complex(r,ones(size(r))),Z,mask), ...
    @() dsc_corr_batch_mex([NaN 0 0; r(2:end,:)],Z,mask), ...
    @() dsc_corr_batch_mex(r,[Inf 0 0; Z(2:end,:)],mask), ...
    @() dsc_corr_batch_mex(r,[NaN 0 0; Z(2:end,:)],mask)};
for i=1:numel(invalid)
    failed=false;
    try
        ignored=invalid{i}(); %#ok<NASGU>
    catch
        failed=true;
    end
    verifyTrue(testCase,failed,sprintf('Invalid-input case %d must fail.',i));
end
end

function [r,C]=random_correlations(T,m)
r=zeros(T,m*(m-1)/2); C=zeros(m,m,T); sel=tril(true(m),-1);
for t=1:T
    X=randn(m); S=X*X'/m+.5*eye(m); sd=sqrt(diag(S));
    C(:,:,t)=S./(sd*sd'); A=logm(C(:,:,t)); r(t,:)=A(sel)';
end
end

function [ll,P,L,logdet]=reference_batch(r,Z,mask)
[T,m]=size(Z); P=zeros(m,m,T); L=zeros(m,m,T);
logdet=zeros(T,1); ll=0; x=[];
for t=1:T
    [C,x]=frozen_original_transform(r(t,:),x); P(:,:,t)=C;
    ix=find(mask(t,:)); n=numel(ix);
    if n==0, continue; end
    factor=chol(C(ix,ix),'lower'); L(1:n,1:n,t)=factor;
    logdet(t)=2*sum(log(diag(factor))); u=factor\Z(t,ix)';
    ll=ll-.5*(logdet(t)+u'*u);
end
end

function [C,x]=frozen_original_transform(r,x)
% Frozen full-product implementation: do not call the optimized MATLAB map.
n=(1+sqrt(1+8*numel(r)))/2;
if isempty(x), x=zeros(n,1); end
sel=tril(true(n),-1); diagonal=1:n+1:n*n;
A=zeros(n); A(sel)=r; A=A+A'; tol=1e-8*sqrt(n); converged=false;
for iteration=1:200
    A(diagonal)=x;
    [Q,D]=eig(A,'vector'); E=Q*diag(exp(D))*Q';
    delta=log(real(diag(E))); x=x-delta;
    if norm(delta)<tol, converged=true; break; end
end
assert(converged,'Reference correlation transform did not converge.');
A(diagonal)=x;
[Q,D]=eig(A,'vector'); C=real(Q*diag(exp(D))*Q');
C=(C+C')/2; C(diagonal)=1;
[~,bad]=chol(C); assert(bad==0,'Reference correlation is not positive definite.');
end
