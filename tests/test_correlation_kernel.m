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

function testNewKernelMatchesNewMatlabTransform(testCase)
assumeTrue(testCase,exist('dsc_corr_batch_mex_new','file')==3, ...
    'Compile the new correlation kernel before running this test.');
rng(2609,'twister'); T=5; m=4; [r,~]=random_correlations(T,m);
Z=randn(T,m); mask=true(T,m);
mask(2,3)=false; mask(4,[1 4])=false; Z(~mask)=NaN;
[ll,P,L,logdet,iterations]=dsc_corr_batch_mex_new(r,Z,mask,2);
[reference_ll,reference_P,reference_L,reference_logdet]= ...
    reference_batch_new(r,Z,mask);
verifyEqual(testCase,ll,reference_ll,'AbsTol',1e-7);
verifyEqual(testCase,P,reference_P,'AbsTol',1e-7);
verifyEqual(testCase,L,reference_L,'AbsTol',1e-7);
verifyEqual(testCase,logdet(:),reference_logdet,'AbsTol',1e-7);
verifyEqual(testCase,numel(iterations),T);
verifyTrue(testCase,all(iterations(:)>=1 & iterations(:)<=200));
end

function testNewKernelRoundTripsAcrossDimensions(testCase)
assumeTrue(testCase,exist('dsc_corr_batch_mex_new','file')==3, ...
    'Compile the new correlation kernel before running this test.');
rng(2716,'twister');
for m=[1 2 3 14]
    % Identity has repeated eigenvalues. Strong positive/negative matrices,
    % interleaved with random matrices, exercise abrupt warm-start changes.
    [r,expected,roundtrip_tolerance]=new_regression_correlations(m);
    T=size(r,1); Z=randn(T,m); mask=true(T,m);
    [ll,P,L,logdet,iterations]=dsc_corr_batch_mex_new(r,Z,mask,1);
    [expected_ll,expected_L,expected_logdet]=expected_batch(expected,Z,mask);
    verifyEqual(testCase,P,expected,'AbsTol',1e-7);
    verifyEqual(testCase,L,expected_L,'AbsTol',2e-7);
    verifyEqual(testCase,ll,expected_ll,'AbsTol',1e-5);
    verifyEqual(testCase,logdet(:),expected_logdet,'AbsTol',2e-6);
    verifyTrue(testCase,all(iterations(:)>=1 & iterations(:)<=200));
    sel=tril(true(m),-1);
    for t=1:T
        verifyEqual(testCase,P(:,:,t),P(:,:,t)','AbsTol',1e-13);
        verifyEqual(testCase,diag(P(:,:,t)),ones(m,1),'AbsTol',1e-13);
        verifyGreaterThan(testCase,min(eig(P(:,:,t))),0);
        verifyEqual(testCase,L(:,:,t)*L(:,:,t)',P(:,:,t),'AbsTol',1e-12);
        A=logm(P(:,:,t));
        verifyEqual(testCase,reshape(A(sel),1,[]),r(t,:),'AbsTol',roundtrip_tolerance(t));
    end
end
end

function testNewKernelStrongCoordinateWarmStarts(testCase)
assumeTrue(testCase,exist('dsc_corr_batch_mex_new','file')==3, ...
    'Compile the new correlation kernel before running this test.');
% This deterministic batch exposed full-Newton overshoot at date 95 before
% adding the initial AH phase. The old solver succeeds on all 100 dates.
rng(441,'twister'); randn(100,91); randn(100,91); %#ok<RAND>
r=randn(100,91); Z=zeros(100,14); mask=true(100,14);
[~,old_P,~,~,old_iterations]=dsc_corr_batch_mex(r,Z,mask,4);
[~,P,~,~,iterations]=dsc_corr_batch_mex_new(r,Z,mask,4);
verifyEqual(testCase,P,old_P,'AbsTol',1e-7);
verifyLessThan(testCase,max(iterations),max(old_iterations));
for t=1:100
    verifyGreaterThan(testCase,min(eig(P(:,:,t))),0);
end
end

function testNewKernelMaskCases(testCase)
assumeTrue(testCase,exist('dsc_corr_batch_mex_new','file')==3, ...
    'Compile the new correlation kernel before running this test.');
rng(8136,'twister'); T=6; m=4; [r,expected]=random_correlations(T,m);
Z=randn(T,m); mask=true(T,m);
mask(2,[1 3])=false; mask(3,:)=false;
mask(4,1:3)=false; mask(5,[2 4])=false;
Z(~mask)=NaN;
[ll,P,L,logdet]=dsc_corr_batch_mex_new(r,Z,mask,2);
[expected_ll,expected_L,expected_logdet]=expected_batch(expected,Z,mask);
verifyEqual(testCase,P,expected,'AbsTol',1e-7);
verifyEqual(testCase,ll,expected_ll,'AbsTol',1e-7);
verifyEqual(testCase,L,expected_L,'AbsTol',1e-7);
verifyEqual(testCase,logdet(:),expected_logdet,'AbsTol',1e-7);
verifyEqual(testCase,L(:,:,3),zeros(m));
verifyEqual(testCase,logdet(3),0);
verifyEqual(testCase,L(1,1,4),1);
verifyEqual(testCase,logdet(4),0);
for t=1:T
    ix=find(mask(t,:)); n=numel(ix);
    verifyEqual(testCase,L(1:n,1:n,t)*L(1:n,1:n,t)',P(ix,ix,t),'AbsTol',1e-12);
    verifyEqual(testCase,L(n+1:end,:,t),zeros(m-n,m));
    verifyEqual(testCase,L(:,n+1:end,t),zeros(m,m-n));
end
Z(~mask)=Inf;
[other_ll,other_P,other_L,other_logdet]=dsc_corr_batch_mex_new(r,Z,mask,2);
verifyEqual(testCase,other_ll,ll);
verifyEqual(testCase,other_P,P);
verifyEqual(testCase,other_L,L);
verifyEqual(testCase,other_logdet,logdet);
[unobserved_ll,unobserved_P,unobserved_L,unobserved_logdet]= ...
    dsc_corr_batch_mex_new(r,nan(T,m),false(T,m),2);
verifyEqual(testCase,unobserved_ll,0);
verifyEqual(testCase,unobserved_P,P);
verifyEqual(testCase,unobserved_L,zeros(m,m,T));
verifyEqual(testCase,unobserved_logdet,zeros(T,1));
end

function testNewKernelThreadDeterminism(testCase)
assumeTrue(testCase,exist('dsc_corr_batch_mex_new','file')==3, ...
    'Compile the new correlation kernel before running this test.');
rng(8137,'twister'); T=17; m=14; [r,expected]=random_correlations(T,m);
Z=randn(T,m); mask=true(T,m);
mask(5,[1 3 7])=false; mask(9,:)=false; mask(13,2:14)=false;
Z(~mask)=NaN;
[expected_ll,expected_L,expected_logdet]=expected_batch(expected,Z,mask);
rng_before=rng;
default=cell(1,5); [default{:}]=dsc_corr_batch_mex_new(r,Z,mask);
for threads=[1 2 4]
    first=cell(1,5); second=cell(1,5);
    [first{:}]=dsc_corr_batch_mex_new(r,Z,mask,threads);
    [second{:}]=dsc_corr_batch_mex_new(r,Z,mask,threads);
    verifyEqual(testCase,second,first);
    verifyEqual(testCase,first{1},expected_ll,'AbsTol',1e-7);
    verifyEqual(testCase,first{2},expected,'AbsTol',1e-7);
    verifyEqual(testCase,first{3},expected_L,'AbsTol',1e-7);
    verifyEqual(testCase,first{4}(:),expected_logdet,'AbsTol',1e-7);
    verifyTrue(testCase,all(first{5}(:)>=1 & first{5}(:)<=200));
    if threads==1, verifyEqual(testCase,first,default); end
end
verifyEqual(testCase,rng,rng_before);
end

function testNewKernelInvalidInputContracts(testCase)
assumeTrue(testCase,exist('dsc_corr_batch_mex_new','file')==3, ...
    'Compile the new correlation kernel before running this test.');
r=zeros(4,3); Z=zeros(4,3); mask=true(4,3);
invalid_inputs={ ...
    {zeros(4,2),Z,mask}, {r,zeros(3,3),mask}, {r,Z,true(4,2)}, ...
    {complex(r,ones(size(r))),Z,mask}, {[NaN 0 0; r(2:end,:)],Z,mask}, ...
    {r,[Inf 0 0; Z(2:end,:)],mask}, {r,[NaN 0 0; Z(2:end,:)],mask}, ...
    {r,Z,double(mask)}, {sparse(r),Z,mask}};
counts={0,-1,1.5,65,NaN,Inf,[1 2],complex(2,1),'two'};
for i=1:numel(counts), invalid_inputs{end+1}={r,Z,mask,counts{i}}; end %#ok<AGROW>
for i=1:numel(invalid_inputs)
    inputs=invalid_inputs{i};
    verifyError(testCase,@()dsc_corr_batch_mex(inputs{:}),'dsc:KernelInput');
    verifyError(testCase,@()dsc_corr_batch_mex_new(inputs{:}),'dsc:KernelInput');
end
end

function testNewKernelThreadedFailureRecovers(testCase)
assumeTrue(testCase,exist('dsc_corr_batch_mex_new','file')==3, ...
    'Compile the new correlation kernel before running this test.');
T=17; m=3; r=zeros(T,3); Z=zeros(T,m); mask=true(T,m);
bad=r; bad(10,1)=1000; rng_before=rng;
for threads=[1 2 4]
    verifyError(testCase,@()dsc_corr_batch_mex(bad,Z,mask,threads), ...
        'dsc:CorrelationNumerical');
    verifyError(testCase,@()dsc_corr_batch_mex_new(bad,Z,mask,threads), ...
        'dsc:CorrelationNumerical');
    [ll,P,L,logdet]=dsc_corr_batch_mex_new(r,Z,mask,threads);
    verifyEqual(testCase,ll,0);
    verifyEqual(testCase,P,repmat(eye(m),1,1,T));
    verifyEqual(testCase,L,P);
    verifyEqual(testCase,logdet,zeros(T,1));
end
verifyEqual(testCase,rng,rng_before);
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

function [r,C,tolerance]=new_regression_correlations(m)
if m==1
    C=ones(1,1,4); r=zeros(4,0); tolerance=1e-7*ones(4,1);
    return
end
[~,random]=random_correlations(2,m);
positive=(1-.95)*eye(m)+.95*ones(m);
rho=-.95/(m-1); negative=(1-rho)*eye(m)+rho*ones(m);
C=cat(3,eye(m),positive,random(:,:,1),negative,random(:,:,2),eye(m),positive);
T=size(C,3); r=zeros(T,m*(m-1)/2); sel=tril(true(m),-1);
for t=1:T
    A=logm(C(:,:,t)); r(t,:)=A(sel)';
end
tolerance=1e-7*ones(T,1);
% Strong correlations amplify logm's inverse sensitivity; all eigenvalues
% remain at least .05 so this still bounds a meaningful coordinate error.
tolerance([2 4 7])=5e-6;
end

function [ll,L,logdet]=expected_batch(P,Z,mask)
[T,m]=size(Z); L=zeros(m,m,T); logdet=zeros(T,1); ll=0;
for t=1:T
    ix=find(mask(t,:)); n=numel(ix);
    if n==0, continue; end
    factor=chol(P(ix,ix,t),'lower'); L(1:n,1:n,t)=factor;
    logdet(t)=2*sum(log(diag(factor))); u=factor\Z(t,ix)';
    ll=ll-.5*(logdet(t)+u'*u);
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

function [ll,P,L,logdet]=reference_batch_new(r,Z,mask)
[T,m]=size(Z); P=zeros(m,m,T); L=zeros(m,m,T);
logdet=zeros(T,1); ll=0; x=[];
for t=1:T
    [C,~,x]=dsc_corr_matrix_new(r(t,:),x); P(:,:,t)=C;
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
