function tests=test_monte_carlo
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
folder=fileparts(mfilename('fullpath')); addpath(folder);
destination=tempname;
testCase.TestData.result=run_monte_carlo(1500,2,destination);
testCase.TestData.old=load(fullfile(destination,'posterior_T3.mat'),'previous');
testCase.TestData.new=load(fullfile(destination,'posterior_T4.mat'),'batch','sequential');
end

function testBatchSequentialEquality(testCase)
data=testCase.TestData.new;
verifyEqual(testCase,data.batch.particles,data.sequential.particles);
verifyEqual(testCase,data.batch.hyper,data.sequential.hyper,'AbsTol',1e-13);
verifyEqual(testCase,data.batch.weights,data.sequential.weights,'AbsTol',1e-13);
verifyEqual(testCase,data.batch.log_evidence,data.sequential.log_evidence,'AbsTol',1e-13);
verifyEqual(testCase,sum(data.batch.weights),1,'AbsTol',1e-13);
end

function testIndependentHyperparameterFormula(testCase)
data=testCase.TestData.new.batch; p=testCase.TestData.result.prior;
for s=[1 17 300]
    B=squeeze(data.particles.B(s,:,:))';
    h=squeeze(data.particles.h(s,:,:))'; r=data.particles.r(s,:)';
    expectedPsi=p.Psi+diff(B)'*diff(B);
    expectedH=p.bh+((h(1,:)-p.mh).^2/(p.ch+1)+sum(diff(h).^2,1))/2;
    expectedR=p.br+((r(1)-p.mr)^2/(p.cr+1)+sum(diff(r).^2))/2;
    actual=data.hyper(s,:);
    verifyEqual(testCase,actual,[p.nu+3 expectedPsi(1,1) expectedPsi(1,2) expectedPsi(2,2) ...
        p.a+2 expectedH p.a+2 expectedR],'AbsTol',1e-13);
end
end

function testLikelihoodMatchesMultivariateNormal(testCase)
data=testCase.TestData.new.batch; Y=testCase.TestData.result.observations;
actual=dsc_rt_loglik(data.particles,Y,1:4);
for s=[1 17 300]
    for t=1:4
        sd=exp(data.particles.h(s,:,t)/2); rho=tanh(data.particles.r(s,t));
        sigma=(sd'*sd).*[1 rho;rho 1];
        expected=log(mvnpdf(Y(t,:),data.particles.B(s,:,t),sigma));
        verifyEqual(testCase,actual(s,t),expected,'AbsTol',1e-10);
    end
end
end

function testUpdateOnlyNeedsNewDataAndUsesIt(testCase)
old=testCase.TestData.old.previous; Y=testCase.TestData.result.observations;
original_rng=rng;
normal=dsc_rt_update(old,Y(4,:));
verifyEqual(testCase,rng,original_rng);
changed=dsc_rt_update(old,Y(4,:)+[2 -2]);
verifyEqual(testCase,normal.particles,changed.particles);
verifyGreaterThan(testCase,max(abs(normal.weights-changed.weights)),1e-5);
verifyEqual(testCase,normal.weights,testCase.TestData.new.batch.weights,'AbsTol',1e-13);
end

function testOldDataAreNotCountedTwice(testCase)
old=testCase.TestData.old.previous; latest=testCase.TestData.new;
Y=testCase.TestData.result.observations;
new_ll=dsc_rt_loglik(latest.sequential.particles,Y(4,:),4);
expected=dsc_rt_normalize(old.log_weights+new_ll);
wrong=dsc_rt_normalize(2*old.log_weights+new_ll);
verifyEqual(testCase,latest.sequential.weights,expected,'AbsTol',1e-13);
verifyGreaterThan(testCase,max(abs(wrong-expected)),1e-5);
end

function testInitialStateFactorsAreIncluded(testCase)
batch=testCase.TestData.new.batch; old=testCase.TestData.old.previous;
verifyEqual(testCase,old.hyper(:,[1 5 8]),repmat([10 7.5 7.5],size(old.hyper,1),1));
verifyEqual(testCase,batch.hyper(:,[1 5 8]),repmat([11 8 8],size(batch.hyper,1),1));
end

function testOldStatesAreSmoothed(testCase)
summary=testCase.TestData.result.state_table;
verifyEqual(testCase,summary.Mean_batch,summary.Mean_sequential,'AbsTol',1e-12);
verifyEqual(testCase,summary.SD_batch,summary.SD_sequential,'AbsTol',1e-12);
verifyGreaterThan(testCase,max(abs(summary.Mean_before_y4(1:6)-summary.Mean_batch(1:6))),1e-4);
end

function testAppendPreservesStaticDrawsAndHistory(testCase)
old=testCase.TestData.old.previous.particles;
new=testCase.TestData.new.sequential.particles;
verifyEqual(testCase,new.B(:,:,1:3),old.B);
verifyEqual(testCase,new.h(:,:,1:3),old.h);
verifyEqual(testCase,new.r(:,1:3),old.r);
verifyEqual(testCase,new.V,old.V);
verifyEqual(testCase,new.sig2h,old.sig2h);
verifyEqual(testCase,new.sig2r,old.sig2r);
end
