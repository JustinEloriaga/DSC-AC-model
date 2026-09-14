function tests=test_sampler
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'core'),fullfile(root,'toolbox'));
rng(100,'twister'); T=32; m=3;
panel=struct(); panel.returns=randn(T,m)*[1 .2 -.1;0 .8 .2;0 0 .5];
panel.dates=(datetime(2003,8,15)+calweeks(0:T-1))';
panel.tickers=["A","B","C"]; panel.observation_mask=true(T,m);
[panel.pair_i,panel.pair_j]=find(tril(true(m),-1)); panel.source_hash='synthetic-test';
cfg=struct('seed',17,'p',0,'prior_weeks',16,'shrinkage',.05, ...
    'calibration_windows',[12 16 20],'calibration_window',16,'ig_shape',10, ...
    'kB',.01,'max_iterations',4,'max_seconds',90,'burnin',1,'thin',1, ...
    'chunk_size',2,'resume',false,'output_root',tempdir,'correlation_backend','auto', ...
    'correlation_threads',4);
testCase.TestData.panel=panel; testCase.TestData.cfg=cfg;
testCase.TestData.priors=dsc_calibrate_priors(panel,cfg);
end

function testBrownianExactFactor(testCase)
rng(41); z=randn(60,1); c=10; sd=.3;
C=c+min((1:60)',1:60);
verifyEqual(testCase,dsc_random_walk_draw(z,c,sd),sd*chol(C,'lower')*z,'AbsTol',1e-12);
end

function testCorrelationRoundTrip(testCase)
X=[1 .4 -.2;.4 1 .1;-.2 .1 1]; A=logm(X);
[C,stats]=dsc_corr_matrix(A(tril(true(3),-1)));
verifyEqual(testCase,C,X,'AbsTol',1e-7);
verifyTrue(testCase,stats.converged);
verifyGreaterThan(testCase,min(eig(C)),0);
end

function testPriorsAndBurninStorage(testCase)
cfg=testCase.TestData.cfg; folder=tempname; mkdir(folder);
result=dsc_sample(testCase.TestData.panel,testCase.TestData.priors,cfg,folder);
verifyEqual(testCase,result.completed_iterations,4);
verifyEqual(testCase,result.saved_draws,3);
verifyFalse(testCase,result.convergence_established);
verifyEqual(testCase,result.T,32);
files=dir(fullfile(folder,'posterior_chunk_*.mat')); iterations=[];
for i=1:numel(files)
    tmp=load(fullfile(files(i).folder,files(i).name),'chunk');
    verifySize(testCase,tmp.chunk.P_pairs,[32 3 numel(tmp.chunk.iterations)]);
    verifyTrue(testCase,all(abs(tmp.chunk.P_pairs(:))<1));
    iterations=[iterations,tmp.chunk.iterations]; %#ok<AGROW>
end
verifyEqual(testCase,iterations,2:4);
end

function testResumeMatchesUninterrupted(testCase)
cfg=testCase.TestData.cfg; panel=testCase.TestData.panel; priors=testCase.TestData.priors;
full=tempname; split=tempname; mkdir(full); mkdir(split);
dsc_sample(panel,priors,cfg,full);
cfg.max_iterations=2; dsc_sample(panel,priors,cfg,split);
cfg.max_iterations=4; cfg.resume=true; dsc_sample(panel,priors,cfg,split);
a=load(fullfile(full,'checkpoint.mat'),'checkpoint');
b=load(fullfile(split,'checkpoint.mat'),'checkpoint');
verifyEqual(testCase,a.checkpoint.state,b.checkpoint.state);
verifyEqual(testCase,a.checkpoint.rng,b.checkpoint.rng);
verifyEqual(testCase,a.checkpoint.diagnostics.sig2h,b.checkpoint.diagnostics.sig2h);
verifyEqual(testCase,a.checkpoint.diagnostics.sig2r,b.checkpoint.diagnostics.sig2r);
end

function testMaskedValuesAreIgnored(testCase)
cfg=testCase.TestData.cfg; cfg.max_iterations=2;
a=testCase.TestData.panel; a.observation_mask(20:21,2)=false;
b=a; a.returns(20:21,2)=NaN; b.returns(20:21,2)=1e12;
priorA=dsc_calibrate_priors(a,cfg); priorB=dsc_calibrate_priors(b,cfg);
verifyEqual(testCase,priorA,priorB);
fa=tempname; fb=tempname; mkdir(fa); mkdir(fb);
dsc_sample(a,priorA,cfg,fa); dsc_sample(b,priorB,cfg,fb);
aa=load(fullfile(fa,'checkpoint.mat'),'checkpoint'); bb=load(fullfile(fb,'checkpoint.mat'),'checkpoint');
verifyEqual(testCase,aa.checkpoint.state,bb.checkpoint.state);
end

function testTimeLimitPreservesCompletedState(testCase)
cfg=testCase.TestData.cfg; cfg.max_seconds=1e-6;
folder=tempname; mkdir(folder);
result=dsc_sample(testCase.TestData.panel,testCase.TestData.priors,cfg,folder);
verifyEqual(testCase,result.status,'time_limit');
verifyEqual(testCase,result.completed_iterations,0);
verifyEqual(testCase,result.saved_draws,0);
tmp=load(fullfile(folder,'checkpoint.mat'),'checkpoint');
verifyEqual(testCase,tmp.checkpoint.completed,0);
verifyEqual(testCase,tmp.checkpoint.state.B,repmat(testCase.TestData.priors.Bbar',32,1));
end

function testInitialMeanPriorIndependentOfEvolutionVariance(testCase)
% With one entirely unobserved date, B1 is a direct prior draw and must not
% depend on V. This catches an erroneous S+V prediction at the first date.
cfg=testCase.TestData.cfg; cfg.max_iterations=1; cfg.burnin=0;
panel=testCase.TestData.panel; panel.dates=panel.dates(1);
panel.returns=nan(1,3); panel.observation_mask=false(1,3);
a=testCase.TestData.priors; b=a; b.V0B=a.V0B*1000;
fa=tempname; fb=tempname; mkdir(fa); mkdir(fb);
dsc_sample(panel,a,cfg,fa); dsc_sample(panel,b,cfg,fb);
aa=load(fullfile(fa,'checkpoint.mat'),'checkpoint'); bb=load(fullfile(fb,'checkpoint.mat'),'checkpoint');
verifyEqual(testCase,aa.checkpoint.completed,1);
verifyEqual(testCase,bb.checkpoint.completed,1);
verifyEqual(testCase,aa.checkpoint.state.B,bb.checkpoint.state.B);
end

function testPerChainBurninAndStorageIdentity(testCase)
% Warm-up is removed independently and every stored chunk keeps its chain.
cfg=testCase.TestData.cfg; cfg.burnin=2; cfg.max_iterations=4;
panel=testCase.TestData.panel; priors=testCase.TestData.priors;
last_states=cell(1,2);
for chain=1:2
    cfg.chain_id=chain; folder=tempname; mkdir(folder);
    result=dsc_sample(panel,priors,cfg,folder);
    verifyEqual(testCase,result.chain_id,chain);
    verifyEqual(testCase,result.seed,cfg.seed+chain-1);
    verifyEqual(testCase,result.saved_draws,2);
    files=dir(fullfile(folder,'posterior_chunk_*.mat')); iterations=[];
    for i=1:numel(files)
        tmp=load(fullfile(files(i).folder,files(i).name),'chunk');
        verifyEqual(testCase,tmp.chunk.chain_id,chain);
        verifyEqual(testCase,tmp.chunk.dates,panel.dates);
        verifyEqual(testCase,tmp.chunk.pair_i,panel.pair_i);
        verifyEqual(testCase,tmp.chunk.pair_j,panel.pair_j);
        verifyGreaterThan(testCase,tmp.chunk.iterations,cfg.burnin);
        verifyEqual(testCase,size(tmp.chunk.P_pairs,3),numel(tmp.chunk.iterations));
        iterations=[iterations,tmp.chunk.iterations]; %#ok<AGROW>
    end
    verifyEqual(testCase,iterations,3:4);
    stored=load(fullfile(folder,'checkpoint.mat'),'checkpoint');
    verifyEqual(testCase,stored.checkpoint.identity.cfg.chain_id,chain);
    verifyEqual(testCase,stored.checkpoint.saved,2);
    last_states{chain}=stored.checkpoint.state;
end
verifyFalse(testCase,isequaln(last_states{1},last_states{2}));
end

function testResumeRejectsChangedIdentity(testCase)
cfg=testCase.TestData.cfg; cfg.max_iterations=2;
panel=testCase.TestData.panel; priors=testCase.TestData.priors;
folder=tempname; mkdir(folder);
dsc_sample(panel,priors,cfg,folder);
before=load(fullfile(folder,'checkpoint.mat'),'checkpoint');
cfg.resume=true; cfg.max_iterations=4;
changed_panel=panel; changed_panel.observation_mask(20,2)=false;
verifyError(testCase,@() dsc_sample(changed_panel,priors,cfg,folder),'dsc:ResumeMismatch');
changed_cfg=cfg; changed_cfg.seed=cfg.seed+1;
verifyError(testCase,@() dsc_sample(panel,priors,changed_cfg,folder),'dsc:ResumeMismatch');
after=load(fullfile(folder,'checkpoint.mat'),'checkpoint');
verifyEqual(testCase,after.checkpoint,before.checkpoint);
end

function testNativeAndMatlabSamplersAgree(testCase)
assumeTrue(testCase,exist('dsc_corr_batch_mex','file')==3,'Build the optional native kernel first.');
cfg=testCase.TestData.cfg; cfg.max_iterations=2;
panel=testCase.TestData.panel; panel.observation_mask(20:21,2)=false;
panel.returns(20:21,2)=NaN;
priors=dsc_calibrate_priors(panel,cfg);
fa=tempname; fb=tempname; mkdir(fa); mkdir(fb);
cfg.correlation_backend='matlab'; dsc_sample(panel,priors,cfg,fa);
cfg.correlation_backend='mex'; dsc_sample(panel,priors,cfg,fb);
a=load(fullfile(fa,'checkpoint.mat'),'checkpoint');
b=load(fullfile(fb,'checkpoint.mat'),'checkpoint');
fields=fieldnames(a.checkpoint.state);
for k=1:numel(fields)
    verifyEqual(testCase,a.checkpoint.state.(fields{k}), ...
        b.checkpoint.state.(fields{k}),'AbsTol',1e-6);
end
verifyEqual(testCase,a.checkpoint.rng,b.checkpoint.rng);
verifyEqual(testCase,a.checkpoint.diagnostics.r_proposals,b.checkpoint.diagnostics.r_proposals);
verifyEqual(testCase,a.checkpoint.diagnostics.h_proposals,b.checkpoint.diagnostics.h_proposals);
verifyEqual(testCase,a.checkpoint.diagnostics.loglik,b.checkpoint.diagnostics.loglik,'AbsTol',1e-6);
end

function testResumeRejectsDifferentImplementation(testCase)
cfg=testCase.TestData.cfg; cfg.max_iterations=1;
folder=tempname; mkdir(folder);
dsc_sample(testCase.TestData.panel,testCase.TestData.priors,cfg,folder);
loaded=load(fullfile(folder,'checkpoint.mat'),'checkpoint'); checkpoint=loaded.checkpoint;
checkpoint.identity.implementation.version=-1;
save(fullfile(folder,'checkpoint.mat'),'checkpoint','-v7.3');
cfg.resume=true; cfg.max_iterations=2;
verifyError(testCase,@()dsc_sample(testCase.TestData.panel,testCase.TestData.priors,cfg,folder), ...
    'dsc:ResumeMismatch');
end
