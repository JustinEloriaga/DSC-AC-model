function tests=test_chain_runner
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'core'),fullfile(root,'toolbox'));
rng(219,'twister'); T=12; m=2;
panel=struct('returns',.3*randn(T,m), ...
    'dates',(datetime(2004,1,2)+calweeks(0:T-1))', ...
    'tickers',["A","B"],'observation_mask',true(T,m), ...
    'pair_i',2,'pair_j',1,'source_hash','synthetic-chain-runner');
priors=struct('Bbar',[0;0],'VBbar',.01*eye(m),'nuB',12,'V0B',.001*eye(m), ...
    'mh0',log([.09;.09]),'mr0',0,'h0_scale',10,'r0_scale',10, ...
    'sig2h_mean',.01,'sig2r_mean',.01,'ig_shape',10, ...
    'ig_scale_h',.09,'ig_scale_r',.09,'tickers',panel.tickers, ...
    'source_hash',panel.source_hash);
cfg=struct('seed',43,'chain_id',3,'p',0,'max_iterations',4,'max_seconds',90, ...
    'burnin',1,'thin',1,'chunk_size',2,'resume',false, ...
    'estimate_parameters',true,'correlation_backend','matlab','correlation_threads',1, ...
    'num_chains',2,'parallel_chains',false,'convergence_mode','report', ...
    'auto_extend',false,'check_every',1,'max_retained_draws',4, ...
    'rhat_threshold',1.01,'min_ess',400,'max_mcse_ratio',.05,'min_diagnostic_draws',100);
testCase.TestData.panel=panel; testCase.TestData.priors=priors; testCase.TestData.cfg=cfg;
end

function testBatchesMatchIndependentChains(testCase)
cfg=testCase.TestData.cfg; folder=tempname;
result=dsc_run_chains(testCase.TestData.panel,testCase.TestData.priors,cfg,folder);
verifyEqual(testCase,result.chain_ids,[3 4]);
verifyEqual(testCase,result.seeds,[45 46]);
verifyEqual(testCase,result.saved_draws_per_chain,[3 3]);
verifyEqual(testCase,result.saved_draws,6);
verifyEqual(testCase,result.parameter_draws,6);
verifyEqual(testCase,result.completed_iterations,8);
verifyEqual(testCase,result.session_completed_iterations,8);
verifyFalse(testCase,isfield(result.paths,'checkpoint'));
verifyFalse(testCase,isfile(fullfile(folder,'checkpoint.mat')));
verifyTrue(testCase,isfile(fullfile(folder,'result.mat')));
for c=1:2
    explicit=sampler_cfg(cfg); explicit.chain_id=cfg.chain_id+c-1;
    standalone=tempname;
    expected=dsc_sample(testCase.TestData.panel,testCase.TestData.priors,explicit,standalone);
    a=load(fullfile(standalone,'checkpoint.mat'),'checkpoint');
    b=load(fullfile(result.chain_dirs{c},'checkpoint.mat'),'checkpoint');
    verifyEqual(testCase,a.checkpoint.state,b.checkpoint.state);
    verifyEqual(testCase,a.checkpoint.rng,b.checkpoint.rng);
    verifyEqual(testCase,a.checkpoint.parameter_sums,b.checkpoint.parameter_sums);
    verifyEqual(testCase,expected.parameter_estimate,result.chain_results{c}.parameter_estimate);
    chunks=dir(fullfile(result.chain_dirs{c},'posterior_chunk_*.mat')); iterations=[];
    for k=1:numel(chunks)
        loaded=load(fullfile(chunks(k).folder,chunks(k).name),'chunk');
        verifyEqual(testCase,loaded.chunk.chain_id,explicit.chain_id);
        iterations=[iterations loaded.chunk.iterations]; %#ok<AGROW>
    end
    verifyEqual(testCase,iterations,2:4);
end
verifyEqual(testCase,result.parameter_estimate.V, ...
    (result.chain_results{1}.parameter_estimate.V+result.chain_results{2}.parameter_estimate.V)/2, ...
    'AbsTol',1e-14);
end

function testResumeExtensionKeepsStatesAndRng(testCase)
cfg=testCase.TestData.cfg; folder=tempname; cfg.max_iterations=2;
first=dsc_run_chains(testCase.TestData.panel,testCase.TestData.priors,cfg,folder);
verifyEqual(testCase,first.saved_draws_per_chain,[1 1]);
cfg.max_iterations=4; cfg.resume=true;
resumed=dsc_run_chains(testCase.TestData.panel,testCase.TestData.priors,cfg,folder);
verifyEqual(testCase,resumed.saved_draws_per_chain,[3 3]);
verifyEqual(testCase,resumed.session_completed_iterations,4);
verifyGreaterThan(testCase,resumed.stage_elapsed_seconds,first.stage_elapsed_seconds);
for c=1:2
    explicit=sampler_cfg(cfg); explicit.chain_id=cfg.chain_id+c-1; explicit.resume=false;
    standalone=tempname;
    dsc_sample(testCase.TestData.panel,testCase.TestData.priors,explicit,standalone);
    a=load(fullfile(standalone,'checkpoint.mat'),'checkpoint');
    b=load(fullfile(resumed.chain_dirs{c},'checkpoint.mat'),'checkpoint');
    verifyEqual(testCase,a.checkpoint.state,b.checkpoint.state);
    verifyEqual(testCase,a.checkpoint.rng,b.checkpoint.rng);
    verifyEqual(testCase,b.checkpoint.saved,3);
    verifyEqual(testCase,b.checkpoint.diagnostics.iteration,(1:4)');
end
end

function testCompletedResumeStillValidatesIdentity(testCase)
cfg=testCase.TestData.cfg; cfg.max_iterations=2; folder=tempname;
dsc_run_chains(testCase.TestData.panel,testCase.TestData.priors,cfg,folder);
cfg.resume=true; cfg.seed=cfg.seed+1;
verifyError(testCase,@() dsc_run_chains(testCase.TestData.panel,testCase.TestData.priors,cfg,folder), ...
    'dsc:ResumeMismatch');
end

function testStaleBatchSummaryRecoversCommittedCheckpoint(testCase)
cfg=testCase.TestData.cfg; cfg.max_iterations=2; folder=tempname;
first=dsc_run_chains(testCase.TestData.panel,testCase.TestData.priors,cfg,folder);
chain_dir=first.chain_dirs{1}; stale=load(fullfile(chain_dir,'result.mat'),'result');
explicit=sampler_cfg(cfg); explicit.chain_id=cfg.chain_id;
explicit.resume=true; explicit.max_iterations=3;
dsc_sample(testCase.TestData.panel,testCase.TestData.priors,explicit,chain_dir);
committed=load(fullfile(chain_dir,'checkpoint.mat'),'checkpoint');
% Simulate interruption after checkpoint commit but before batch summary.
result=stale.result; save(fullfile(chain_dir,'result.mat'),'result','-v7.3');
cfg.resume=true;
recovered=dsc_run_chains(testCase.TestData.panel,testCase.TestData.priors,cfg,folder);
verifyEqual(testCase,recovered.saved_draws_per_chain,[2 1]);
verifyEqual(testCase,recovered.completed_iterations_per_chain,[3 2]);
after=load(fullfile(chain_dir,'checkpoint.mat'),'checkpoint');
verifyEqual(testCase,after.checkpoint.state,committed.checkpoint.state);
verifyEqual(testCase,after.checkpoint.rng,committed.checkpoint.rng);
verifyEqual(testCase,recovered.parameter_estimate.V, ...
    (2*recovered.chain_results{1}.parameter_estimate.V+recovered.chain_results{2}.parameter_estimate.V)/3, ...
    'AbsTol',1e-14);
end

function testOneChainKeepsDirectLayout(testCase)
cfg=testCase.TestData.cfg; cfg.num_chains=1; cfg.max_iterations=2;
folder=tempname;
result=dsc_run_chains(testCase.TestData.panel,testCase.TestData.priors,cfg,folder);
verifyEqual(testCase,result.chain_dirs,{folder});
verifyEqual(testCase,result.paths.checkpoint,fullfile(folder,'checkpoint.mat'));
verifyTrue(testCase,isfile(result.paths.checkpoint));
verifyNotEmpty(testCase,dir(fullfile(folder,'posterior_chunk_*.mat')));
verifyFalse(testCase,isfolder(fullfile(folder,'chains')));
verifyFalse(testCase,result.randomize_initial_state);
verifyFalse(testCase,result.convergence_established);
cfg.resume=true; cfg.max_iterations=3;
resumed=dsc_run_chains(testCase.TestData.panel,testCase.TestData.priors,cfg,folder);
verifyEqual(testCase,resumed.saved_draws,2);
end

function testAutomaticExtensionStopsAtPerChainCap(testCase)
cfg=testCase.TestData.cfg; cfg.max_iterations=2; cfg.auto_extend=true;
cfg.max_retained_draws=3; cfg.check_every=2;
result=dsc_run_chains(testCase.TestData.panel,testCase.TestData.priors,cfg,tempname);
verifyEqual(testCase,result.saved_draws_per_chain,[3 3]);
verifyEqual(testCase,result.status,'retained_draw_limit');
verifyEqual(testCase,result.requested_retained_draws,1);
verifyFalse(testCase,result.convergence_established);
for c=1:2
    loaded=load(fullfile(result.chain_dirs{c},'checkpoint.mat'),'checkpoint');
    verifyEqual(testCase,loaded.checkpoint.diagnostics.iteration,(1:4)');
end
end

function testCumulativeTimeBudgetDoesNotRestart(testCase)
cfg=testCase.TestData.cfg; cfg.max_seconds=1e-6; cfg.auto_extend=true;
folder=tempname;
first=dsc_run_chains(testCase.TestData.panel,testCase.TestData.priors,cfg,folder);
verifyEqual(testCase,first.status,'time_limit');
verifyEqual(testCase,first.saved_draws_per_chain,[0 0]);
cfg.resume=true;
again=dsc_run_chains(testCase.TestData.panel,testCase.TestData.priors,cfg,folder);
verifyEqual(testCase,again.status,'time_limit');
verifyEqual(testCase,again.completed_iterations,0);
verifyGreaterThan(testCase,again.stage_elapsed_seconds,first.stage_elapsed_seconds);
verifyEmpty(testCase,dir(fullfile(folder,'chains','*','posterior_chunk_*.mat')));
end

function testDispersedStartsAreSeededValidAndPartOfResumeIdentity(testCase)
cfg=sampler_cfg(testCase.TestData.cfg); cfg.max_seconds=1e-6;
panel=testCase.TestData.panel; priors=testCase.TestData.priors;
first=tempname; second=tempname; repeated=tempname;
dsc_sample(panel,priors,cfg,first); dsc_sample(panel,priors,cfg,repeated);
cfg.chain_id=cfg.chain_id+1; dsc_sample(panel,priors,cfg,second);
a=load(fullfile(first,'checkpoint.mat'),'checkpoint');
b=load(fullfile(second,'checkpoint.mat'),'checkpoint');
c=load(fullfile(repeated,'checkpoint.mat'),'checkpoint');
verifyEqual(testCase,a.checkpoint.completed,0);
verifyEqual(testCase,a.checkpoint.state,c.checkpoint.state);
verifyEqual(testCase,a.checkpoint.rng,c.checkpoint.rng);
for name={'B','h','r','V','sig2h','sig2r'}
    verifyFalse(testCase,isequal(a.checkpoint.state.(name{1}),b.checkpoint.state.(name{1})));
    verifyTrue(testCase,all(isfinite(a.checkpoint.state.(name{1})),'all'));
end
verifyGreaterThan(testCase,min(eig(a.checkpoint.state.V)),0);
verifyGreaterThan(testCase,a.checkpoint.state.sig2h,0);
verifyGreaterThan(testCase,a.checkpoint.state.sig2r,0);
cfg.chain_id=cfg.chain_id-1; cfg.resume=true; cfg.randomize_initial_state=false;
verifyError(testCase,@() dsc_sample(panel,priors,cfg,first),'dsc:ResumeMismatch');
end

function testUnsupportedConvergenceModeAndInvalidCapFailBeforeWrites(testCase)
cfg=testCase.TestData.cfg; cfg.convergence_mode='required'; folder=tempname;
verifyError(testCase,@() dsc_run_chains(testCase.TestData.panel,testCase.TestData.priors,cfg,folder), ...
    'MATLAB:unrecognizedStringChoice');
verifyFalse(testCase,isfolder(folder));
cfg=testCase.TestData.cfg; cfg.auto_extend=true; cfg.max_retained_draws=2;
verifyError(testCase,@() dsc_run_chains(testCase.TestData.panel,testCase.TestData.priors,cfg,folder), ...
    'dsc:ChainConfig');
verifyFalse(testCase,isfolder(folder));
end

function testMissingToolboxFailsBeforeWrites(testCase)
assumeFalse(testCase,dsc_parallel_available(),'Toolbox is installed; missing-toolbox branch is unavailable.');
cfg=testCase.TestData.cfg; cfg.parallel_chains=true; folder=tempname;
verifyError(testCase,@() dsc_run_chains(testCase.TestData.panel,testCase.TestData.priors,cfg,folder), ...
    'dsc:ParallelUnavailable');
verifyFalse(testCase,isfolder(folder));
end

function cfg=sampler_cfg(cfg)
names={'num_chains','parallel_chains','convergence_mode','auto_extend', ...
    'check_every','max_retained_draws','rhat_threshold','min_ess','max_mcse_ratio','min_diagnostic_draws'};
cfg=rmfield(cfg,names); cfg.randomize_initial_state=true;
end
