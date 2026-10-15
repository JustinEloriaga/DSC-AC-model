function tests=test_prior_smoothing
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'toolbox'));
saved_rng=rng; restore=onCleanup(@()rng(saved_rng)); %#ok<NASGU>
rng(219,'twister'); T=10; m=2;
panel=struct('returns',.3*randn(T,m), ...
    'dates',(datetime(2004,1,2)+calweeks(0:T-1))', ...
    'tickers',["A","B"],'observation_mask',true(T,m), ...
    'pair_i',2,'pair_j',1,'source_hash','synthetic-prior-smoothing');
priors=struct('Bbar',[0;0],'VBbar',.01*eye(m),'nuB',12,'V0B',.001*eye(m), ...
    'mh0',log([.09;.09]),'mr0',0,'h0_scale',10,'r0_scale',10, ...
    'sig2h_mean',.01,'sig2r_mean',.01,'ig_shape',10, ...
    'ig_scale_h',.09,'ig_scale_r',.09,'tickers',panel.tickers, ...
    'source_hash',panel.source_hash);
cfg=struct('seed',43,'chain_id',3,'p',0,'max_iterations',3,'max_seconds',90, ...
    'burnin',1,'thin',1,'chunk_size',2,'resume',false, ...
    'estimate_parameters',false,'parameter_source','prior','prior_parameter_draws',2, ...
    'correlation_backend','matlab','correlation_threads',1, ...
    'num_chains',2,'parallel_chains',false,'convergence_mode','report', ...
    'auto_extend',false,'check_every',2,'max_retained_draws',4, ...
    'rhat_threshold',1.01,'min_ess',400,'max_mcse_ratio',.05,'min_diagnostic_draws',6);
testCase.TestData.panel=panel; testCase.TestData.priors=priors; testCase.TestData.cfg=cfg;
end

function testSeedRngAndPositiveDefiniteDraws(testCase)
before=rng;
a=dsc_draw_prior_parameters(testCase.TestData.priors,2,8,217);
verifyEqual(testCase,rng,before);
b=dsc_draw_prior_parameters(testCase.TestData.priors,2,8,217);
c=dsc_draw_prior_parameters(testCase.TestData.priors,2,8,218);
verifyEqual(testCase,a,b); verifyNotEqual(testCase,a.V,c.V);
verifySize(testCase,a.V,[2 2 8]);
verifySize(testCase,a.sig2h,[8 2]); verifySize(testCase,a.sig2r,[8 1]);
verifyEqual(testCase,a.weights,ones(1,8)/8); verifyEqual(testCase,a.draw_ids,1:8);
verifyTrue(testCase,all(a.sig2h(:)>0)); verifyTrue(testCase,all(a.sig2r(:)>0));
for d=1:8, verifyGreaterThan(testCase,min(eig(a.V(:,:,d))),0); end
end

function testPriorMeansUseIgScaleAndIwScale(testCase)
priors=testCase.TestData.priors;
draws=dsc_draw_prior_parameters(priors,2,12000,20261015);
% Broad Monte Carlo tolerances verify scale conventions, not exact sequences.
verifyEqual(testCase,mean(draws.sig2h,1),[.01 .01],'RelTol',.025);
verifyEqual(testCase,mean(draws.sig2r,1),.01,'RelTol',.025);
expected=priors.V0B/(priors.nuB-2-1);
verifyEqual(testCase,diag(mean(draws.V,3)),diag(expected),'RelTol',.035);
end

function testInvalidPriorFailsBeforeTouchingRng(testCase)
priors=testCase.TestData.priors; before=rng;
priors.V0B=[1 2;2 1];
verifyError(testCase,@()dsc_draw_prior_parameters(priors,2,1,0),'dsc:PriorDrawConfig');
verifyEqual(testCase,rng,before);
priors=testCase.TestData.priors; priors.nuB=1;
verifyError(testCase,@()dsc_draw_prior_parameters(priors,2,1,0),'dsc:PriorDrawConfig');
priors=testCase.TestData.priors; priors.ig_scale_r=0;
verifyError(testCase,@()dsc_draw_prior_parameters(priors,2,1,0),'dsc:PriorDrawConfig');
verifyError(testCase,@()dsc_draw_prior_parameters(testCase.TestData.priors,2,0,0),'dsc:PriorDrawConfig');
verifyError(testCase,@()dsc_draw_prior_parameters(testCase.TestData.priors,2,1,2^32),'dsc:PriorDrawConfig');
verifyEqual(testCase,rng,before);
end

function testFixedThetaPerComponentWarmupIdentityAndMixture(testCase)
cfg=testCase.TestData.cfg; panel=testCase.TestData.panel; priors=testCase.TestData.priors;
folder=new_run_folder(testCase);
result=dsc_run_prior_smoothing(panel,priors,cfg,folder);
verifyTrue(testCase,result.prior_sampling.mixture_complete);
verifyEqual(testCase,result.prior_sampling.completed_parameter_draws,2);
verifyEqual(testCase,result.chain_ids,3:6); verifyEqual(testCase,result.seeds,45:48);
verifyEqual(testCase,result.saved_draws_per_chain,[2 2 2 2]);
verifyEqual(testCase,result.completed_iterations_per_chain,[3 3 3 3]);
verifyEqual(testCase,result.saved_draws,8);
verifyEqual(testCase,result.prior_sampling.parameter_draw_ids_per_chain,[1 1 2 2]);
verifyEqual(testCase,result.prior_sampling.mixture_weights,[.5 .5]);
verifyEqual(testCase,result.inference_mode,'prior_parameter_smoothing');
verifyEqual(testCase,result.parameter_smoothing,'prior_draws');
verifyFalse(testCase,result.estimate_parameters); verifyEmpty(testCase,result.parameter_estimate);
verifyEqual(testCase,result.convergence.status,'not_checked');
verifyFalse(testCase,result.convergence_established);
verifyFalse(testCase,result.prior_sampling.all_conditional_convergence_passed);
verifyNumElements(testCase,result.prior_sampling.conditional_convergence,2);
verifyEqual(testCase,result.stage_seconds.mean_variance,0);
verifyEqual(testCase,result.stage_seconds.correlation_variance,0);
verifyEqual(testCase,result.stage_seconds.volatility_variance,0);
loaded=load(result.paths.prior_parameters,'prior_parameters'); bundle=loaded.prior_parameters;
expected=dsc_draw_prior_parameters(priors,2,2,mod(cfg.seed+104729,2^32));
verifyEqual(testCase,bundle.prior_draws,expected);
verifyEqual(testCase,bundle.parameter_source,'prior');
verifyFalse(testCase,isfile(fullfile(folder,'latest_parameters.mat')));
verifyTrue(testCase,isfile(result.paths.prior_sampling_metadata));
verifyTrue(testCase,isfile(result.paths.result));
for k=1:4
    d=ceil(k/2); loaded=load(fullfile(result.chain_dirs{k},'checkpoint.mat'),'checkpoint');
    cp=loaded.checkpoint;
    verifyFalse(testCase,cp.identity.cfg.estimate_parameters);
    verifyEqual(testCase,cp.state.V,expected.V(:,:,d));
    verifyEqual(testCase,cp.state.sig2h,expected.sig2h(d,:));
    verifyEqual(testCase,cp.state.sig2r,expected.sig2r(d,:));
    verifyEqual(testCase,cp.identity.cfg.chain_id,2+k);
    files=dir(fullfile(result.chain_dirs{k},'posterior_chunk_*.mat')); iterations=[];
    for f=1:numel(files)
        loaded=load(fullfile(files(f).folder,files(f).name),'chunk'); chunk=loaded.chunk;
        n=numel(chunk.iterations); iterations=[iterations chunk.iterations]; %#ok<AGROW>
        verifyEqual(testCase,chunk.V,repmat(expected.V(:,:,d),1,1,n));
        verifyEqual(testCase,chunk.sig2h,repmat(expected.sig2h(d,:),n,1));
        verifyEqual(testCase,chunk.sig2r,repmat(expected.sig2r(d,:),n,1));
    end
    verifyEqual(testCase,iterations,2:3);
end
component=load(fullfile(folder,'prior_component_001','result.mat'),'result');
verifyLessThanOrEqual(testCase,component.result.stage_elapsed_seconds,result.total_seconds);
verifyEqual(testCase,result.prior_sampling.elapsed_component_wall_seconds, ...
    sum(result.prior_sampling.conditional_component_elapsed_seconds));
end

function testSharedTimeCapPreservesAllPriorDrawsAndIncompleteMixture(testCase)
cfg=testCase.TestData.cfg; cfg.max_seconds=1e-9; folder=new_run_folder(testCase);
result=dsc_run_prior_smoothing(testCase.TestData.panel,testCase.TestData.priors,cfg,folder);
verifyEqual(testCase,result.status,'time_limit');
verifyFalse(testCase,result.prior_sampling.mixture_complete);
verifyEqual(testCase,result.prior_sampling.completed_parameter_draws,0);
verifyEqual(testCase,result.saved_draws_per_chain,[0 0 0 0]);
verifyTrue(testCase,isfile(result.paths.prior_parameters));
loaded=load(result.paths.prior_parameters,'prior_parameters');
verifyEqual(testCase,size(loaded.prior_parameters.prior_draws.V,3),2);
verifyTrue(testCase,isempty(dir(fullfile(folder,'prior_component_*'))));
verifyEqual(testCase,result.prior_sampling.mixture_weights,[.5 .5]);
cfg.resume=true;
verifyError(testCase,@()dsc_run_prior_smoothing(testCase.TestData.panel, ...
    testCase.TestData.priors,cfg,folder),'dsc:PriorResumeUnsupported');
end

function testInvalidConfigsNeverCreateRunDirectory(testCase)
cfg=testCase.TestData.cfg; folder=new_run_folder(testCase); panel=testCase.TestData.panel; priors=testCase.TestData.priors;
cfg.auto_extend=true;
verifyError(testCase,@()dsc_run_prior_smoothing(panel,priors,cfg,folder),'dsc:PriorSmoothingConfig');
verifyFalse(testCase,isfolder(folder));
cfg=testCase.TestData.cfg; cfg.num_chains=0;
verifyError(testCase,@()dsc_run_prior_smoothing(panel,priors,cfg,folder),'dsc:PriorSmoothingConfig');
cfg=testCase.TestData.cfg; cfg.prior_parameter_draws=1.5;
verifyError(testCase,@()dsc_run_prior_smoothing(panel,priors,cfg,folder),'dsc:PriorSmoothingConfig');
cfg=testCase.TestData.cfg; cfg.max_iterations=cfg.burnin;
verifyError(testCase,@()dsc_run_prior_smoothing(panel,priors,cfg,folder),'dsc:PriorSmoothingConfig');
cfg=testCase.TestData.cfg; cfg.seed=2^32-2;
verifyError(testCase,@()dsc_run_prior_smoothing(panel,priors,cfg,folder),'dsc:PriorSmoothingConfig');
verifyFalse(testCase,isfolder(folder));
end

function testFailurePreservesDrawsWithoutRetryOrDroppingComponent(testCase)
cfg=testCase.TestData.cfg; cfg.correlation_backend='deliberately_invalid';
folder=new_run_folder(testCase);
result=dsc_run_prior_smoothing(testCase.TestData.panel,testCase.TestData.priors,cfg,folder);
verifyEqual(testCase,result.status,'failed');
verifyFalse(testCase,result.prior_sampling.mixture_complete);
verifyEqual(testCase,result.prior_sampling.completed_parameter_draws,0);
verifyEqual(testCase,result.saved_draws_per_chain,[0 0 0 0]);
verifyNotEmpty(testCase,result.failure_identifier);
loaded=load(result.paths.prior_parameters,'prior_parameters');
expected=dsc_draw_prior_parameters(testCase.TestData.priors,2,2,mod(cfg.seed+104729,2^32));
verifyEqual(testCase,loaded.prior_parameters.prior_draws,expected);
verifyEqual(testCase,result.prior_sampling.conditional_convergence{1}.status,'failed');
verifyEqual(testCase,result.prior_sampling.conditional_convergence{2}.status,'not_started');
verifyTrue(testCase,isfile(result.paths.result));
verifyTrue(testCase,isfile(result.paths.prior_sampling_metadata));
end

function folder=new_run_folder(testCase)
folder=tempname;
testCase.addTeardown(@()cleanup_run(folder));
end

function cleanup_run(folder)
if isfolder(folder), rmdir(folder,'s'); end
end
