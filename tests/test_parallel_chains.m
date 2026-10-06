function tests=test_parallel_chains
% Opt-in integration test: runtests('tests/test_parallel_chains.m').
% Kept outside ordinary acceptance because it starts two MATLAB processes.
tests=functiontests(localfunctions);
end

function testProcessChainsMatchSequentialDraws(testCase)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root,fullfile(root,'toolbox'));
assumeTrue(testCase,dsc_parallel_available(),'Parallel Computing Toolbox is required for this opt-in test.');
folder=tempname; mkdir(folder);
testCase.addTeardown(@()rmdir(folder,'s'));
pool=gcp('nocreate');
if isempty(pool)
    storage=fullfile(folder,'pool_jobs'); mkdir(storage);
    cluster=parcluster('Processes'); cluster.JobStorageLocation=storage;
    pool=parpool(cluster,2);
    testCase.addTeardown(@()delete(pool));
else
    assumeFalse(testCase,isa(pool,'parallel.ThreadPool'),'Existing user thread pool is left intact.');
    assumeGreaterThanOrEqual(testCase,pool.NumWorkers,2,'Existing user pool has fewer than two workers.');
end
rng(721,'twister'); T=8; m=2;
panel=struct('returns',.3*randn(T,m), ...
    'dates',(datetime(2005,1,7)+calweeks(0:T-1))','tickers',["A","B"], ...
    'observation_mask',true(T,m),'pair_i',2,'pair_j',1,'source_hash','synthetic-process-chains');
priors=struct('Bbar',[0;0],'VBbar',.01*eye(m),'nuB',12,'V0B',.001*eye(m), ...
    'mh0',log([.09;.09]),'mr0',0,'h0_scale',10,'r0_scale',10, ...
    'sig2h_mean',.01,'sig2r_mean',.01,'ig_shape',10, ...
    'ig_scale_h',.09,'ig_scale_r',.09,'tickers',panel.tickers,'source_hash',panel.source_hash);
backend='matlab'; if exist('dsc_corr_batch_mex','file')==3, backend='mex'; end
cfg=struct('seed',912,'chain_id',1,'p',0,'max_iterations',3,'max_seconds',120, ...
    'burnin',1,'thin',1,'chunk_size',2,'resume',false,'estimate_parameters',true, ...
    'correlation_backend',backend,'correlation_threads',1,'num_chains',2, ...
    'parallel_chains',false,'convergence_mode','report','check_every',1);
serial=dsc_run_chains(panel,priors,cfg,fullfile(folder,'serial'));
cfg.parallel_chains=true;
parallel=dsc_run_chains(panel,priors,cfg,fullfile(folder,'parallel'));
verifyTrue(testCase,parallel.parallel_chains);
verifyEqual(testCase,parallel.saved_draws_per_chain,[2 2]);
verifyEqual(testCase,parallel.saved_draws_per_chain,serial.saved_draws_per_chain);
verifyEqual(testCase,parallel.chain_ids,serial.chain_ids);
verifyEqual(testCase,parallel.seeds,[912 913]);
for c=1:2
    a=load(fullfile(serial.chain_dirs{c},'checkpoint.mat'),'checkpoint');
    b=load(fullfile(parallel.chain_dirs{c},'checkpoint.mat'),'checkpoint');
    verifyEqual(testCase,a.checkpoint.state,b.checkpoint.state);
    verifyEqual(testCase,a.checkpoint.rng,b.checkpoint.rng);
    verifyEqual(testCase,a.checkpoint.parameter_sums,b.checkpoint.parameter_sums);
    verifyEqual(testCase,a.checkpoint.diagnostics.loglik,b.checkpoint.diagnostics.loglik);
    serial_files=dir(fullfile(serial.chain_dirs{c},'posterior_chunk_*.mat'));
    parallel_files=dir(fullfile(parallel.chain_dirs{c},'posterior_chunk_*.mat'));
    verifyEqual(testCase,numel(serial_files),numel(parallel_files));
    for k=1:numel(serial_files)
        a=load(fullfile(serial_files(k).folder,serial_files(k).name),'chunk');
        b=load(fullfile(parallel_files(k).folder,parallel_files(k).name),'chunk');
        verifyEqual(testCase,a.chunk,b.chunk);
    end
end
% Verify the process workers also start from different, reproducible states.
initial_cfg=rmfield(cfg,{'num_chains','parallel_chains','convergence_mode','check_every'});
initial_cfg.randomize_initial_state=true; initial_cfg.max_seconds=1e-6;
initial_dirs={fullfile(folder,'initial_1'),fullfile(folder,'initial_2')};
futures=cell(1,2);
for c=1:2
    local_cfg=initial_cfg; local_cfg.chain_id=c;
    futures{c}=parfeval(pool,@dsc_sample,1,panel,priors,local_cfg,initial_dirs{c});
    future=futures{c}; testCase.addTeardown(@()cancel(future));
end
for c=1:2
    result=fetchOutputs(futures{c}); verifyEqual(testCase,result.completed_iterations,0);
end
a=load(fullfile(initial_dirs{1},'checkpoint.mat'),'checkpoint');
b=load(fullfile(initial_dirs{2},'checkpoint.mat'),'checkpoint');
verifyFalse(testCase,isequaln(a.checkpoint.state,b.checkpoint.state));
verifyFalse(testCase,isequal(a.checkpoint.state.B,b.checkpoint.state.B));
verifyFalse(testCase,isequal(a.checkpoint.state.r,b.checkpoint.state.r));
end
