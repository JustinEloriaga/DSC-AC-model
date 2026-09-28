function tests=test_parameter_files
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'toolbox'));
rng(712,'twister'); T=32; m=3;
panel=struct('returns',randn(T,m),'dates',(datetime(2003,8,15)+calweeks(0:T-1))', ...
    'tickers',["A","B","C"],'observation_mask',true(T,m),'source_hash','original-source');
[panel.pair_i,panel.pair_j]=find(tril(true(m),-1));
cfg=weekly_config(); cfg.prior_weeks=16;
cfg.calibration_windows=[12 16 20]; cfg.calibration_window=16;
cfg.max_iterations=4; cfg.burnin=1; cfg.thin=1;
priors=dsc_calibrate_priors(panel,cfg);
estimate=struct('V',[.03 .004 0;.004 .02 0;0 0 .01], ...
    'sig2h',[.01 .02 .03],'sig2r',[.004 .003 .002]);
run_dir=tempname; mkdir(run_dir);
training=subset_panel(panel,cfg.prior_weeks+1:T);
write_parameter_chunks(run_dir,estimate,3,training,priors,cfg);
result=struct('status','iteration_limit','parameter_estimate',estimate, ...
    'parameter_draws',3,'saved_draws',3,'run_dir',run_dir);
testCase.addTeardown(@() rmdir(run_dir,'s'));
testCase.TestData.panel=panel; testCase.TestData.priors=priors;
testCase.TestData.calibration=subset_panel(panel,1:cfg.prior_weeks);
testCase.TestData.training=subset_panel(panel,cfg.prior_weeks+1:T);
testCase.TestData.cfg=cfg; testCase.TestData.result=result;
end

function setup(testCase)
folder=tempname; mkdir(folder);
testCase.TestData.folder=folder;
testCase.addTeardown(@() rmdir(folder,'s'));
end

function testRoundTripAndAtomicPublish(testCase)
path=write_bundle(testCase);
b=dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel);
verifyEqual(testCase,b.parameter_estimate,testCase.TestData.result.parameter_estimate,'AbsTol',1e-14);
verifyEqual(testCase,b.priors,testCase.TestData.priors);
verifyEqual(testCase,b.schema_version,4);
verifyEqual(testCase,b.parameter_draws,load_parameter_draws(testCase.TestData.result.run_dir));
verifyEqual(testCase,b.training_returns,testCase.TestData.training.returns);
verifyEqual(testCase,b.training_mask,testCase.TestData.training.observation_mask);
verifyEqual(testCase,b.training_dates,testCase.TestData.training.dates);
verifyEqual(testCase,b.calibration_returns,testCase.TestData.calibration.returns);
verifyEqual(testCase,b.calibration_mask,testCase.TestData.calibration.observation_mask);
verifyEqual(testCase,b.calibration_dates,testCase.TestData.calibration.dates);
verifyEqual(testCase,b.calibration_weeks,16);
verifyEqual(testCase,b.estimation_start,b.calibration_end+calweeks(1));
verifyEqual(testCase,b.retained_draws,3);
verifyFalse(testCase,b.convergence_established);
verifyEqual(testCase,b.path,path);
verifyEqual(testCase,b.estimation_end,testCase.TestData.panel.dates(end));
verifyNotEmpty(testCase,regexp(path,'dsc_parameters_\d{8}T\d{9}Z_.*\.mat$','once'));
verifyEmpty(testCase,dir(fullfile(testCase.TestData.folder,'parameters_pending_*.mat')));
end

function testAppendPermitsChangedSourceHashKeepsPriors(testCase)
write_bundle(testCase);
p=testCase.TestData.panel;
p.dates(end+1)=p.dates(end)+calweeks(1);
p.returns(end+1,:)=[.4 -.2 .1]; p.observation_mask(end+1,:)=true;
p.source_hash='source-with-new-week';
b=dsc_load_parameters(testCase.TestData.folder,p);
verifyEqual(testCase,b.priors,testCase.TestData.priors);
verifyEqual(testCase,b.source_hash,'original-source');
verifyEqual(testCase,b.estimation_end,p.dates(end-1));
end

function testLatestUsesSavedTimeAndExplicitCanChooseOlder(testCase)
newer=write_bundle(testCase); older=write_bundle(testCase);
% Rewrite the older estimate last: its file modification time is newest.
shift_timestamp(newer,100);
shift_timestamp(older,-100);
loaded=dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel);
verifyEqual(testCase,loaded.path,newer);
chosen=dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel,older);
verifyEqual(testCase,chosen.path,older);
end

function testRejectsRevisedHistoricalReturnsDatesOrMask(testCase)
write_bundle(testCase);
for row=[2 18] % Revisions in either calibration or likelihood history fail.
    p=testCase.TestData.panel; p.returns(row,1)=p.returns(row,1)+1e-10;
    verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,p),'dsc:ParametersHistory');
    p=testCase.TestData.panel; p.dates(row)=p.dates(row)+days(1);
    verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,p),'dsc:ParametersHistory');
    p=testCase.TestData.panel; p.observation_mask(row,1)=false;
    verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,p),'dsc:ParametersHistory');
end
end

function testRejectsSeriesOrPairOrderMismatch(testCase)
write_bundle(testCase); p=testCase.TestData.panel;
p.tickers=p.tickers([2 1 3]);
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,p),'dsc:ParametersIdentity');
p=testCase.TestData.panel; p.pair_i=flip(p.pair_i);
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,p),'dsc:ParametersIdentity');
end

function testRejectsTruncatedHistory(testCase)
write_bundle(testCase); p=testCase.TestData.panel;
p.dates=p.dates(1:end-1); p.returns=p.returns(1:end-1,:);
p.observation_mask=p.observation_mask(1:end-1,:);
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,p),'dsc:ParametersPanel');
end

function testNewestInvalidBundleDoesNotFallback(testCase)
write_bundle(testCase); latest=write_bundle(testCase);
shift_timestamp(latest,100);
data=load(latest,'parameters'); parameters=data.parameters;
parameters=rmfield(parameters,'priors');
save(latest,'parameters','-v7');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel), ...
    'dsc:ParametersSchema');
end

function testRejectsUnknownSchemaInvalidCovarianceAndVariances(testCase)
path=write_bundle(testCase); data=load(path,'parameters'); original=data.parameters;
parameters=original; parameters.schema_version=5;
save(path,'parameters','-v7');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel), ...
    'dsc:ParametersSchema');
parameters=original; parameters.parameter_estimate.V(1,1)=-1;
save(path,'parameters','-v7');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel), ...
    'dsc:ParametersValues');
parameters=original; parameters.parameter_estimate.sig2r(1)=0;
save(path,'parameters','-v7');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel), ...
    'dsc:ParametersValues');
parameters=original; parameters.priors.source_hash='wrong-prior-source';
save(path,'parameters','-v7');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel), ...
    'dsc:ParametersPriors');
end

function testIgnoresOtherMatFilesAndRejectsLegacyExplicitFile(testCase)
parameters=struct('last_state',struct('V',eye(3))); %#ok<NASGU>
legacy=fullfile(testCase.TestData.folder,'checkpoint.mat'); save(legacy,'parameters','-v7');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel), ...
    'dsc:ParametersMissing');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel,legacy), ...
    'dsc:ParametersSchema');
end

function testCannotSaveWarmupTooFewOrFailedDraws(testCase)
result=testCase.TestData.result; result.parameter_draws=0; result.saved_draws=0;
verifyError(testCase,@()save_result(testCase,result),'dsc:ParametersDraws');
result.parameter_draws=1; result.saved_draws=1;
verifyError(testCase,@()save_result(testCase,result),'dsc:ParametersDraws');
result=testCase.TestData.result; result.status='failed';
verifyError(testCase,@()save_result(testCase,result),'dsc:ParametersFailedRun');
result=rmfield(testCase.TestData.result,'parameter_estimate');
verifyError(testCase,@()save_result(testCase,result),'dsc:ParametersResult');
verifyEmpty(testCase,dir(fullfile(testCase.TestData.folder,'dsc_parameters_*.mat')));
end

function testRejectsInvalidDrawsBeforePublishing(testCase)
path=fullfile(testCase.TestData.result.run_dir,'posterior_chunk_000001.mat');
loaded=load(path,'chunk'); original=loaded.chunk; chunk=original;
chunk.sig2h(1)=NaN;
save(path,'chunk','-v7');
testCase.addTeardown(@() restore_chunk(path,original));
result=testCase.TestData.result;
verifyError(testCase,@()save_result(testCase,result),'dsc:ParametersValues');
verifyEmpty(testCase,dir(fullfile(testCase.TestData.folder,'*.mat')));
end

function testSmoothingCannotPublishParameters(testCase)
result=testCase.TestData.result; result.estimate_parameters=false;
verifyError(testCase,@()save_result(testCase,result),'dsc:ParametersFixedMode');
result=testCase.TestData.result; result.inference_mode='fixed_parameter_smoothing';
verifyError(testCase,@()save_result(testCase,result),'dsc:ParametersFixedMode');
cfg=testCase.TestData.cfg; cfg.estimate_parameters=false;
verifyError(testCase,@()dsc_save_parameters(testCase.TestData.result,testCase.TestData.training, ...
    testCase.TestData.priors,cfg,testCase.TestData.folder,testCase.TestData.calibration),'dsc:ParametersFixedMode');
verifyEmpty(testCase,dir(fullfile(testCase.TestData.folder,'*.mat')));
end

function testVarianceVectorsAreCanonicalRows(testCase)
result=testCase.TestData.result;
result.parameter_estimate.sig2h=result.parameter_estimate.sig2h';
result.parameter_estimate.sig2r=result.parameter_estimate.sig2r';
path=save_result(testCase,result);
data=load(path,'parameters');
verifySize(testCase,data.parameters.parameter_estimate.sig2h,[1 3]);
verifySize(testCase,data.parameters.parameter_estimate.sig2r,[1 3]);
parameters=data.parameters;
parameters.parameter_estimate.sig2h=parameters.parameter_estimate.sig2h';
parameters.parameter_estimate.sig2r=parameters.parameter_estimate.sig2r';
save(path,'parameters','-v7');
b=dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel);
verifySize(testCase,b.parameter_estimate.sig2h,[1 3]);
verifySize(testCase,b.parameter_estimate.sig2r,[1 3]);
end

function testSchemaOneRequiresReestimation(testCase)
path=write_bundle(testCase); loaded=load(path,'parameters'); parameters=loaded.parameters;
parameters.schema_version=1;
parameters=rmfield(parameters,{'calibration_dates','calibration_returns','calibration_mask', ...
    'calibration_weeks','calibration_start','calibration_end'});
save(path,'parameters','-v7');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel), ...
    'dsc:ParametersLegacyCalibration');
try
    dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel);
catch problem
    verifySubstring(testCase,problem.message,'included calibration weeks in the likelihood');
    verifySubstring(testCase,problem.message,'Re-estimate');
end
end

function testSchemaTwoMeanOnlyRequiresReestimation(testCase)
path=write_bundle(testCase); loaded=load(path,'parameters'); parameters=loaded.parameters;
parameters.schema_version=2;
parameters.estimator='posterior_mean';
parameters=rmfield(parameters,'parameter_draws');
save(path,'parameters','-v7');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel), ...
    'dsc:ParametersLegacyDraws');
end

function testSchemaThreeLoadsAsNotChecked(testCase)
path=write_bundle(testCase); loaded=load(path,'parameters'); parameters=loaded.parameters;
parameters.schema_version=3;
parameters=rmfield(parameters,{'chain_ids','chain_dirs','chain_seeds','retained_draws_per_chain','convergence'});
parameters.parameter_draws=rmfield(parameters.parameter_draws,{'chain_id','iteration'});
save(path,'parameters','-v7');
b=dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel,path);
verifyEqual(testCase,b.schema_version,3); verifyFalse(testCase,b.convergence_established);
verifyEqual(testCase,b.convergence.status,'not_checked'); verifyFalse(testCase,b.convergence.passed);
verifyEqual(testCase,b.parameter_estimate,testCase.TestData.result.parameter_estimate,'AbsTol',1e-14);
parameters.convergence_established=true; save(path,'parameters','-v7');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel,path), ...
    'dsc:ParametersSchema');
end

function testMultichainParameterPoolingPreservesDrawIdentity(testCase)
result=multichain_fixture(testCase);
path=save_result(testCase,result);
b=dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel,path);
verifyEqual(testCase,b.schema_version,4);
verifyEqual(testCase,b.chain_ids,[7 9]); verifyEqual(testCase,b.retained_draws_per_chain,[2 4]);
verifyEqual(testCase,b.chain_seeds,testCase.TestData.cfg.seed+[6 8]);
verifyEqual(testCase,b.parameter_draws.chain_id,[7;7;9;9;9;9]);
verifyEqual(testCase,b.parameter_draws.iteration,[2;3;2;3;4;5]);
verifyEqual(testCase,b.retained_draws,6);
verifyEqual(testCase,b.parameter_estimate,result.parameter_estimate,'AbsTol',1e-14);
verifyEqual(testCase,b.parameter_estimate.sig2h,mean(b.parameter_draws.sig2h,1),'AbsTol',1e-14);
verifyFalse(testCase,b.convergence_established);
end

function testSaveRejectsDuplicateChainsAndAlteredSourceIdentity(testCase)
result=multichain_fixture(testCase); duplicated=result;
duplicated.chain_dirs=result.chain_dirs([1 1]);
verifyError(testCase,@()save_result(testCase,duplicated),'dsc:ParametersIdentity');
cp_path=fullfile(result.chain_dirs{2},'checkpoint.mat'); loaded=load(cp_path,'checkpoint'); original=loaded.checkpoint;
checkpoint=original; checkpoint.identity.priors.ig_scale_h=2*checkpoint.identity.priors.ig_scale_h;
save(cp_path,'checkpoint','-v7');
verifyError(testCase,@()save_result(testCase,result),'dsc:ParametersIdentity');
checkpoint=original; save(cp_path,'checkpoint','-v7');
chunk_path=fullfile(result.chain_dirs{2},'posterior_chunk_000001.mat'); loaded=load(chunk_path,'chunk'); chunk=loaded.chunk;
chunk.tickers=chunk.tickers([2 1 3]); save(chunk_path,'chunk','-v7');
verifyError(testCase,@()save_result(testCase,result),'dsc:ParametersIdentity');
verifyEmpty(testCase,dir(fullfile(testCase.TestData.folder,'dsc_parameters_*.mat')));
end

function testRejectsForgedConvergenceAndPerDrawProvenance(testCase)
result=multichain_fixture(testCase);
result.convergence=struct('status','passed','passed',true);
verifyError(testCase,@()save_result(testCase,result),'dsc:ConvergenceMetadata');
result=rmfield(result,'convergence'); path=save_result(testCase,result);
loaded=load(path,'parameters'); original=loaded.parameters;
parameters=original; parameters.parameter_draws.chain_id(1)=9; save(path,'parameters','-v7');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel,path), ...
    'dsc:ParametersIdentity');
parameters=original; parameters.parameter_draws.iteration(2)=parameters.parameter_draws.iteration(1);
save(path,'parameters','-v7');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel,path), ...
    'dsc:ParametersIdentity');
parameters=original; parameters.convergence_established=true; save(path,'parameters','-v7');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel,path), ...
    'dsc:ConvergenceMetadata');
end

function testRejectsInvalidSeedMetadata(testCase)
path=write_bundle(testCase); loaded=load(path,'parameters'); original=loaded.parameters;
for seed=[-1 2^32 .5 original.chain_seeds+1]
    parameters=original; parameters.chain_seeds=seed; save(path,'parameters','-v7');
    verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel,path), ...
        'dsc:ParametersIdentity');
end
for seed=[-1 2^32 .5]
    parameters=original; parameters.estimation_config.seed=seed;
    save(path,'parameters','-v7');
    verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel,path), ...
        'dsc:ParametersIdentity');
end
parameters=original; parameters.estimation_config.chain_id=0;
save(path,'parameters','-v7');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel,path), ...
    'dsc:ParametersIdentity');
end

function testRejectsInvalidWarmupThinningAndConfiguredChainBase(testCase)
path=write_bundle(testCase); loaded=load(path,'parameters'); original=loaded.parameters;
for thin=[0 -1 .5 NaN]
    parameters=original; parameters.estimation_config.thin=thin;
    % A zero thinning value previously accepted repeated retained iterations.
    if thin==0, parameters.parameter_draws.iteration(:)=parameters.estimation_config.burnin; end
    save(path,'parameters','-v7');
    verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel,path), ...
        'dsc:ParametersIdentity');
end
for burnin=[-1 .5 NaN]
    parameters=original; parameters.estimation_config.burnin=burnin;
    save(path,'parameters','-v7');
    verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel,path), ...
        'dsc:ParametersIdentity');
end
parameters=original; parameters.estimation_config.num_chains=2;
save(path,'parameters','-v7');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel,path), ...
    'dsc:ParametersIdentity');
end

function testCalibrationBlockIsMandatoryAndExcluded(testCase)
verifyError(testCase,@()dsc_save_parameters(testCase.TestData.result,testCase.TestData.training, ...
    testCase.TestData.priors,testCase.TestData.cfg,testCase.TestData.folder),'dsc:ParametersCalibration');
path=write_bundle(testCase);
b=dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel,path);
verifyEmpty(testCase,intersect(b.calibration_dates,b.training_dates));
verifyEqual(testCase,b.training_dates(1),testCase.TestData.panel.dates(17));
verifyEqual(testCase,b.calibration_dates(end),testCase.TestData.panel.dates(16));
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.training), ...
    'dsc:ParametersPanel');
end

function testRejectsGapOverlapAndIncorrectCalibrationCount(testCase)
cal=testCase.TestData.calibration; training=testCase.TestData.training;
cfg=testCase.TestData.cfg; priors=testCase.TestData.priors;
for offset=[-1 1]
    changed=training; changed.dates=training.dates+calweeks(offset);
    verifyError(testCase,@()dsc_save_parameters(testCase.TestData.result,changed,priors,cfg, ...
        testCase.TestData.folder,cal),'dsc:ParametersIdentity');
end
changed=cal; changed.dates(2)=changed.dates(2)+days(1);
verifyError(testCase,@()dsc_save_parameters(testCase.TestData.result,training,priors,cfg, ...
    testCase.TestData.folder,changed),'dsc:ParametersCalibration');
cfg.prior_weeks=15;
verifyError(testCase,@()dsc_save_parameters(testCase.TestData.result,training,priors,cfg, ...
    testCase.TestData.folder,cal),'dsc:ParametersCalibration');
cfg=testCase.TestData.cfg; priors.initial_weeks=15;
verifyError(testCase,@()dsc_save_parameters(testCase.TestData.result,training,priors,cfg, ...
    testCase.TestData.folder,cal),'dsc:ParametersCalibration');
verifyEmpty(testCase,dir(fullfile(testCase.TestData.folder,'*.mat')));
end

function testRejectsInvalidStoredBlocksAndBoundary(testCase)
path=write_bundle(testCase); data=load(path,'parameters'); original=data.parameters;
parameters=original; parameters.calibration_returns(1,1)=NaN;
save(path,'parameters','-v7');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel), ...
    'dsc:ParametersSchema');
parameters=original; parameters.training_mask=parameters.training_mask(1:end-1,:);
save(path,'parameters','-v7');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel), ...
    'dsc:ParametersSchema');
parameters=original; parameters.training_dates=parameters.training_dates+calweeks(1);
parameters.estimation_start=parameters.training_dates(1); parameters.estimation_end=parameters.training_dates(end);
save(path,'parameters','-v7');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel), ...
    'dsc:ParametersCalibration');
end

function path=write_bundle(testCase)
path=save_result(testCase,testCase.TestData.result);
end

function path=save_result(testCase,result)
path=dsc_save_parameters(result,testCase.TestData.training,testCase.TestData.priors, ...
    testCase.TestData.cfg,testCase.TestData.folder,testCase.TestData.calibration);
end

function out=subset_panel(panel,rows)
out=panel; out.dates=panel.dates(rows); out.returns=panel.returns(rows,:);
out.observation_mask=panel.observation_mask(rows,:);
end

function shift_timestamp(path,seconds_to_add)
data=load(path,'parameters'); parameters=data.parameters;
parameters.estimated_at_posix=parameters.estimated_at_posix+seconds_to_add;
date=datetime(parameters.estimated_at_posix,'ConvertFrom','posixtime','TimeZone','UTC');
parameters.estimated_at=char(string(date,"yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
save(path,'parameters','-v7');
end

function restore_chunk(path,chunk)
save(path,'chunk','-v7');
end

function write_parameter_chunks(run_dir,estimate,count,panel,priors,cfg)
scale=linspace(.9,1.1,count);
chunk=struct();
chunk.iterations=cfg.burnin+cfg.thin*(1:count);
chunk.dates=panel.dates; chunk.tickers=panel.tickers;
chunk.chain_id=cfg.chain_id; chunk.pair_i=panel.pair_i; chunk.pair_j=panel.pair_j;
chunk.V=zeros(size(estimate.V,1),size(estimate.V,2),count);
chunk.sig2h=zeros(count,numel(estimate.sig2h));
chunk.sig2r=zeros(count,numel(estimate.sig2r));
for k=1:count
    chunk.V(:,:,k)=estimate.V*scale(k);
    chunk.sig2h(k,:)=estimate.sig2h*scale(k);
    chunk.sig2r(k,:)=estimate.sig2r*scale(k);
end
save(fullfile(run_dir,'posterior_chunk_000001.mat'),'chunk','-v7');
identity=struct('dates',panel.dates,'tickers',panel.tickers,'returns',panel.returns, ...
    'mask',logical(panel.observation_mask)&isfinite(panel.returns),'source_hash',panel.source_hash, ...
    'cfg',cfg,'priors',priors);
checkpoint=struct('identity',identity,'saved',count,'completed',chunk.iterations(end));
save(fullfile(run_dir,'checkpoint.mat'),'checkpoint','-v7');
end

function draws=load_parameter_draws(run_dir)
loaded=load(fullfile(run_dir,'posterior_chunk_000001.mat'),'chunk');
draws=struct('V',loaded.chunk.V,'sig2h',loaded.chunk.sig2h,'sig2r',loaded.chunk.sig2r, ...
    'chain_id',repmat(loaded.chunk.chain_id,numel(loaded.chunk.iterations),1),'iteration',loaded.chunk.iterations(:));
end

function result=multichain_fixture(testCase)
result=testCase.TestData.result; result.run_dir=fullfile(testCase.TestData.folder,'source_chains');
mkdir(result.run_dir); result.chain_dirs=cell(1,2);
counts=[2 4]; ids=[7 9]; factors=[1 2]; base=result.parameter_estimate;
for c=1:2
    result.chain_dirs{c}=fullfile(result.run_dir,sprintf('chain_%d',ids(c))); mkdir(result.chain_dirs{c});
    cfg=testCase.TestData.cfg; cfg.chain_id=ids(c); cfg.max_iterations=cfg.burnin+counts(c)*cfg.thin;
    estimate=struct('V',base.V*factors(c),'sig2h',base.sig2h*factors(c),'sig2r',base.sig2r*factors(c));
    write_parameter_chunks(result.chain_dirs{c},estimate,counts(c),testCase.TestData.training,testCase.TestData.priors,cfg);
end
factor=sum(counts.*factors)/sum(counts);
result.parameter_estimate=struct('V',base.V*factor,'sig2h',base.sig2h*factor,'sig2r',base.sig2r*factor);
result.parameter_draws=sum(counts); result.saved_draws=sum(counts);
end
