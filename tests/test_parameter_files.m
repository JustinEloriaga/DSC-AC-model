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
write_parameter_chunks(run_dir,estimate,3);
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
verifyEqual(testCase,b.schema_version,3);
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
parameters=original; parameters.schema_version=4;
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

function testSchemaTwoMeanOnlyRequiresReestimation(testCase)
path=write_bundle(testCase); loaded=load(path,'parameters'); parameters=loaded.parameters;
parameters.schema_version=2;
parameters.estimator='posterior_mean';
parameters=rmfield(parameters,'parameter_draws');
save(path,'parameters','-v7');
verifyError(testCase,@()dsc_load_parameters(testCase.TestData.folder,testCase.TestData.panel), ...
    'dsc:ParametersLegacyDraws');
end
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
        testCase.TestData.folder,cal),'dsc:ParametersCalibration');
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

function write_parameter_chunks(run_dir,estimate,count)
scale=[.9 1 1.1];
chunk=struct();
chunk.iterations=2:count+1;
chunk.V=zeros(size(estimate.V,1),size(estimate.V,2),count);
chunk.sig2h=zeros(count,numel(estimate.sig2h));
chunk.sig2r=zeros(count,numel(estimate.sig2r));
for k=1:count
    chunk.V(:,:,k)=estimate.V*scale(k);
    chunk.sig2h(k,:)=estimate.sig2h*scale(k);
    chunk.sig2r(k,:)=estimate.sig2r*scale(k);
end
save(fullfile(run_dir,'posterior_chunk_000001.mat'),'chunk','-v7');
end

function draws=load_parameter_draws(run_dir)
loaded=load(fullfile(run_dir,'posterior_chunk_000001.mat'),'chunk');
draws=struct('V',loaded.chunk.V,'sig2h',loaded.chunk.sig2h,'sig2r',loaded.chunk.sig2r);
end
