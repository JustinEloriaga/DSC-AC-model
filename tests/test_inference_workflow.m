function tests=test_inference_workflow
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'toolbox'));
rng(773,'twister'); T=32; m=3;
panel=struct('returns',randn(T,m),'dates',(datetime(2020,1,3)+calweeks(0:T-1))', ...
    'tickers',["A","B","C"],'observation_mask',true(T,m),'source_hash','workflow-synthetic');
[panel.pair_i,panel.pair_j]=find(tril(true(m),-1));
cfg=weekly_config(); cfg.output_root=tempname; mkdir(cfg.output_root);
cfg.prior_weeks=16; cfg.calibration_windows=[12 16 20]; cfg.calibration_window=16;
cfg.burnin=1; cfg.max_iterations=4; cfg.max_seconds=90;
cfg.chunk_size=2; cfg.correlation_backend='matlab';
testCase.TestData.panel=panel; testCase.TestData.cfg=cfg;
end

function testCutoffPriorsIgnoreLaterReturns(testCase)
cfg=testCase.TestData.cfg; a=testCase.TestData.panel;
verifyEqual(testCase,cfg.parameter_smoothing,'mean');
cfg.estimation_end_date=char(string(a.dates(24),'yyyy-MM-dd'));
b=a; b.returns(25:end,:)=1000*b.returns(25:end,:);
[pa,ia]=dsc_prepare_inference(a,cfg); [pb,ib]=dsc_prepare_inference(b,cfg);
verifyEqual(testCase,pa,pb); verifyEqual(testCase,ia.training_panel,ib.training_panel);
verifyEqual(testCase,numel(ia.training_panel.dates),8);
verifyEqual(testCase,ia.training_panel.dates,a.dates(17:24));
verifyEqual(testCase,ia.training_panel.returns,a.returns(17:24,:));
verifyEqual(testCase,ia.calibration_weeks,16);
verifyEqual(testCase,ia.calibration_panel.dates,a.dates(1:16));
verifyEqual(testCase,ia.calibration_panel.returns,a.returns(1:16,:));
verifyEqual(testCase,ia.smoothing_panel.dates,a.dates(17:32));
cfg_without_smoothing=rmfield(cfg,'parameter_smoothing');
[~,default_inference]=dsc_prepare_inference(a,cfg_without_smoothing);
verifyEqual(testCase,default_inference.parameter_smoothing,'mean');
cfg.parameter_smoothing='';
[~,empty_inference]=dsc_prepare_inference(a,cfg);
verifyEqual(testCase,empty_inference.parameter_smoothing,'mean');
verifyEqual(testCase,pa.initial_dates,a.dates([1 16]));
verifyEqual(testCase,pa.excluded_initial_weeks,16);
verifyEqual(testCase,pa.likelihood_start_date,a.dates(17));
verifyError(testCase,@()dsc_sample(a,pa,cfg,tempname),'dsc:CalibrationOverlap');
end

function testEstimateThenExtendWithFrozenParameters(testCase)
cfg=testCase.TestData.cfg; panel=testCase.TestData.panel;
cfg.output_root=tempname; mkdir(cfg.output_root);
cfg.estimation_end_date=char(string(panel.dates(24),'yyyy-MM-dd'));
% Missing calibration cells must not leave those dates in the state path.
% Missing later cells remain on their original inference dates.
panel.observation_mask(2:3,2)=false; panel.returns(2:3,2)=NaN;
panel.observation_mask(20:21,2)=false; panel.returns(20:21,2)=NaN;
[priors,inference]=dsc_prepare_inference(panel,cfg);
first=dsc_run_inference(panel,priors,cfg,fullfile(cfg.output_root,'first'),inference);
verifyEqual(testCase,first.inference_mode,'fixed_parameter_smoothing');
verifyEqual(testCase,first.T,16);
verifyEqual(testCase,first.dates,panel.dates(17:end));
verifyEqual(testCase,first.inference.parameter_estimation_end,cfg.estimation_end_date);
verifyEqual(testCase,first.inference.parameter_estimation_start,date_text(panel.dates(17)));
verifyEqual(testCase,first.inference.calibration_weeks,16);
verifyEqual(testCase,first.inference.calibration_start,date_text(panel.dates(1)));
verifyEqual(testCase,first.inference.calibration_end,date_text(panel.dates(16)));
verifyEqual(testCase,first.inference.full_data_start,date_text(panel.dates(1)));
verifyEqual(testCase,first.inference.smoothing_start,date_text(panel.dates(17)));
verifyEqual(testCase,first.inference.parameter_estimation_draws,3);
verifyEqual(testCase,first.inference.parameter_smoothing,'mean');
verifyFalse(testCase,first.inference.parameter_uncertainty_in_bands);
verifyTrue(testCase,isfile(first.paths.parameters));
cp=verify_likelihood_rows(testCase, ...
    fullfile(cfg.output_root,'first','estimation','checkpoint.mat'),panel,17:24);
verifyEqual(testCase,cp.identity.priors,priors);
verifyEqual(testCase,nnz(cp.identity.mask),22);
cp=verify_likelihood_rows(testCase,first.paths.checkpoint,panel,17:32);
verifyEqual(testCase,nnz(cp.identity.mask),46);
stored=load(first.paths.parameters,'parameters');
verifyEqual(testCase,stored.parameters.schema_version,4);
verifyTrue(testCase,isfield(stored.parameters,'parameter_draws'));
verifyEqual(testCase,stored.parameters.training_dates,panel.dates(17:24));
verifyEqual(testCase,cp.state.V,stored.parameters.parameter_estimate.V);
verifyEqual(testCase,cp.state.sig2h,stored.parameters.parameter_estimate.sig2h);
verifyEqual(testCase,cp.state.sig2r,stored.parameters.parameter_estimate.sig2r);

extended=panel; extended.dates(end+1)=panel.dates(end)+calweeks(1);
extended.returns(end+1,:)=[.2 -.3 .4]; extended.observation_mask(end+1,:)=true;
extended.source_hash='appended-file-hash';
cfg.estimate_parameters=false; cfg.estimation_end_date='';
cfg=rmfield(cfg,'parameter_smoothing'); % Loading also defaults to mean smoothing.
% These impossible calibration settings must be ignored when loading.
cfg.prior_weeks=1; cfg.calibration_window=999;
[loaded_priors,loaded_inference]=dsc_prepare_inference(extended,cfg);
verifyEqual(testCase,loaded_priors.V0B,priors.V0B);
verifyEqual(testCase,loaded_priors.estimation_source_hash,panel.source_hash);
verifyEqual(testCase,loaded_inference.calibration_weeks,16);
verifyEqual(testCase,loaded_inference.calibration_panel.dates,panel.dates(1:16));
verifyEqual(testCase,loaded_inference.smoothing_panel.dates,extended.dates(17:end));
second=dsc_run_inference(extended,loaded_priors,cfg,fullfile(cfg.output_root,'second'),loaded_inference);
verifyEqual(testCase,second.T,17);
verifyEqual(testCase,second.dates,extended.dates(17:end));
verifyEqual(testCase,second.paths.parameters,first.paths.parameters);
verifyEqual(testCase,second.inference.parameters_estimated_at,first.inference.parameters_estimated_at);
verifyEqual(testCase,second.inference.smoothing_start,date_text(panel.dates(17)));
verifyEqual(testCase,second.inference.smoothing_end,date_text(extended.dates(end)));
verifyEqual(testCase,second.inference.parameter_smoothing,'mean');
verifyFalse(testCase,second.inference.parameter_uncertainty_in_bands);
cp=verify_likelihood_rows(testCase,second.paths.checkpoint,extended,17:33);
verifyEqual(testCase,cp.state.V,loaded_inference.parameters.parameter_estimate.V);
verifyEqual(testCase,cp.state.sig2h,loaded_inference.parameters.parameter_estimate.sig2h);
verifyEqual(testCase,cp.state.sig2r,loaded_inference.parameters.parameter_estimate.sig2r);
verifyEqual(testCase,numel(dir(fullfile(inference.parameter_dir,'dsc_parameters_*.mat'))),1);
manifest=jsondecode(fileread(fullfile(second.run_dir,'run_manifest.json')));
verifyEqual(testCase,manifest.inference.parameter_file,first.paths.parameters);
verifyEqual(testCase,manifest.inference.calibration_weeks,16);
verifyEqual(testCase,manifest.inference.smoothing_start,date_text(panel.dates(17)));
end

function testLatestDateDefaultsToMeanSmoothing(testCase)
cfg=testCase.TestData.cfg; cfg.output_root=tempname; mkdir(cfg.output_root);
panel=testCase.TestData.panel;
[priors,inference]=dsc_prepare_inference(panel,cfg);
run_dir=fullfile(cfg.output_root,'mean-default');
result=dsc_run_inference(panel,priors,cfg,run_dir,inference);
verifyEqual(testCase,result.inference_mode,'fixed_parameter_smoothing');
verifyFalse(testCase,result.inference.parameter_uncertainty_in_bands);
verifyEqual(testCase,result.inference.parameter_smoothing,'mean');
verifyEqual(testCase,result.inference.parameter_estimation_run,fullfile(run_dir,'estimation'));
verifyEqual(testCase,result.T,16);
verifyEqual(testCase,result.dates,panel.dates(17:32));
verifyEqual(testCase,result.inference.parameter_estimation_start,date_text(panel.dates(17)));
verifyEqual(testCase,result.inference.parameter_estimation_end,date_text(panel.dates(32)));
verify_likelihood_rows(testCase,result.paths.checkpoint,panel,17:32);
estimated=load(fullfile(run_dir,'estimation','checkpoint.mat'),'checkpoint');
verifyEqual(testCase,estimated.checkpoint.saved,3);
stored=load(result.paths.parameters,'parameters');
verifyTrue(testCase,isfield(stored.parameters,'parameter_draws'));
verifyEqual(testCase,stored.parameters.retained_draws,3);
smoothed=load(result.paths.checkpoint,'checkpoint');
verifyEqual(testCase,smoothed.checkpoint.state.V,stored.parameters.parameter_estimate.V);
verifyEqual(testCase,smoothed.checkpoint.state.sig2h,stored.parameters.parameter_estimate.sig2h);
verifyEqual(testCase,smoothed.checkpoint.state.sig2r,stored.parameters.parameter_estimate.sig2r);
verifyTrue(testCase,isfolder(fullfile(run_dir,'estimation')));
end

function testLatestDateExplicitDrawsUsesJointEstimation(testCase)
cfg=testCase.TestData.cfg; cfg.output_root=tempname; mkdir(cfg.output_root);
cfg.parameter_smoothing='draws';
panel=testCase.TestData.panel;
[priors,inference]=dsc_prepare_inference(panel,cfg);
result=dsc_run_inference(panel,priors,cfg,fullfile(cfg.output_root,'joint'),inference);
verifyEqual(testCase,result.inference_mode,'parameter_estimation');
verifyTrue(testCase,result.inference.parameter_uncertainty_in_bands);
verifyEqual(testCase,result.inference.parameter_smoothing,'draws');
verifyEqual(testCase,result.inference.parameter_estimation_run,result.run_dir);
verifyFalse(testCase,isfolder(fullfile(result.run_dir,'estimation')));
end

function testWarmupOnlyCannotSaveParameters(testCase)
cfg=testCase.TestData.cfg; cfg.output_root=tempname; mkdir(cfg.output_root);
cfg.burnin=cfg.max_iterations;
panel=testCase.TestData.panel;
[priors,inference]=dsc_prepare_inference(panel,cfg);
verifyError(testCase,@()dsc_run_inference(panel,priors,cfg, ...
    fullfile(cfg.output_root,'warmup'),inference),'dsc:TooFewParameterDraws');
verifyEmpty(testCase,dir(fullfile(inference.parameter_dir,'dsc_parameters_*.mat')));
end

function testInvalidModeOptions(testCase)
cfg=testCase.TestData.cfg; panel=testCase.TestData.panel;
cfg.estimate_parameters=false; cfg.estimation_end_date='2020-04-03';
verifyError(testCase,@()dsc_prepare_inference(panel,cfg),'dsc:InferenceConfig');
cfg.estimate_parameters=true; cfg.estimation_end_date='2099-01-01';
verifyError(testCase,@()dsc_prepare_inference(panel,cfg),'dsc:EstimationDate');
cfg.estimation_end_date=date_text(panel.dates(16));
verifyError(testCase,@()dsc_prepare_inference(panel,cfg),'dsc:EstimationDate');
cfg.estimation_end_date=date_text(panel.dates(15));
verifyError(testCase,@()dsc_prepare_inference(panel,cfg),'dsc:EstimationDate');
cfg.estimation_end_date=''; cfg.parameter_file='arbitrary.mat';
verifyError(testCase,@()dsc_prepare_inference(panel,cfg),'dsc:InferenceConfig');
end

function testTwoChainReportEstimateSmoothReloadPosteriorMeans(testCase)
cfg=convergence_workflow_config(testCase); panel=testCase.TestData.panel;
cfg.num_chains=2; cfg.convergence_mode='report'; cfg.parameter_smoothing='mean';
cfg.estimation_end_date=date_text(panel.dates(24));
[priors,inference]=dsc_prepare_inference(panel,cfg);
first=dsc_run_inference(panel,priors,cfg,fullfile(cfg.output_root,'two_chain_first'),inference);
stored=load(first.paths.parameters,'parameters'); parameters=stored.parameters;
verifyEqual(testCase,parameters.schema_version,4);
verifyEqual(testCase,parameters.retained_draws,6);
verifyEqual(testCase,parameters.retained_draws_per_chain,[3 3]);
verifyEqual(testCase,parameters.parameter_draws.chain_id,[1;1;1;2;2;2]);
verifyEqual(testCase,first.saved_draws_per_chain,[3 3]);
verifyEqual(testCase,first.inference.retained_draws_per_chain,[3 3]);
verifyEqual(testCase,first.inference.parameter_estimation_draws,6);
verifyEqual(testCase,first.inference.estimation_convergence.status,'insufficient_draws');
verifyEqual(testCase,first.inference.state_convergence.status,'insufficient_draws');
verifyEqual(testCase,first.inference.estimation_convergence.draws_per_chain,3);
verifyEqual(testCase,first.inference.state_convergence.draws_per_chain,3);
verifyEqual(testCase,first.inference.estimation_convergence.quantities_checked,48);
verifyEqual(testCase,first.inference.state_convergence.quantities_checked,96);
verifyFalse(testCase,first.inference.estimation_convergence.passed);
verifyFalse(testCase,first.inference.state_convergence.passed);
verifyFalse(testCase,first.inference.parameter_uncertainty_in_bands);
for c=1:2
    cp=verify_likelihood_rows(testCase,fullfile(first.chain_dirs{c},'checkpoint.mat'),panel,17:32);
    verifyEqual(testCase,cp.state.V,parameters.parameter_estimate.V);
    verifyEqual(testCase,cp.state.sig2h,parameters.parameter_estimate.sig2h);
    verifyEqual(testCase,cp.state.sig2r,parameters.parameter_estimate.sig2r);
end
cfg.estimate_parameters=false; cfg.estimation_end_date=''; cfg.parameter_file=first.paths.parameters;
[loaded_priors,loaded_inference]=dsc_prepare_inference(panel,cfg);
second=dsc_run_inference(panel,loaded_priors,cfg,fullfile(cfg.output_root,'two_chain_reload'),loaded_inference);
verifyEqual(testCase,second.paths.parameters,first.paths.parameters);
verifyEqual(testCase,second.inference.parameters_estimated_at,first.inference.parameters_estimated_at);
verifyEqual(testCase,second.inference.estimation_convergence,parameters.convergence);
verifyEqual(testCase,second.inference.state_convergence.status,'insufficient_draws');
verifyEqual(testCase,second.saved_draws_per_chain,[3 3]);
for c=1:2
    cp=verify_likelihood_rows(testCase,fullfile(second.chain_dirs{c},'checkpoint.mat'),panel,17:32);
    verifyEqual(testCase,cp.state.V,parameters.parameter_estimate.V);
    verifyEqual(testCase,cp.state.sig2h,parameters.parameter_estimate.sig2h);
    verifyEqual(testCase,cp.state.sig2r,parameters.parameter_estimate.sig2r);
end
verifyEqual(testCase,numel(dir(fullfile(inference.parameter_dir,'dsc_parameters_*.mat'))),1);
end

function testReportModeContinuesAfterFailedEstimationDiagnostics(testCase)
cfg=convergence_workflow_config(testCase); panel=testCase.TestData.panel;
cfg.num_chains=4; cfg.convergence_mode='report'; cfg.parameter_smoothing='mean';
[priors,inference]=dsc_prepare_inference(panel,cfg);
run_dir=fullfile(cfg.output_root,'report_failed_estimation');
result=dsc_run_inference(panel,priors,cfg,run_dir,inference);
verifyTrue(testCase,isfile(result.paths.parameters));
diagnostic=jsondecode(fileread(fullfile(run_dir,'estimation','convergence','convergence_diagnostics.json')));
verifyEqual(testCase,diagnostic.status,'insufficient_draws'); verifyFalse(testCase,diagnostic.passed);
verifyEqual(testCase,diagnostic.actual_draws_per_chain(:),[3;3;3;3]);
verifyEqual(testCase,result.inference.estimation_convergence.status,'insufficient_draws');
verifyEqual(testCase,result.inference.state_convergence.status,'insufficient_draws');
verifyEqual(testCase,result.saved_draws_per_chain,[3 3 3 3]);
for c=1:4
    cp_path=fullfile(run_dir,'estimation','chains',sprintf('chain_%03d',c),'checkpoint.mat');
    verifyTrue(testCase,isfile(cp_path)); loaded=load(cp_path,'checkpoint');
    verifyEqual(testCase,loaded.checkpoint.completed,4); verifyEqual(testCase,loaded.checkpoint.saved,3);
end
verifyTrue(testCase,isfolder(fullfile(run_dir,'chains')));
verifyTrue(testCase,isfile(fullfile(run_dir,'inference_metadata.json')));
verifyTrue(testCase,isfile(fullfile(run_dir,'result.mat')));
end

function testOffModeAndReportModeContinueWithUndiagnosedLoad(testCase)
cfg=convergence_workflow_config(testCase); panel=testCase.TestData.panel;
cfg.convergence_mode='off'; cfg.parameter_smoothing='mean';
[priors,inference]=dsc_prepare_inference(panel,cfg);
result=dsc_run_inference(panel,priors,cfg,fullfile(cfg.output_root,'off'),inference);
verifyTrue(testCase,isfile(result.paths.parameters));
stored=load(result.paths.parameters,'parameters'); original=stored.parameters;
verifyEqual(testCase,original.schema_version,4); verifyFalse(testCase,original.convergence_established);
verifyEqual(testCase,original.convergence.status,'not_checked');
verifyEqual(testCase,result.inference.estimation_convergence.status,'not_checked');
verifyEqual(testCase,result.inference.state_convergence.status,'not_checked');
verifyEqual(testCase,original.retained_draws,3); verifyEqual(testCase,result.saved_draws_per_chain,3);
cfg.convergence_mode='report'; cfg.num_chains=4; cfg.estimate_parameters=false;
cfg.parameter_file=result.paths.parameters;
[loaded_priors,loaded_inference]=dsc_prepare_inference(panel,cfg);
loaded_dir=fullfile(cfg.output_root,'report_undiagnosed_load');
loaded_result=dsc_run_inference(panel,loaded_priors,cfg,loaded_dir,loaded_inference);
verifyTrue(testCase,isfile(loaded_result.paths.result));
verifyEqual(testCase,loaded_result.inference.estimation_convergence.status,'not_checked');
verifyEqual(testCase,loaded_result.inference.state_convergence.status,'insufficient_draws');
stored=load(result.paths.parameters,'parameters'); verifyEqual(testCase,stored.parameters,original);
verifyEqual(testCase,numel(dir(fullfile(inference.parameter_dir,'dsc_parameters_*.mat'))),1);
% Conditional draw cycling is reported as not_checked but still completes.
cfg.parameter_smoothing='draws';
[loaded_priors,loaded_inference]=dsc_prepare_inference(panel,cfg);
draws_dir=fullfile(cfg.output_root,'report_undiagnosed_draws');
draws_result=dsc_run_inference(panel,loaded_priors,cfg,draws_dir,loaded_inference);
verifyTrue(testCase,isfile(draws_result.paths.result));
verifyEqual(testCase,draws_result.inference.estimation_convergence.status,'not_checked');
verifyEqual(testCase,draws_result.inference.state_convergence.status,'not_checked');
verifyTrue(testCase,contains(draws_result.inference.state_convergence.reason,'conditional transition kernel'));
end

function testUnsupportedConvergenceModeRejectedBeforeInferenceWrites(testCase)
cfg=convergence_workflow_config(testCase); panel=testCase.TestData.panel;
cfg.convergence_mode='required';
[priors,inference]=dsc_prepare_inference(panel,cfg);
run_dir=fullfile(cfg.output_root,'unsupported_mode');
verifyError(testCase,@()dsc_run_inference(panel,priors,cfg,run_dir,inference), ...
    'MATLAB:unrecognizedStringChoice');
verifyFalse(testCase,isfolder(run_dir));
verifyEmpty(testCase,dir(fullfile(inference.parameter_dir,'dsc_parameters_*.mat')));
end

function testPublicEntryRejectsRemovedModeBeforeDataAccess(testCase)
verifyError(testCase,@()run_weekly_model_report(ConvergenceMode="required", ...
    SourceFile="/missing/source.csv",RunTests=false,BuildReport=false,VerifyReport=false), ...
    'MATLAB:validators:mustBeMember');
end

function cfg=convergence_workflow_config(testCase)
cfg=testCase.TestData.cfg; cfg.output_root=tempname; mkdir(cfg.output_root);
testCase.addTeardown(@()rmdir(cfg.output_root,'s'));
cfg.parallel_chains=false; cfg.auto_extend=false;
end

function checkpoint=verify_likelihood_rows(testCase,path,panel,rows)
loaded=load(path,'checkpoint'); checkpoint=loaded.checkpoint;
verifyEqual(testCase,checkpoint.identity.dates,panel.dates(rows));
verifyTrue(testCase,isequaln(checkpoint.identity.returns,panel.returns(rows,:)));
verifyEqual(testCase,checkpoint.identity.mask,logical(panel.observation_mask(rows,:)));
verifyEqual(testCase,size(checkpoint.state.B,1),numel(rows));
verifyEqual(testCase,size(checkpoint.state.h,1),numel(rows));
verifyEqual(testCase,size(checkpoint.state.r,1),numel(rows));
verifyFalse(testCase,any(ismember(checkpoint.identity.dates,panel.dates(1:16))));
chunks=dir(fullfile(fileparts(path),'posterior_chunk_*.mat'));
for k=1:numel(chunks)
    saved=load(fullfile(chunks(k).folder,chunks(k).name),'chunk');
    verifyEqual(testCase,saved.chunk.dates,panel.dates(rows));
    verifyEqual(testCase,size(saved.chunk.P_pairs,1),numel(rows));
end
end

function value=date_text(date)
value=char(string(date,'yyyy-MM-dd'));
end
