function tests=test_prior_entry_point
tests=functiontests(localfunctions);
end

function testPriorThenEstimationOnlyPublicOptions(testCase)
% A small CSV exercises the actual entry point and its exported artifacts.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'toolbox'),fullfile(root,'scripts'));
folder=tempname; mkdir(folder);
testCase.addTeardown(@()rmdir(folder,'s'));
prior_path=path;
testCase.addTeardown(@()path(prior_path));
prior_warning=warning('off','econ:adftest:InvalidStatistic');
testCase.addTeardown(@()warning(prior_warning));
for name={'run_weekly_model_report.m','weekly_config.m','run_weekly_research.m'}
    copyfile(fullfile(root,name{1}),fullfile(folder,name{1}));
end
copyfile(fullfile(root,'toolbox'),fullfile(folder,'toolbox'));
copyfile(fullfile(root,'scripts'),fullfile(folder,'scripts'));
addpath(folder,fullfile(folder,'toolbox'),fullfile(folder,'scripts'));
cfg=weekly_config(); original=weekly_config_from_root(root);
header=readlines(original.source_file);
dates=(datetime(2003,8,4):days(1):datetime(2005,10,7))';
dates=dates(~ismember(weekday(dates),[1 7]));
rng(139,'twister'); levels=100*exp(cumsum(.002*randn(numel(dates),14),1));
source=fullfile(folder,'synthetic_levels.csv'); fid=fopen(source,'w');
guard=onCleanup(@()fclose(fid));
for k=1:3, fprintf(fid,'%s\n',header(k)); end
for k=1:numel(dates)
    fprintf(fid,'%s',string(dates(k),'dd/MM/yyyy'));
    fprintf(fid,',%.16g',levels(k,:)); fprintf(fid,'\n');
end
clear guard
variables=["NKYTR Index","SPXT Index","SX5T Index","TUKXG Index"];
prior=run_weekly_model_report(SourceFile=source,VariableNames=variables, ...
    ParameterSource="prior",PriorParameterDraws=2, ...
    WarmupIterations=1,RetainedDraws=2,MaxHours=.1, ...
    CorrelationBackend="mex_new",CorrelationThreads=1,RunTests=false);
verifyEqual(testCase,prior.inference_mode,'prior_parameter_smoothing');
verifyEqual(testCase,prior.saved_draws,4);
verifyEqual(testCase,prior.saved_draws_per_chain,[2 2]);
verifyTrue(testCase,prior.prior_sampling.mixture_complete);
verifyFalse(testCase,prior.inference.estimate_parameters_requested);
verifyTrue(testCase,prior.inference.parameter_uncertainty_in_bands);
verifyEqual(testCase,prior.inference.parameter_source,'prior');
verifyEmpty(testCase,prior.paths.parameters);
verifyTrue(testCase,isfile(prior.paths.prior_parameters));
verifyTrue(testCase,all(isfile(prior.figure_paths)));
verifyTrue(testCase,isfile(fullfile(prior.run_dir,'inference_metadata.json')));
verifyFalse(testCase,prior.convergence_established);
verifyEmpty(testCase,dir(fullfile(cfg.output_root,'parameters','parameters_*.mat')));
estimated=run_weekly_model_report(SourceFile=source,VariableNames=variables, ...
    SmoothStates=false,WarmupIterations=1,RetainedDraws=2,MaxHours=.1, ...
    CorrelationBackend="mex_new",CorrelationThreads=1,RunTests=false);
verifyEqual(testCase,estimated.inference_mode,'parameter_estimation');
verifyEqual(testCase,estimated.inference.workflow,'estimate_only');
verifyFalse(testCase,estimated.inference.smoothing_stage_performed);
verifyEmpty(testCase,estimated.figure_paths);
verifyFalse(testCase,isfile(fullfile(estimated.run_dir,'convergence','chain_mean_paths.pdf')));
verifyTrue(testCase,isfile(estimated.paths.parameters));
verifyTrue(testCase,isfile(fullfile(estimated.run_dir,'result.mat')));
verifyTrue(testCase,isfile(estimated.convergence_artifacts.json));
end

function cfg=weekly_config_from_root(root)
% The copied entry point overrides the data path, but the CSV header is read
% from the actual repository, which is deliberately not copied into the test.
cfg=struct('source_file',fullfile(root,'data','yad_tickers_no_usdcnh_from_20030804.csv'));
end
