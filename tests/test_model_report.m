function tests=test_model_report
tests=functiontests(localfunctions);
end

function testCsvToEstimateAndLoadWorkflow(testCase)
% Exercise the public options on a small synthetic 14-variable CSV. No
% real-market estimation or shared publication output is generated here.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'toolbox'),fullfile(root,'scripts'));
prior_warning=warning('off','econ:adftest:InvalidStatistic');
warning_guard=onCleanup(@()warning(prior_warning)); %#ok<NASGU>
folder=tempname; mkdir(folder);
testCase.addTeardown(@()rmdir(folder,'s'));
cfg=weekly_config(); header=readlines(cfg.source_file);
% Run the fixed-output entry point from an isolated repository copy.
prior_path=path;
path_guard=onCleanup(@()path(prior_path)); %#ok<NASGU>
for name={'run_weekly_model_report.m','weekly_config.m','run_weekly_research.m'}
    copyfile(fullfile(root,name{1}),fullfile(folder,name{1}));
end
copyfile(fullfile(root,'toolbox'),fullfile(folder,'toolbox'));
copyfile(fullfile(root,'scripts'),fullfile(folder,'scripts'));
addpath(folder,fullfile(folder,'toolbox'),fullfile(folder,'scripts'));
output_folder=fullfile(folder,'outputs');
dates=(datetime(2003,8,4):days(1):datetime(2005,10,7))';
dates=dates(~ismember(weekday(dates),[1 7]));
rng(139,'twister');
levels=100*exp(cumsum(.002*randn(numel(dates),14),1));
source=fullfile(folder,'synthetic_levels.csv'); fid=fopen(source,'w');
guard=onCleanup(@()fclose(fid));
for i=1:3, fprintf(fid,'%s\n',header(i)); end
for i=1:numel(dates)
    fprintf(fid,'%s',string(dates(i),'dd/MM/yyyy'));
    fprintf(fid,',%.16g',levels(i,:)); fprintf(fid,'\n');
end
clear guard
result=run_weekly_model_report(SourceFile=source, ...
    EstimateParameters=true,EstimationEndDate="2005-08-26", ...
    WarmupIterations=1,RetainedDraws=2,MaxHours=.1, ...
    CorrelationThreads=1,RunTests=false,BuildReport=false,VerifyReport=false);
verifyEqual(testCase,result.saved_draws,2);
[~,run_name]=fileparts(result.run_dir);
verifyTrue(testCase,~isempty(regexp(run_name, ...
    '^\d{8}-\d{6}-chain1-w1-r2-thin1$','once')));
verifyEqual(testCase,result.inference_mode,'fixed_parameter_smoothing');
verifyEqual(testCase,result.inference.parameter_smoothing,'mean');
verifyFalse(testCase,result.inference.parameter_uncertainty_in_bands);
verifyEqual(testCase,result.T,9);
verifyEqual(testCase,result.dates,(datetime(2005,8,12)+calweeks(0:8))');
verifyEqual(testCase,result.inference.calibration_weeks,104);
verifyEqual(testCase,result.inference.calibration_start,'2003-08-15');
verifyEqual(testCase,result.inference.calibration_end,'2005-08-05');
verifyEqual(testCase,result.inference.parameter_estimation_start,'2005-08-12');
verifyEqual(testCase,result.inference.parameter_estimation_end,'2005-08-26');
verifyEqual(testCase,result.inference.smoothing_start,'2005-08-12');
verifyEqual(testCase,result.inference.smoothing_end,'2005-10-07');
verifyTrue(testCase,isfile(result.paths.parameters));
verifyTrue(testCase,all(isfile(result.figure_paths)));
verifyFalse(testCase,isfile(fullfile(result.run_dir,'prior_summary.json')));
verifyFalse(testCase,isfile(fullfile(result.run_dir,'result.mat')));
estimated=load(result.paths.parameters,'parameters');
verifyEqual(testCase,size(estimated.parameters.training_returns),[3 14]);
verifyEqual(testCase,estimated.parameters.training_dates, ...
    (datetime(2005,8,12)+calweeks(0:2))');
smoothed=load(fullfile(result.run_dir,'posterior_chunk_000001.mat'),'chunk');
verifyEqual(testCase,size(smoothed.chunk.h,1),9);
verifyEqual(testCase,smoothed.chunk.dates,result.dates);
verifyTrue(testCase,all(estimated.parameters.training_mask,'all'));
verifyEqual(testCase,fileparts(fileparts(result.run_dir)),output_folder);
verifyFalse(testCase,isfolder(fullfile(output_folder,'data')));
loaded=run_weekly_model_report(SourceFile=source, ...
    EstimateParameters=false,WarmupIterations=1,RetainedDraws=2,MaxHours=.1, ...
    CorrelationThreads=1,RunTests=false,BuildReport=false,VerifyReport=false);
verifyEqual(testCase,loaded.paths.parameters,result.paths.parameters);
verifyEqual(testCase,loaded.inference.parameters_estimated_at,result.inference.parameters_estimated_at);
verifyEqual(testCase,loaded.T,9);
verifyEqual(testCase,loaded.dates,result.dates);
verifyEqual(testCase,loaded.inference.calibration_weeks,104);
verifyEqual(testCase,loaded.inference.parameter_smoothing,'mean');
verifyFalse(testCase,loaded.inference.parameter_uncertainty_in_bands);
cp=load(fullfile(loaded.run_dir,'posterior_chunk_000001.mat'),'chunk');
verifyEqual(testCase,cp.chunk.dates,result.dates);
verifyEqual(testCase,size(cp.chunk.h,1),9);
verifyEqual(testCase,numel(dir(fullfile(output_folder,'parameters','parameters_*.mat'))),1);
latest=jsondecode(fileread(fullfile(output_folder,'latest_model_report.json')));
verifyEqual(testCase,latest.inference.parameter_file,result.paths.parameters);
verifyEqual(testCase,latest.actual_retained_draws_per_chain,2);
verifyEqual(testCase,latest.run_result.saved_draws_total,2);
verifyEqual(testCase,latest.run_result.status,loaded.status);
end
