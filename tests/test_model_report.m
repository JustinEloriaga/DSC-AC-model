function tests=test_model_report
tests=functiontests(localfunctions);
end

function testCsvToEstimateAndLoadWorkflow(testCase)
% Exercise the public options on a small synthetic 14-variable CSV. No
% real-market estimation or shared publication output is generated here.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root,fullfile(root,'toolbox'),fullfile(root,'scripts'));
prior_warning=warning('off','econ:adftest:InvalidStatistic');
warning_guard=onCleanup(@()warning(prior_warning)); %#ok<NASGU>
folder=tempname; mkdir(folder);
cfg=weekly_config(); header=readlines(cfg.source_file);
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
result=run_weekly_model_report(SourceFile=source,OutputRoot=folder, ...
    EstimateParameters=true,EstimationEndDate="2005-08-26", ...
    WarmupIterations=1,RetainedDraws=2,MaxHours=.1, ...
    RunTests=false,BuildReport=false,VerifyReport=false);
verifyEqual(testCase,result.saved_draws,2);
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
verifyTrue(testCase,isfile(fullfile(result.run_dir,'prior_summary.json')));
estimated=load(fullfile(result.run_dir,'estimation','checkpoint.mat'),'checkpoint');
verifyEqual(testCase,size(estimated.checkpoint.identity.returns),[3 14]);
verifyEqual(testCase,estimated.checkpoint.identity.dates, ...
    (datetime(2005,8,12)+calweeks(0:2))');
smoothed=load(result.paths.checkpoint,'checkpoint');
verifyEqual(testCase,size(smoothed.checkpoint.identity.returns),[9 14]);
verifyEqual(testCase,smoothed.checkpoint.identity.dates,result.dates);
verifyTrue(testCase,all(estimated.checkpoint.identity.mask,'all'));
verifyTrue(testCase,all(smoothed.checkpoint.identity.mask,'all'));
loaded=run_weekly_model_report(SourceFile=source,OutputRoot=folder, ...
    EstimateParameters=false,WarmupIterations=1,RetainedDraws=2,MaxHours=.1, ...
    RunTests=false,BuildReport=false,VerifyReport=false);
verifyEqual(testCase,loaded.paths.parameters,result.paths.parameters);
verifyEqual(testCase,loaded.inference.parameters_estimated_at,result.inference.parameters_estimated_at);
verifyEqual(testCase,loaded.T,9);
verifyEqual(testCase,loaded.dates,result.dates);
verifyEqual(testCase,loaded.inference.calibration_weeks,104);
verifyEqual(testCase,loaded.inference.parameter_smoothing,'mean');
verifyFalse(testCase,loaded.inference.parameter_uncertainty_in_bands);
cp=load(loaded.paths.checkpoint,'checkpoint');
verifyEqual(testCase,cp.checkpoint.identity.dates,result.dates);
verifyEqual(testCase,size(cp.checkpoint.state.r,1),9);
verifyEqual(testCase,numel(dir(fullfile(folder,'parameters','dsc_parameters_*.mat'))),1);
latest=jsondecode(fileread(fullfile(folder,'latest_model_report.json')));
verifyEqual(testCase,latest.inference.parameter_file,result.paths.parameters);
end
