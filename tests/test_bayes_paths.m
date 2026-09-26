function tests=test_bayes_paths
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'scripts'));
testCase.TestData.root=root;
end

function setup(testCase)
folder=tempname; mkdir(folder);
testCase.addTeardown(@()rmdir(folder,'s'));
dates=(datetime(2020,1,3)+calweeks(0:15))';
[ii,jj]=find(tril(true(3),-1));
paths=[linspace(-.2,.4,16)' linspace(.4,.1,16)' .1*sin((1:16)'/3)];
chunk=struct('dates',dates,'tickers',["A","B","C"],'pair_i',ii,'pair_j',jj, ...
    'chain_id',1,'P_pairs',cat(3,paths-.02,paths,paths+.02),'iterations',2:4);
save(fullfile(folder,'posterior_chunk_000001.mat'),'chunk','-v7');
summary=struct('saved_draws',3,'burnin',1,'convergence_established',false, ...
    'first_date','2020-01-03','last_date','2020-04-17','inference_mode','fixed_parameter_smoothing');
metadata=struct('inference_mode','fixed_parameter_smoothing','parameter_file','fixture_parameters.mat', ...
    'parameters_estimated_at','2026-09-26T11:22:33.000Z','parameter_estimation_start','2020-01-03', ...
    'parameter_estimation_end','2020-02-21','parameter_estimation_draws',5, ...
    'parameter_estimation_run','fixture_estimation','smoothing_end','2020-04-17', ...
    'parameter_smoothing','mean','parameter_uncertainty_in_bands',false);
write_fixture(fullfile(folder,'pilot_summary.json'),summary);
write_fixture(fullfile(folder,'inference_metadata.json'),metadata);
write_fixture(fullfile(folder,'run_manifest.json'),struct('source_hash','fixture','inference',metadata));
testCase.TestData.folder=folder; testCase.TestData.metadata=metadata;
testCase.TestData.summary=summary;
end

function testConditionalPathsRecordSavedProvenance(testCase)
metadata=testCase.TestData.metadata;
metadata.calibration_weeks=104;
metadata.calibration_start='2018-01-05'; metadata.calibration_end='2019-12-27';
metadata.full_data_start='2018-01-05'; metadata.smoothing_start='2020-01-03';
summary=testCase.TestData.summary; summary.excluded_initial_weeks=104;
write_fixture(fullfile(testCase.TestData.folder,'pilot_summary.json'),summary);
write_fixture(fullfile(testCase.TestData.folder,'inference_metadata.json'),metadata);
write_fixture(fullfile(testCase.TestData.folder,'run_manifest.json'), ...
    struct('source_hash','fixture','inference',metadata));
paths=plot_bayes_correlation_paths(testCase.TestData.folder);
verifyEqual(testCase,numel(paths),2);
verifyTrue(testCase,all(isfile(paths)));
manifest=jsondecode(fileread(fullfile(testCase.TestData.folder,'figures','bayes_correlation_paths_manifest.json')));
verifyEqual(testCase,manifest.inference,metadata);
verifyEqual(testCase,manifest.first_date,'2020-01-03');
verifyEqual(testCase,manifest.last_date,'2020-04-17');
bands=load(fullfile(testCase.TestData.folder,'posterior_correlation_bands.mat'),'inference');
verifyEqual(testCase,bands.inference,metadata);
% Optional local visual QA preserves a copy before fixture teardown.
qa_dir=getenv('DSC_TEST_PDF_QA_DIR');
if ~isempty(qa_dir)
    if ~isfolder(qa_dir), mkdir(qa_dir); end
    copyfile(paths(1),fullfile(qa_dir,'correlation_paths_calibration.pdf'));
end
end

function testStaleRunMetadataRejectsBeforeExport(testCase)
meta=testCase.TestData.metadata; meta.parameter_estimation_end='2020-02-14';
write_fixture(fullfile(testCase.TestData.folder,'inference_metadata.json'),meta);
verifyError(testCase,@()plot_bayes_correlation_paths(testCase.TestData.folder),'DSC:InferenceMetadata');
verifyFalse(testCase,isfolder(fullfile(testCase.TestData.folder,'figures')));
end

function testSampleMismatchRejectsBeforeExport(testCase)
meta=testCase.TestData.metadata; meta.smoothing_end='2020-04-10';
write_fixture(fullfile(testCase.TestData.folder,'inference_metadata.json'),meta);
write_fixture(fullfile(testCase.TestData.folder,'run_manifest.json'),struct('source_hash','fixture','inference',meta));
verifyError(testCase,@()plot_bayes_correlation_paths(testCase.TestData.folder),'DSC:InferenceMetadata');
end

function testMissingFixedMetadataRejects(testCase)
delete(fullfile(testCase.TestData.folder,'inference_metadata.json'));
verifyError(testCase,@()plot_bayes_correlation_paths(testCase.TestData.folder),'DSC:InferenceMetadata');
end

function testLegacyPathsRemainSupported(testCase)
delete(fullfile(testCase.TestData.folder,'inference_metadata.json'));
write_fixture(fullfile(testCase.TestData.folder,'run_manifest.json'), ...
    struct('source_hash','fixture','inference','historical full-sample smoothing; empirical-Bayes priors'));
summary=rmfield(testCase.TestData.summary,'inference_mode');
write_fixture(fullfile(testCase.TestData.folder,'pilot_summary.json'),summary);
paths=plot_bayes_correlation_paths(testCase.TestData.folder);
verifyTrue(testCase,all(isfile(paths)));
manifest=jsondecode(fileread(fullfile(testCase.TestData.folder,'figures','bayes_correlation_paths_manifest.json')));
verifyFalse(testCase,isfield(manifest,'inference'));
end

function write_fixture(path,value)
fid=fopen(path,'w'); assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value));
end
