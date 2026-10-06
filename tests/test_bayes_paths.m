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

function testDiagnosticMetadataPreservedWithoutExtraJson(testCase)
folder=testCase.TestData.folder; metadata=testCase.TestData.metadata;
metadata.estimation_convergence=struct('status','not_checked','passed',false, ...
    'max_rhat',NaN,'chunk_metadata',{{struct('file','one_chunk.mat')}});
write_fixture(fullfile(folder,'inference_metadata.json'),metadata);
write_fixture(fullfile(folder,'run_manifest.json'),struct('source_hash','fixture','inference',metadata));
plot_bayes_correlation_paths(folder);
original=jsondecode(fileread(fullfile(folder,'inference_metadata.json')));
bands=load(fullfile(folder,'posterior_correlation_bands.mat'),'inference');
verifyEqual(testCase,bands.inference,original);
verifyFalse(testCase,isfile(fullfile(folder,'figures','bayes_correlation_paths_manifest.json')));
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
verifyEqual(testCase,numel(paths),1);
verifyTrue(testCase,all(isfile(paths)));
bands=load(fullfile(testCase.TestData.folder,'posterior_correlation_bands.mat'));
verifyEqual(testCase,bands.inference,metadata);
verifyEqual(testCase,char(string(bands.dates(1),'yyyy-MM-dd')),'2020-01-03');
verifyEqual(testCase,char(string(bands.dates(end),'yyyy-MM-dd')),'2020-04-17');
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
bands=load(fullfile(testCase.TestData.folder,'posterior_correlation_bands.mat'),'inference');
verifyEmpty(testCase,bands.inference);
end

function testMultipleChainsPoolBandsAndPreserveDrawProvenance(testCase)
[dirs,expected]=multichain_fixture(testCase);
paths=plot_bayes_correlation_paths(testCase.TestData.folder);
verifyTrue(testCase,all(isfile(paths)));
bands=load(fullfile(testCase.TestData.folder,'posterior_correlation_bands.mat'));
quantiles=prctile(expected,[16 50 84],3);
verifyEqual(testCase,bands.lower,quantiles(:,:,1));
verifyEqual(testCase,bands.median_path,quantiles(:,:,2));
verifyEqual(testCase,bands.upper,quantiles(:,:,3));
verifyEqual(testCase,bands.chain_dirs,dirs);
verifyEqual(testCase,bands.chain_ids,[3 7]);
verifyEqual(testCase,bands.draw_chain_ids,[3 3 3 7 7]);
verifyEqual(testCase,bands.iterations,[2 3 4 2 3]);
verifyEqual(testCase,bands.retained_draws_per_chain,[3 2]);
verifyEqual(testCase,bands.warmup_per_chain,[1 1]);
verifyFalse(testCase,isfile(fullfile(testCase.TestData.folder,'figures','bayes_correlation_paths_manifest.json')));
end

function testDuplicateChainDirectoriesRejectBeforeExport(testCase)
[dirs,~]=multichain_fixture(testCase);
summary=jsondecode(fileread(fullfile(testCase.TestData.folder,'pilot_summary.json')));
summary.chain_dirs={dirs{1},dirs{1}};
write_fixture(fullfile(testCase.TestData.folder,'pilot_summary.json'),summary);
verifyError(testCase,@()plot_bayes_correlation_paths(testCase.TestData.folder),'DSC:PathIdentity');
verifyFalse(testCase,isfolder(fullfile(testCase.TestData.folder,'figures')));
end

function testDuplicateChainIdsRejectBeforeExport(testCase)
[dirs,~]=multichain_fixture(testCase);
path=fullfile(dirs{2},'posterior_chunk_000001.mat'); loaded=load(path,'chunk');
chunk=loaded.chunk; chunk.chain_id=3; save(path,'chunk','-v7');
verifyError(testCase,@()plot_bayes_correlation_paths(testCase.TestData.folder),'DSC:PathIdentity');
end

function testPerChainCountsMustMatchSummary(testCase)
multichain_fixture(testCase);
summary=jsondecode(fileread(fullfile(testCase.TestData.folder,'pilot_summary.json')));
summary.saved_draws_per_chain=[2 3];
write_fixture(fullfile(testCase.TestData.folder,'pilot_summary.json'),summary);
verifyError(testCase,@()plot_bayes_correlation_paths(testCase.TestData.folder),'DSC:PathDraws');
end

function testIterationsMustIncreaseWithinEachChain(testCase)
[dirs,~]=multichain_fixture(testCase);
path=fullfile(dirs{2},'posterior_chunk_000001.mat'); loaded=load(path,'chunk');
chunk=loaded.chunk; chunk.iterations=[3 2]; save(path,'chunk','-v7');
verifyError(testCase,@()plot_bayes_correlation_paths(testCase.TestData.folder),'DSC:PathDraws');
end

function testMultichainDateIdentityMustMatch(testCase)
[dirs,~]=multichain_fixture(testCase);
path=fullfile(dirs{2},'posterior_chunk_000001.mat'); loaded=load(path,'chunk');
chunk=loaded.chunk; chunk.dates=chunk.dates+calweeks(1); save(path,'chunk','-v7');
verifyError(testCase,@()plot_bayes_correlation_paths(testCase.TestData.folder),'DSC:PathIdentity');
end

function testCannotLabelConvergedWithoutSavedDiagnostics(testCase)
multichain_fixture(testCase);
path=fullfile(testCase.TestData.folder,'pilot_summary.json');
summary=jsondecode(fileread(path)); summary.convergence_established=true;
write_fixture(path,summary);
verifyError(testCase,@()plot_bayes_correlation_paths(testCase.TestData.folder),'DSC:ConvergenceMetadata');
verifyFalse(testCase,isfolder(fullfile(testCase.TestData.folder,'figures')));
end

function [dirs,expected]=multichain_fixture(testCase)
folder=testCase.TestData.folder;
dirs={fullfile(folder,'chains','chain_001'),fullfile(folder,'chains','chain_002')};
for c=1:2, mkdir(dirs{c}); end
loaded=load(fullfile(folder,'posterior_chunk_000001.mat'),'chunk');
chunk=loaded.chunk; chunk.P_pairs=chunk.P_pairs(:,1,:);
chunk.tickers=["A","B"]; chunk.pair_i=2; chunk.pair_j=1;
chunk.chain_id=3; expected=chunk.P_pairs;
save(fullfile(dirs{1},'posterior_chunk_000001.mat'),'chunk','-v7');
chunk.chain_id=7; chunk.iterations=[2 3]; chunk.P_pairs=chunk.P_pairs(:,:,1:2)+.1;
expected=cat(3,expected,chunk.P_pairs);
save(fullfile(dirs{2},'posterior_chunk_000001.mat'),'chunk','-v7');
% The coordinator summary is authoritative; root chunks belong to the legacy
% fixture and must not accidentally be included in this pooled run.
summary=testCase.TestData.summary;
summary.chain_dirs=dirs; summary.chain_ids=[3 7]; summary.saved_draws=5;
summary.saved_draws_per_chain=[3 2]; summary.num_chains=2;
write_fixture(fullfile(folder,'pilot_summary.json'),summary);
end

function write_fixture(path,value)
fid=fopen(path,'w'); assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value));
end
