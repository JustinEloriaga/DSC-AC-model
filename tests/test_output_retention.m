function tests=test_output_retention
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
root=fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.previous_path=path;
addpath(fullfile(root,'toolbox'));
end

function teardownOnce(testCase)
path(testCase.TestData.previous_path);
end

function testKeepsThreeNewestRunsAndThreeOfEachPdf(testCase)
root=tempname;
mkdir(root);
testCase.addTeardown(@()rmdir(root,'s'));
runs=fullfile(root,'outputs','runs'); mkdir(runs);
reports=fullfile(root,'reports'); mkdir(reports);
ids=["20260101-000000-chain1-w10-r100-thin1", ...
    "20260102-000000-chain1-w10-r100-thin1", ...
    "20260103-000000-chain1-w10-r100-thin1", ...
    "20260104-000000-chain1-w10-r100-thin1"];
for id=ids, mkdir(fullfile(runs,id)); end
for day=1:4
    touch(fullfile(reports,sprintf('report_2026010%d-000000-chain1-w10-r100-thin1.pdf',day)));
    touch(fullfile(reports,sprintf('figure_2026010%d-000000-chain1-w10-r100-thin1.pdf',day)));
end
current_run=fullfile(runs,ids(end));
dsc_prune_outputs(fullfile(root,'outputs'),reports,current_run);
entries=dir(runs);
remaining=string({entries([entries.isdir]).name});
remaining=remaining(~ismember(remaining,[".",".."]));
verifyEqual(testCase,sort(remaining),sort(ids(end-2:end)));
verifyEqual(testCase,numel(dir(fullfile(reports,'report_*.pdf'))),3);
verifyEqual(testCase,numel(dir(fullfile(reports,'figure_*.pdf'))),3);
end

function touch(file)
fid=fopen(file,'w');
assert(fid>=0);
fclose(fid);
end
