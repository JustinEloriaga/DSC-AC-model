function checks = run_weekly_acceptance(cfg)
%RUN_WEEKLY_ACCEPTANCE Run synthetic/data checks and retain machine-readable evidence.
% Does not run the full-panel MCMC pilot or change the source data.
if nargin<1 || isempty(cfg), cfg=weekly_config(); end
root=fileparts(mfilename('fullpath'));
checks=runtests(fullfile(root,'tests'));
report=struct('created_at',char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss')), ...
    'matlab_version',version,'passed',sum([checks.Passed]), ...
    'failed',sum([checks.Failed]),'incomplete',sum([checks.Incomplete]), ...
    'duration_seconds',sum([checks.Duration]),'tests',[]);
report.tests=arrayfun(@(x)struct('name',x.Name,'passed',x.Passed, ...
    'failed',x.Failed,'incomplete',x.Incomplete,'duration_seconds',x.Duration),checks);
if ~exist(cfg.output_root,'dir'), mkdir(cfg.output_root); end
fid=fopen(fullfile(cfg.output_root,'matlab_test_results.json'),'w');
assert(fid>=0,'DSC:WriteFailed','Cannot write acceptance-test results.');
guard=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
assertSuccess(checks);
end
