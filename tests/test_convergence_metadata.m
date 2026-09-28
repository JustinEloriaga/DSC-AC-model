function tests=test_convergence_metadata
tests=functiontests(localfunctions);
end

function testPassingReportRequiresFullCoverageAndIdentity(testCase)
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'toolbox'));
report=passing_report();
dsc_validate_convergence(report,1:4,1000*ones(1,4),72);
bad=report; bad.quantities_checked=1;
verifyError(testCase,@()dsc_validate_convergence(bad,1:4,1000*ones(1,4),72),'dsc:ConvergenceMetadata');
bad=report; bad.chain_ids=[1 2 3 5];
verifyError(testCase,@()dsc_validate_convergence(bad,1:4,1000*ones(1,4),72),'dsc:ConvergenceMetadata');
bad=report; bad.max_rhat=-1;
verifyError(testCase,@()dsc_validate_convergence(bad,1:4,1000*ones(1,4),72),'dsc:ConvergenceMetadata');
bad=report; bad.max_mcse_sd_ratio=-.01;
verifyError(testCase,@()dsc_validate_convergence(bad,1:4,1000*ones(1,4),72),'dsc:ConvergenceMetadata');
bad=report; bad.thresholds.min_draws=4;
verifyError(testCase,@()dsc_validate_convergence(bad,1:4,1000*ones(1,4),72),'dsc:ConvergenceMetadata');
bad=report; bad.thresholds.rhat_threshold=1.1;
verifyError(testCase,@()dsc_validate_convergence(bad,1:4,1000*ones(1,4),72),'dsc:ConvergenceMetadata');
bad=report; bad.failing_quantities={'B(1)[2020-01-03]'};
verifyError(testCase,@()dsc_validate_convergence(bad,1:4,1000*ones(1,4),72),'dsc:ConvergenceMetadata');
end

function report=passing_report()
report=struct('status','passed','passed',true,'chains',4,'chain_ids',1:4, ...
    'actual_draws_per_chain',1000*ones(1,4),'draws_per_chain',1000, ...
    'quantities_checked',72,'scope','All joint quantities at all dates', ...
    'max_rhat',1.005,'min_ess_bulk',900,'min_ess_tail',750,'max_mcse_sd_ratio',.04, ...
    'thresholds',struct('rhat',1.01,'ess',400,'mcse_ratio',.05,'min_draws',100), ...
    'failing_quantities',{{}},'unavailable_quantities',{{}});
end
