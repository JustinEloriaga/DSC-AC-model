function report = verify_weekly_pilot(run_dir)
%VERIFY_WEEKLY_PILOT Check a completed bounded pilot without advancing its chain.
% Reconstructs the last warm-up state's matrices only for numerical checks;
% they are not posterior estimates and are not published as such.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'toolbox'));
loaded=load(fullfile(run_dir,'checkpoint.mat'),'checkpoint');
cp=loaded.checkpoint;
summary=jsondecode(fileread(fullfile(run_dir,'pilot_summary.json')));
assert(ismember(string(summary.status),["time_limit","iteration_limit"]), ...
    'DSC:PilotFailure','Pilot ended with a numerical failure.');
assert(cp.completed==summary.completed_iterations && cp.completed<=20);
assert(cp.saved==0 && summary.saved_draws==0 && summary.burnin>=cp.completed);
assert(~summary.convergence_established);
[T,m]=size(cp.identity.returns);
nr=m*(m-1)/2;
assert(T>=1 && m>=2 && summary.T==T && summary.m==m && summary.pairs==nr, ...
    'DSC:PilotDimensions','Sampler summary and checkpoint dimensions differ.');
dates=cp.identity.dates;
assert(isdatetime(dates) && numel(dates)==T && ~any(isnat(dates)) && ...
    all(diff(dates)>seconds(0)) && numel(cp.identity.tickers)==m, ...
    'DSC:PilotDates','Checkpoint dates and series must match the observation arrays.');
first_date=char(string(dates(1),'yyyy-MM-dd'));
last_date=char(string(dates(end),'yyyy-MM-dd'));
assert(strcmp(summary.first_date,first_date) && strcmp(summary.last_date,last_date), ...
    'DSC:PilotDates','Summary dates differ from the checkpoint likelihood dates.');
assert(isequal(size(cp.identity.mask),[T,m]), ...
    'DSC:PilotDimensions','Observation mask dimensions differ from the likelihood.');
assert(isequal(size(cp.state.B),[T,m]) && isequal(size(cp.state.h),[T,m]));
assert(isequal(size(cp.state.r),[T,nr]));
excluded=0; calibration_start=''; calibration_end='';
priors=cp.identity.priors;
exclusion_fields={'likelihood_start_date','excluded_initial_weeks','calibration_dates'};
if any(isfield(priors,exclusion_fields))
    assert(all(isfield(priors,exclusion_fields)), ...
        'DSC:PilotCalibration','Checkpoint has incomplete calibration-exclusion metadata.');
    excluded=priors.excluded_initial_weeks;
    calibration_dates=priors.calibration_dates;
    assert(isnumeric(excluded) && isscalar(excluded) && isfinite(excluded) && ...
        excluded>=1 && excluded==floor(excluded) && ...
        isdatetime(calibration_dates) && numel(calibration_dates)==2 && ...
        ~any(isnat(calibration_dates)) && calibration_dates(1)<=calibration_dates(2) && ...
        calibration_dates(2)<dates(1) && isequal(priors.likelihood_start_date,dates(1)), ...
        'DSC:PilotCalibration','Calibration dates overlap or disagree with the likelihood sample.');
    if isfield(priors,'initial_weeks')
        assert(priors.initial_weeks==excluded,'DSC:PilotCalibration', ...
            'Calibration exclusion count differs from the initial-prior window.');
    end
    if isfield(priors,'initial_dates')
        assert(isequal(priors.initial_dates(:),calibration_dates(:)), ...
            'DSC:PilotCalibration','Calibration dates differ from the initial-prior window.');
    end
    calibration_start=char(string(calibration_dates(1),'yyyy-MM-dd'));
    calibration_end=char(string(calibration_dates(2),'yyyy-MM-dd'));
end
assert(all(isfinite(cp.diagnostics.loglik)));
assert(all(cp.state.sig2h>0) && all(cp.state.sig2r>0));
[~,bad]=chol(cp.state.V); assert(bad==0);
minimum_eigenvalue=Inf; maximum_symmetry_error=0;
maximum_diagonal_error=0; maximum_coordinate_error=0;
sel=tril(true(m),-1);
for t=1:T
    C=dsc_corr_matrix(cp.state.r(t,:)');
    [Q,D]=eig(C,'vector');
    minimum_eigenvalue=min(minimum_eigenvalue,min(D));
    maximum_symmetry_error=max(maximum_symmetry_error,max(abs(C-C'),[],'all'));
    maximum_diagonal_error=max(maximum_diagonal_error,max(abs(diag(C)-1)));
    A=Q*diag(log(D))*Q';
    maximum_coordinate_error=max(maximum_coordinate_error,max(abs(A(sel)-cp.state.r(t,:)')));
end
assert(minimum_eigenvalue>0 && maximum_symmetry_error<1e-12);
assert(maximum_diagonal_error<1e-12 && maximum_coordinate_error<1e-6);
report=struct('status','passed','purpose','Numerical validation of last completed warm-up state; not posterior inference', ...
    'completed_iterations',cp.completed,'saved_posterior_draws',cp.saved, ...
    'matrices_checked',T,'variables',m,'distinct_pairs',nnz(sel), ...
    'first_date',first_date,'last_date',last_date, ...
    'excluded_initial_weeks',excluded,'calibration_start',calibration_start, ...
    'calibration_end',calibration_end, ...
    'minimum_eigenvalue',minimum_eigenvalue, ...
    'maximum_symmetry_error',maximum_symmetry_error, ...
    'maximum_diagonal_error',maximum_diagonal_error, ...
    'maximum_coordinate_roundtrip_error',maximum_coordinate_error, ...
    'source_hash',cp.identity.source_hash);
fid=fopen(fullfile(run_dir,'validation.json'),'w');
assert(fid>=0,'DSC:WriteFailed','Cannot write pilot validation evidence.');
guard=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
end
