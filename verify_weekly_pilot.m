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
assert(T==1204 && m==14);
assert(isequal(size(cp.state.B),[T,m]) && isequal(size(cp.state.h),[T,m]));
assert(isequal(size(cp.state.r),[T,91]));
assert(isequal(cp.identity.dates([1 end]),[datetime(2003,8,15);datetime(2026,9,4)]));
assert(~any(contains(cp.identity.tickers,'USDCNH')));
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
