function report=benchmark_dsc_kernel(checkpoint_file,output_dir)
% Fixed-input likelihood benchmark; does not advance any MCMC chain.
root=fileparts(mfilename('fullpath')); addpath(fullfile(root,'toolbox'));
if nargin<1||isempty(checkpoint_file)
    checkpoint_file=fullfile(root,'outputs','weekly_research','runs', ...
        '20260911-175022-407-chain1','checkpoint.mat');
end
if nargin<2, output_dir=fullfile(root,'outputs','weekly_research','performance'); end
if ~exist(output_dir,'dir'), mkdir(output_dir); end
loaded=load(checkpoint_file,'checkpoint'); cp=loaded.checkpoint;
r=cp.state.r; Z=(cp.identity.returns-cp.state.B).*exp(-cp.state.h/2);
mask=cp.identity.mask; Z(~mask)=0;
assert(exist('dsc_corr_batch_mex','file')==3,'dsc:MissingKernel','Run build_dsc_mex first.');
% Warm up both paths before wall-clock measurements.
reference=matlab_likelihood(r,Z,mask);
[native,P,L,logdet,iterations]=dsc_corr_batch_mex(r,Z,mask,4); %#ok<ASGLU>
assert(abs(reference-native)<1e-6,'dsc:BenchmarkMismatch','Native/reference likelihood mismatch.');
thread_counts=[1 2 4]; seconds_matlab=zeros(5,1); seconds_by_threads=zeros(5,numel(thread_counts));
for k=1:5
    timer=tic; matlab_likelihood(r,Z,mask); seconds_matlab(k)=toc(timer);
    for j=1:numel(thread_counts)
        timer=tic; value=dsc_corr_batch_mex(r,Z,mask,thread_counts(j));
        seconds_by_threads(k,j)=toc(timer);
        assert(abs(value-reference)<1e-6,'dsc:BenchmarkMismatch','Threaded/reference likelihood mismatch.');
    end
end
seconds_native=seconds_by_threads(:,end);
profile clear; profile on
matlab_likelihood(r,Z,mask);
profile off; baseline_profile=profile('info');
save(fullfile(output_dir,'matlab_likelihood_profile.mat'),'baseline_profile');
f=baseline_profile.FunctionTable; [~,order]=sort([f.TotalTime],'descend');
hotspots=arrayfun(@(x)struct('function',x.FunctionName,'seconds',x.TotalTime, ...
    'calls',x.NumCalls),f(order(1:min(12,numel(order)))));
report=struct('checkpoint',checkpoint_file,'source_hash',cp.identity.source_hash, ...
    'baseline_sweep_seconds',cp.diagnostics.sweep_seconds, ...
    'baseline_r_proposals_by_sweep',cp.diagnostics.r_proposals, ...
    'baseline_h_proposals_by_sweep',cp.diagnostics.h_proposals, ...
    'T',size(Z,1),'m',size(Z,2),'pairs',size(r,2), ...
    'reference_loglik_without_constants',reference,'native_loglik_without_constants',native, ...
    'absolute_likelihood_error',abs(reference-native), ...
    'matlab_seconds',seconds_matlab,'native_seconds',seconds_native, ...
    'thread_counts',thread_counts,'native_seconds_by_thread_count',seconds_by_threads, ...
    'median_matlab_seconds',median(seconds_matlab),'median_native_seconds',median(seconds_native), ...
    'kernel_speedup',median(seconds_matlab)/median(seconds_native), ...
    'mean_fixed_point_iterations',mean(iterations),'max_fixed_point_iterations',max(iterations), ...
    'implementation',dsc_sampler_signature('mex'),'matlab_profile_hotspots',hotspots, ...
    'scope','Same saved warm-up state and likelihood; no MCMC, no convergence claim.');
fid=fopen(fullfile(output_dir,'kernel_benchmark.json'),'w');
assert(fid>=0); guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('Kernel benchmark: MATLAB %.4fs, native %.4fs, %.2fx; likelihood error %.3g\n', ...
    report.median_matlab_seconds,report.median_native_seconds,report.kernel_speedup,report.absolute_likelihood_error);
end

function ll=matlab_likelihood(r,Z,mask)
ll=0; x=[];
for t=1:size(Z,1)
    [C,~,x]=dsc_corr_matrix(r(t,:),x); ix=find(mask(t,:));
    if isempty(ix), continue; end
    L=chol(C(ix,ix),'lower'); u=L\Z(t,ix)';
    ll=ll-.5*(2*sum(log(diag(L)))+u'*u);
end
end
