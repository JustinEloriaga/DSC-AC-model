function results = benchmark_correlation_backends(output_file)
%BENCHMARK_CORRELATION_BACKENDS Time the same 100 correlation draws in MEX.
% No estimation, file output or directory creation unless output_file is set.
% Timings include five-output transforms and masked Gaussian likelihoods.
% An optional dsc_corr_batch_mex_prototype on the path adds the prior build.
if nargin<1, output_file=''; end
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'toolbox'));
if strlength(string(output_file))>0
    output_parent=fileparts(char(output_file));
    if isempty(output_parent), output_parent=pwd; end
    assert(isfolder(output_parent),'Benchmark output directory must already exist.');
end
saved_rng=rng; restore_rng=onCleanup(@()rng(saved_rng)); %#ok<NASGU>
rng(20261012,'twister');
T=100; m=14; q=m*(m-1)/2; selector=tril(true(m),-1);
r=zeros(T,q);
for t=1:T
    X=randn(m); S=X*X'/m+.5*eye(m); sd=sqrt(diag(S));
    A=logm(S./(sd*sd')); r(t,:)=A(selector)';
end
Z=randn(T,m); mask=true(T,m);
mask(10:10:T,1:3)=false; mask(25:25:T,:)=false; Z(~mask)=NaN;
names={'mex','mex_new'};
kernels={@dsc_corr_batch_mex,@dsc_corr_batch_mex_new};
if exist('dsc_corr_batch_mex_prototype','file')==3
    names{end+1}='prototype'; kernels{end+1}=@dsc_corr_batch_mex_prototype;
end
results=struct('draws',T,'dimension',m,'seed',20261012, ...
    'matlab',version,'architecture',computer('arch'), ...
    'description','Identical synthetic correlation matrices; masked likelihood; five outputs; median of seven alternating trials.', ...
    'measurements',[]);
fprintf('100 correlation draws, m=14; seconds per complete five-output call\n');
fprintf('threads backend seconds speedup_vs_mex max_abs_P abs_ll_diff median_iter max_iter\n');
for threads=[1 4]
    count=numel(kernels); outputs=cell(count,5); calls=zeros(1,count);
    % Warm each kernel and measure enough calls to avoid short-timer noise.
    for k=1:count
        f=kernels{k}; [outputs{k,:}]=f(r,Z,mask,threads);
        started=tic; [outputs{k,:}]=f(r,Z,mask,threads); elapsed=toc(started);
        calls(k)=min(200,max(1,ceil(.15/max(elapsed,1e-6))));
    end
    times=zeros(7,count);
    for trial=1:7
        order=1:count; if mod(trial,2)==0, order=fliplr(order); end
        for k=order
            f=kernels{k}; started=tic;
            for rep=1:calls(k), [outputs{k,:}]=f(r,Z,mask,threads); end
            times(trial,k)=toc(started)/calls(k);
        end
    end
    medians=median(times,1);
    for k=1:count
        diffP=outputs{k,2}-outputs{1,2};
        deltaP=max(abs(diffP(:))); deltaLL=abs(outputs{k,1}-outputs{1,1});
        assert(deltaP<2e-7,'Correlation matrices disagree.');
        assert(deltaLL<2e-5,'Likelihoods disagree.');
        row=struct('threads',threads,'backend',names{k}, ...
            'seconds_per_100_draws',medians(k),'speedup_vs_mex',medians(1)/medians(k), ...
            'max_abs_P_difference',deltaP,'abs_likelihood_difference',deltaLL, ...
            'median_iterations',median(outputs{k,5}),'max_iterations',max(outputs{k,5}), ...
            'trials_seconds',times(:,k)','calls_per_trial',calls(k));
        results.measurements=[results.measurements; row]; %#ok<AGROW>
        fprintf('%d %-10s %.6f %.2f %.3g %.3g %.1f %.0f\n', ...
            threads,names{k},medians(k),row.speedup_vs_mex,deltaP,deltaLL, ...
            row.median_iterations,row.max_iterations);
    end
end
if strlength(string(output_file))>0
    fid=fopen(char(output_file),'w'); assert(fid~=-1,'Could not open benchmark output.');
    close_file=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s\n',jsonencode(results,'PrettyPrint',true));
end
end
