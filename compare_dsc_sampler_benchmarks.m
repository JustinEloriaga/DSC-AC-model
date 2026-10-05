function comparison=compare_dsc_sampler_benchmarks(baseline_run,optimized_run)
% Reconcile timings and every retained draw, not just the final MCMC state.
arguments
    baseline_run (1,1) string
    optimized_run (1,1) string
end
a=load(fullfile(baseline_run,'benchmark_report.mat'),'report'); baseline=a.report;
b=load(fullfile(optimized_run,'benchmark_report.mat'),'report'); optimized=b.report;
for field={'source_hash','tickers','dates','seed','warmup_iterations','retained_draws'}
    name=field{1}; assert(isequal(baseline.(name),optimized.(name)), ...
        'dsc:BenchmarkMismatch','Benchmarks differ in %s.',name);
end
assert(~baseline.profile_enabled&&~optimized.profile_enabled, ...
    'dsc:BenchmarkMismatch','Profiled timings are not comparable to unprofiled timings.');
a=load(fullfile(baseline_run,'checkpoint.mat'),'checkpoint');
b=load(fullfile(optimized_run,'checkpoint.mat'),'checkpoint');
for field={'dates','tickers','returns','mask','source_hash','priors'}
    name=field{1}; assert(isequaln(a.checkpoint.identity.(name),b.checkpoint.identity.(name)), ...
        'dsc:BenchmarkMismatch','Sampler inputs differ in %s.',name);
end
fields=fieldnames(a.checkpoint.state); state_errors=struct();
for k=1:numel(fields)
    name=fields{k}; state_errors.(name)=max_error(a.checkpoint.state.(name),b.checkpoint.state.(name));
end
assert(isequal(a.checkpoint.rng,b.checkpoint.rng),'dsc:BenchmarkMismatch','Random-generator states differ.');
for field={'iteration','r_proposals','h_proposals','r_invalid_proposals','h_invalid_proposals'}
    name=field{1}; assert(isequal(a.checkpoint.diagnostics.(name),b.checkpoint.diagnostics.(name)), ...
        'dsc:BenchmarkMismatch','Diagnostics differ in %s.',name);
end
likelihood_error=max_error(a.checkpoint.diagnostics.loglik,b.checkpoint.diagnostics.loglik);
files_a=dir(fullfile(baseline_run,'posterior_chunk_*.mat'));
files_b=dir(fullfile(optimized_run,'posterior_chunk_*.mat'));
assert(isequal({files_a.name},{files_b.name}),'dsc:BenchmarkMismatch','Retained chunks differ.');
draw_errors=struct(); verified_draws=0;
for k=1:numel(files_a)
    x=load(fullfile(files_a(k).folder,files_a(k).name),'chunk');
    y=load(fullfile(files_b(k).folder,files_b(k).name),'chunk');
    for field={'iterations','dates','tickers','pair_i','pair_j','chain_id'}
        name=field{1}; assert(isequal(x.chunk.(name),y.chunk.(name)), ...
            'dsc:BenchmarkMismatch','Chunk metadata differs in %s.',name);
    end
    verified_draws=verified_draws+numel(x.chunk.iterations);
    for field={'P_pairs','h','B','V','sig2h','sig2r'}
        name=field{1}; present=isfield(x.chunk,name);
        assert(present==isfield(y.chunk,name),'dsc:BenchmarkMismatch','Retained field %s differs.',name);
        if ~present, continue; end
        if ~isfield(draw_errors,name), draw_errors.(name)=0; end
        draw_errors.(name)=max(draw_errors.(name),max_error(x.chunk.(name),y.chunk.(name)));
    end
end
assert(verified_draws==baseline.retained_draws,'dsc:BenchmarkMismatch','Not all draws were verified.');
comparison=struct('baseline',baseline,'optimized',optimized, ...
    'seconds_saved',baseline.retained_seconds-optimized.retained_seconds, ...
    'percent_saved',100*(1-optimized.retained_seconds/baseline.retained_seconds), ...
    'speedup',baseline.retained_seconds/optimized.retained_seconds, ...
    'verified_retained_draws',verified_draws,'state_max_absolute_errors',state_errors, ...
    'retained_draw_max_absolute_errors',draw_errors,'max_absolute_loglik_error',likelihood_error, ...
    'rng_identical',true,'proposal_counts_identical',true);
save(fullfile(optimized_run,'performance_comparison.mat'),'comparison');
fid=fopen(fullfile(optimized_run,'performance_comparison.json'),'w'); assert(fid>=0);
guard=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(comparison,PrettyPrint=true));
fprintf('%d retained draws: %.3fs -> %.3fs; saved %.3fs (%.2f%%), %.2fx faster\n', ...
    verified_draws,baseline.retained_seconds,optimized.retained_seconds, ...
    comparison.seconds_saved,comparison.percent_saved,comparison.speedup);
end

function error=max_error(a,b)
assert(isequal(size(a),size(b)),'dsc:BenchmarkMismatch','Output shapes differ.');
assert(all(isfinite(a(:)))&&all(isfinite(b(:))),'dsc:BenchmarkMismatch','Output contains nonfinite values.');
error=max(abs(a(:)-b(:)),[],'omitnan'); if isempty(error), error=0; end
assert(error<=1e-6,'dsc:BenchmarkMismatch','Output difference %.12g exceeds 1e-6.',error);
end
