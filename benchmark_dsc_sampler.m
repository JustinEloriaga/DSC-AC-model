function report=benchmark_dsc_sampler(checkpoint_file,options)
% Time retained draws on an identical saved panel, prior, seed and warm-up.
arguments
    checkpoint_file (1,1) string
    options.SourceRoot (1,1) string = string(fileparts(mfilename('fullpath')))
    options.OutputRoot (1,1) string = ""
    options.Label (1,1) string = "sampler"
    options.RetainedDraws (1,1) double {mustBeInteger,mustBePositive} = 100
    options.WarmupIterations (1,1) double {mustBeInteger,mustBeNonnegative} = 10
    options.Profile (1,1) logical = false
    options.CorrelationBackend (1,1) string {mustBeMember(options.CorrelationBackend,["auto","mex","matlab"])} = "auto"
    options.CorrelationThreads (1,1) double {mustBeInteger,mustBePositive} = 4
end
root=fileparts(mfilename('fullpath'));
if strlength(options.OutputRoot)==0
    options.OutputRoot=fullfile(root,'outputs','weekly_research','performance','sampler_optimization');
end
previous_path=path; restore=onCleanup(@()path(previous_path));
addpath(fullfile(options.SourceRoot,'toolbox'));
clear dsc_sample dsc_sampler_signature dsc_corr_batch_mex dsc_mean_mex
clear dsc_corr_matrix dsc_random_walk_draw dsc_volatility_likelihood
loaded=load(checkpoint_file,'checkpoint'); identity=loaded.checkpoint.identity;
panel=struct('returns',identity.returns,'observation_mask',identity.mask, ...
    'dates',identity.dates,'tickers',identity.tickers,'source_hash',identity.source_hash);
[panel.pair_i,panel.pair_j]=find(tril(true(numel(panel.tickers)),-1));
cfg=identity.cfg; cfg.burnin=options.WarmupIterations; cfg.thin=1;
cfg.chunk_size=10; cfg.resume=false; cfg.max_seconds=12*3600;
cfg.max_iterations=cfg.burnin; cfg.correlation_backend=char(options.CorrelationBackend);
cfg.correlation_threads=options.CorrelationThreads;
stamp=char(datetime('now','Format','yyyyMMdd-HHmmss-SSS'));
run_dir=fullfile(options.OutputRoot,options.Label+"-"+stamp); mkdir(run_dir);
fprintf('BENCHMARK %s: %d dates, %d variables, %d warm-up + %d retained draws\n', ...
    options.Label,size(panel.returns,1),size(panel.returns,2),cfg.burnin,options.RetainedDraws);
warmup_clock=tic;
if cfg.burnin>0
    warmup=dsc_sample(panel,identity.priors,cfg,run_dir);
else
    warmup=struct('completed_iterations',0,'stage_seconds',struct());
end
warmup_seconds=toc(warmup_clock);
assert(warmup.completed_iterations==cfg.burnin,'dsc:BenchmarkIncomplete','Warm-up did not finish.');
cfg.resume=cfg.burnin>0; cfg.max_iterations=cfg.burnin+options.RetainedDraws;
if options.Profile, profile clear; profile on; end
draw_clock=tic;
result=dsc_sample(panel,identity.priors,cfg,run_dir);
retained_seconds=toc(draw_clock);
if options.Profile
    profile off; sampler_profile=profile('info');
    save(fullfile(run_dir,'sampler_profile.mat'),'sampler_profile');
    functions=sampler_profile.FunctionTable; [~,order]=sort([functions.TotalTime],'descend');
    hotspots=arrayfun(@(f)struct('function',f.FunctionName,'seconds',f.TotalTime, ...
        'calls',f.NumCalls),functions(order(1:min(20,numel(order)))));
else
    hotspots=[];
end
assert(result.saved_draws==options.RetainedDraws,'dsc:BenchmarkIncomplete','Draw target did not finish.');
stages=fieldnames(result.stage_seconds); retained_stage_seconds=struct();
for k=1:numel(stages)
    name=stages{k}; previous=0;
    if isfield(warmup.stage_seconds,name), previous=warmup.stage_seconds.(name); end
    retained_stage_seconds.(name)=result.stage_seconds.(name)-previous;
end
report=struct('label',options.Label,'input_checkpoint',checkpoint_file, ...
    'source_root',options.SourceRoot,'run_dir',run_dir,'source_hash',panel.source_hash, ...
    'tickers',panel.tickers,'dates',numel(panel.dates),'seed',result.seed,'chain_id',cfg.chain_id, ...
    'warmup_iterations',cfg.burnin,'retained_draws',result.saved_draws, ...
    'correlation_threads',cfg.correlation_threads,'profile_enabled',options.Profile, ...
    'warmup_seconds',warmup_seconds,'retained_seconds',retained_seconds, ...
    'seconds_per_retained_draw',retained_seconds/result.saved_draws, ...
    'retained_stage_seconds',retained_stage_seconds, ...
    'retained_sweep_seconds',result.sweep_seconds(cfg.burnin+1:end), ...
    'implementation',result.implementation,'hotspots',hotspots, ...
    'scope','Sampler wall time including retained chunks, checkpoints and final sampler files; excludes warm-up, data analysis and figures.');
save(fullfile(run_dir,'benchmark_report.mat'),'report');
fid=fopen(fullfile(run_dir,'benchmark_report.json'),'w'); assert(fid>=0);
guard=onCleanup(@()fclose(fid)); fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('BENCHMARK %s COMPLETE: %d retained draws in %.3fs (%.3fs/draw)\n', ...
    options.Label,result.saved_draws,retained_seconds,report.seconds_per_retained_draw);
end
