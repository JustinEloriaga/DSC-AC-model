function result = run_weekly_model_report(options)
%RUN_WEEKLY_MODEL_REPORT Fit constant means and save joint draws and reports.
% Example:
%   result = run_weekly_model_report(WarmupIterations=10,RetainedDraws=100);
%
% RetainedDraws is the number of draws per chain after warm-up and thinning.
% The sampler therefore performs WarmupIterations + RetainedDraws*Thin
% completed sweeps unless MaxHours stops it first.
arguments
    options.WarmupIterations (1,1) double {mustBeInteger,mustBeNonnegative} = 10
    options.RetainedDraws (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(options.RetainedDraws,2)} = 100
    options.Thin (1,1) double {mustBeInteger,mustBePositive} = 1
    options.MaxHours (1,1) double {mustBePositive,mustBeFinite} = 12
    options.ChunkSize (1,1) double = NaN
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 20260914
    options.ChainID (1,1) double {mustBeInteger,mustBePositive} = 1
    options.NumChains (1,1) double {mustBeInteger,mustBePositive} = 1
    options.ParallelChains (1,1) logical = false
    options.CorrelationBackend (1,1) string {mustBeMember(options.CorrelationBackend,["auto","mex","matlab","matlab_new","mex_new"])} = "mex_new"
    options.CorrelationThreads = "auto"
    options.MeanBackend (1,1) string {mustBeMember(options.MeanBackend,["auto","native","matlab"])} = "auto"
    options.CreatePDF (1,1) logical = true
    options.OutputDirectory (1,1) string = ""
    options.ParameterFile (1,1) string = ""
    options.SourceFile (1,1) string = ""
    options.VariableNames (1,:) string = ["NKYTR Index","SPXT Index","SX5T Index","TUKXG Index"]
    options.EstimateParameters (1,1) logical = true
    options.EstimationEndDate (1,1) string = ""
    options.KeepNewestParameterFiles (1,1) double {mustBeInteger,mustBeNonnegative} = 3
    options.BuildReport (1,1) logical = false
    options.VerifyReport (1,1) logical = false
    options.PythonExecutable (1,1) string = "python3"
end

if isnan(options.ChunkSize)
    options.ChunkSize=options.RetainedDraws;
else
    validateattributes(options.ChunkSize,{'numeric'},{'scalar','integer','positive','finite'});
end


% Compact constant means are the production model. Source CSVs and existing
% historical runs remain untouched. All new artifacts share an existing
% directory and use a unique prefix rather than creating run directories.
if options.BuildReport||options.VerifyReport
    error('DSC:LegacyReport','The historical methodology report describes dynamic means. Run with BuildReport=false to save the current model parameters and correlation PDFs.');
end
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'toolbox'),fullfile(root,'scripts'));
cfg=weekly_config();
if strlength(options.OutputDirectory)>0
    cfg.output_root=char(options.OutputDirectory);
else
    cfg.output_root=fullfile(root,'outputs');
end
assert(isfolder(cfg.output_root),'DSC:OutputDirectory','OutputDirectory must already exist.');
if strlength(options.SourceFile)>0, cfg.source_file=char(options.SourceFile); end
cfg.variable_names=options.VariableNames;
stamp=char(datetime('now','TimeZone','Europe/Madrid','Format','yyyyMMdd''T''HHmmssSSS'));
cfg.file_prefix=['constant_mu_' stamp '_'];
assert(~isempty(regexp(cfg.file_prefix,'^[A-Za-z0-9_-]+$','once')), ...
    'DSC:RunPrefix','RunPrefix may contain only letters, numbers, underscores and hyphens.');
if ~isempty(dir(fullfile(cfg.output_root,[cfg.file_prefix '*'])))
    error('DSC:ExistingRunPrefix','Generated run prefix already exists; retry the run.');
end
cfg.burnin=options.WarmupIterations; cfg.thin=options.Thin;
cfg.max_iterations=cfg.burnin+options.RetainedDraws*cfg.thin;
cfg.max_seconds=options.MaxHours*3600; cfg.chunk_size=options.ChunkSize;
cfg.seed=options.Seed; cfg.chain_id=options.ChainID;
cfg.num_chains=options.NumChains; cfg.parallel_chains=options.ParallelChains;
cfg.resume=false;
cfg.mean_backend=char(options.MeanBackend);
cfg.correlation_backend=char(options.CorrelationBackend);
cfg.correlation_warm_start=true;
cfg.correlation_threads=options.CorrelationThreads;
cfg.estimate_parameters=options.EstimateParameters;
cfg.estimation_end_date=char(options.EstimationEndDate);
cfg.create_pdf=options.CreatePDF;
cfg.parameter_dir=cfg.output_root;
if ~isempty(dir(fullfile(cfg.parameter_dir, ...
        [cfg.file_prefix 'parameters.mat'])))
    error('DSC:ExistingRunPrefix','Generated run prefix already exists; retry the run.');
end
cfg.parameter_file=char(options.ParameterFile);
if cfg.parallel_chains
    existing=gcp('nocreate');
    assert(~isempty(existing),'DSC:ParallelPool', ...
        'ParallelChains requires an existing process pool; starting one creates job directories. Use sequential chains or start the pool explicitly.');
end
panel=prepare_weekly_panel(cfg);
[priors,inference]=dsc_prepare_inference(panel,cfg);
% Select the native thread count once from this workload and CPU budget.
% Coordinated batches receive a fixed numeric count for exact continuation.
cfg.correlation_threads=dsc_resolve_correlation_threads( ...
    cfg.correlation_threads,cfg.num_chains,cfg.parallel_chains, ...
    inference.smoothing_panel,priors,cfg.correlation_backend);
manifest=struct('model','joint_dsc_constant_mean','mean_model','constant', ...
    'equation','y_t = mu + D_t P_t^(1/2) z_t; z_t iid N(0,I)', ...
    'innovation_variances','sig2h governs h (delta) random walks; sig2r governs r (gamma) random walks', ...
    'source_file',panel.source_file,'source_hash',panel.source_hash, ...
    'tickers',panel.tickers,'config',cfg,'matlab_version',version, ...
    'model_start',char(string(inference.smoothing_panel.dates(1),'yyyy-MM-dd')), ...
    'model_end',char(string(inference.smoothing_panel.dates(end),'yyyy-MM-dd')));
result=dsc_run_inference(panel,priors,cfg,cfg.output_root,inference);
if result.saved_draws<2
    error('DSC:TooFewRetainedDraws','Only %d retained draws; completed checkpoints were preserved.',result.saved_draws);
end
manifest.inference=result.inference;
result.figure_paths=strings(0,1);
result.report_pdf=''; result.figure_pdf='';
result.publication_checked=false; result.run_prefix=cfg.file_prefix;
result.model='joint_dsc_constant_mean'; result.mean_model='constant';
result.delta_definition='delta is the log-volatility state h; sig2h is its random-walk innovation variance';
result.gamma_definition='gamma is the matrix-log correlation state r; sig2r is its random-walk innovation variance';
result=dsc_finalize_outputs(result,cfg,manifest,inference.smoothing_panel);
end
