function result = run_weekly_model_report(options)
%RUN_WEEKLY_MODEL_REPORT Run the weekly DSC model and build its PDF report.
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
    options.ChunkSize (1,1) double {mustBeInteger,mustBePositive} = 10
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 20260914
    options.ChainID (1,1) double {mustBeInteger,mustBePositive} = 1
    options.NumChains (1,1) double {mustBeInteger,mustBePositive} = 1
    options.ParallelChains (1,1) logical = false
    options.ConvergenceMode (1,1) string {mustBeMember(options.ConvergenceMode,["off","report"])} = "report"
    options.AutoExtend (1,1) logical = false
    options.CheckEvery (1,1) double {mustBeInteger,mustBePositive} = 250
    options.MaxRetainedDraws (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(options.MaxRetainedDraws,2)} = 2000
    options.MinDiagnosticDraws (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(options.MinDiagnosticDraws,6)} = 100
    options.RhatThreshold (1,1) double {mustBeGreaterThan(options.RhatThreshold,1),mustBeFinite} = 1.01
    options.MinESS (1,1) double {mustBePositive,mustBeFinite} = 400
    options.MaxMCSERatio (1,1) double {mustBePositive,mustBeFinite} = 0.05
    options.CorrelationBackend (1,1) string {mustBeMember(options.CorrelationBackend,["auto","mex","matlab"])} = "auto"
    options.CorrelationThreads (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(options.CorrelationThreads,1),mustBeLessThanOrEqual(options.CorrelationThreads,64)} = 4
    options.SourceFile (1,1) string = ""
    options.OutputRoot (1,1) string = ""
    options.VariableNames (1,:) string = "all"
    options.EstimateParameters (1,1) logical = true
    options.EstimationEndDate (1,1) string = ""
    options.ParameterDirectory (1,1) string = ""
    options.ParameterFile (1,1) string = ""
    options.ParameterSmoothing (1,1) string {mustBeMember(options.ParameterSmoothing,["draws","mean"])} = "mean"
    options.RunTests (1,1) logical = true
    options.BuildReport (1,1) logical = true
    options.VerifyReport (1,1) logical = true
    options.PythonExecutable (1,1) string = "python3"
end

if options.VerifyReport && ~options.BuildReport
    error('DSC:PipelineOptions','VerifyReport requires BuildReport=true.');
end
root = fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'toolbox'),fullfile(root,'core'),fullfile(root,'scripts'));

cfg = weekly_config();
if strlength(options.SourceFile)>0, cfg.source_file=char(options.SourceFile); end
if strlength(options.OutputRoot)>0, cfg.output_root=char(options.OutputRoot); end
cfg.variable_names = options.VariableNames;
if ~exist(cfg.source_file,'file')
    error('DSC:SourceFile','Current-data source file does not exist: %s',cfg.source_file);
end
cfg.burnin = options.WarmupIterations;
cfg.thin = options.Thin;
cfg.max_iterations = options.WarmupIterations + options.RetainedDraws*options.Thin;
cfg.max_seconds = options.MaxHours*3600;
cfg.chunk_size = options.ChunkSize;
cfg.seed = options.Seed;
cfg.chain_id = options.ChainID;
cfg.num_chains = options.NumChains;
cfg.parallel_chains = options.ParallelChains;
cfg.convergence_mode = char(options.ConvergenceMode);
cfg.auto_extend = options.AutoExtend;
cfg.check_every = options.CheckEvery;
cfg.max_retained_draws = options.MaxRetainedDraws;
cfg.min_diagnostic_draws = options.MinDiagnosticDraws;
cfg.rhat_threshold = options.RhatThreshold;
cfg.min_ess = options.MinESS;
cfg.max_mcse_ratio = options.MaxMCSERatio;
cfg.resume = false;
cfg.correlation_backend = char(options.CorrelationBackend);
cfg.correlation_threads = options.CorrelationThreads;
cfg.estimate_parameters = options.EstimateParameters;
cfg.estimation_end_date = char(options.EstimationEndDate);
cfg.parameter_dir = char(options.ParameterDirectory);
cfg.parameter_file = char(options.ParameterFile);
cfg.parameter_smoothing = char(options.ParameterSmoothing);

% Reject unsupported execution settings before data preparation or sampling.
if cfg.auto_extend&&(strcmp(cfg.convergence_mode,'off')||cfg.max_retained_draws<options.RetainedDraws)
    error('dsc:ConvergenceConfig','AutoExtend needs convergence checking and MaxRetainedDraws >= RetainedDraws.');
end
if cfg.parallel_chains&&~dsc_parallel_available()
    error('dsc:ParallelUnavailable','ParallelChains=true requires Parallel Computing Toolbox and an available license. Use ParallelChains=false for sequential chains.');
end

if options.RunTests
    run_weekly_acceptance(cfg);
end
analysis = run_weekly_research('analysis',cfg);

stamp = char(datetime('now','Format','yyyyMMdd-HHmmss-SSS'));
run_name = sprintf('%s-chain%d-w%d-r%d-thin%d',stamp,cfg.chain_id,cfg.burnin,options.RetainedDraws,cfg.thin);
run_dir = fullfile(cfg.output_root,'runs',run_name);
if ~exist(run_dir,'dir'), mkdir(run_dir); end
run_manifest = analysis.manifest;
run_manifest.run_kind = 'configurable sampler and report pipeline';
run_manifest.requested_warmup_iterations = options.WarmupIterations;
run_manifest.requested_retained_draws = options.RetainedDraws;
run_manifest.requested_total_sweeps = cfg.max_iterations;
run_manifest.requested_num_chains = cfg.num_chains;
run_manifest.parallel_chains = cfg.parallel_chains;
run_manifest.convergence_mode = cfg.convergence_mode;
run_manifest.requested_max_hours = options.MaxHours;
run_manifest.convergence_established = false;
write_json_local(fullfile(run_dir,'run_manifest.json'),run_manifest);

result = dsc_run_inference(analysis.panel,analysis.priors,cfg,run_dir,analysis.inference);
if result.saved_draws < 2
    error('DSC:TooFewRetainedDraws', ...
        ['The run was preserved at %s, but it retained only %d draw(s). ' ...
         'Increase MaxHours or reduce WarmupIterations.'],run_dir,result.saved_draws);
end
if any(result.saved_draws_per_chain < options.RetainedDraws)
    warning('DSC:IncompleteDrawTarget', ...
        'Some chains retained fewer than %d requested draws each; %d draws were retained in total.', ...
        options.RetainedDraws,result.saved_draws);
end

result.figure_paths = plot_bayes_correlation_paths(run_dir);
result.report_pdf = '';
result.publication_checked = false;
if options.BuildReport
    builder = fullfile(root,'scripts','build_publication.py');
    summary_path = fullfile(run_dir,'pilot_summary.json');
    command = sprintf('%s %s --results %s --pilot %s --posterior-run %s --report-only', ...
        shell_arg(options.PythonExecutable),shell_arg(builder),shell_arg(cfg.output_root), ...
        shell_arg(summary_path),shell_arg(run_dir));
    [status,output] = system(command,'-echo');
    if status~=0
        error('DSC:PublicationBuild','Report builder failed with status %d:\n%s',status,output);
    end
    result.report_pdf = fullfile(root,'output','pdf','weekly_correlation_report.pdf');
    if options.VerifyReport
        checker = fullfile(root,'scripts','check_publication.py');
        command = sprintf('%s %s --report-only',shell_arg(options.PythonExecutable),shell_arg(checker));
        [status,output] = system(command,'-echo');
        if status~=0
            error('DSC:PublicationCheck','Report verification failed with status %d:\n%s',status,output);
        end
        result.publication_checked = true;
    end
end

latest = struct('run_dir',run_dir,'summary',fullfile(run_dir,'pilot_summary.json'), ...
    'figures',fullfile(run_dir,'figures'),'report_pdf',result.report_pdf, ...
    'requested_retained_draws',options.RetainedDraws,'actual_retained_draws',result.saved_draws, ...
    'warmup_iterations',cfg.burnin,'thin',cfg.thin,'completed_iterations',result.completed_iterations, ...
    'status',result.status,'created_at',char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss')));
latest.inference=result.inference;
latest.num_chains=cfg.num_chains;
latest.actual_retained_draws_per_chain=result.saved_draws_per_chain;
latest.convergence=result.convergence;
write_json_local(fullfile(cfg.output_root,'latest_model_report.json'),latest);
end

function value = shell_arg(value)
value = char(value);
quote = char(39);
value = strrep(value,quote,[quote '"' quote '"' quote]);
value = [quote value quote];
end

function write_json_local(path,value)
fid=fopen(path,'w');
assert(fid>=0,'DSC:WriteFailed','Cannot write %s.',path);
guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
